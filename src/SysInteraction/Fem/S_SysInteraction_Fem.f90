! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Implementation submodule for FEM-grid interactions.
!>
!> @details Nodal source, multiplication and conjugation; dispatch to the
!> interaction model.
submodule(M_SysInteraction_Fem) S_SysInteraction_Fem

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Fabricate the FEM-grid interaction back-end.
  module subroutine SysInteraction_Fem_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysInteraction
    use M_SysInteraction_Fem_Coulomb

    call Say_Fabricate("sysInteraction.fem")

    !------------------------------------
    ! validate grid requirements
    !------------------------------------

    if (.not. Json_GetExistence("grid.fem")) then
      error stop "grid.fem is required for sysInteraction.fem"
    end if

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysInteraction_FillInteractionSrc => FillInteractionSrc
    SysInteraction_MultiplyWithInteractionPotential => MultiplyWithInteractionPotential
    SysInteraction_ConjugateInteractionPotential => ConjugateInteractionPotential

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("sysInteraction.fem.coulomb")) then
      call SysInteraction_Fem_Coulomb_Fabricate

    else
      error stop "sysInteraction.fem is missing one of: coulomb"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Weighted pair density src(i) = conjg(orbConjg(i)) orb(i) w(i), the
  !> Galerkin load of rho = conjg(orbConjg) orb.
  subroutine FillInteractionSrc(src, orbConjg, orb)
    !$ use omp_lib, only: omp_in_parallel
    use M_Grid
    use M_Grid_Fem

    complex(R64), intent(out), allocatable :: src(:)
    complex(R64), intent(in), contiguous   :: orbConjg(:)
    complex(R64), intent(in), contiguous   :: orb(:)

    integer(I32) :: i

    if (size(orbConjg) .ne. Grid_nPoints .or. size(orb) .ne. Grid_nPoints) then
      error stop "sysInteraction.fem: orbitals must have Grid_nPoints entries"
    end if

    allocate (src(Grid_nPoints))

    !$omp parallel do default(shared) private(i) schedule(static) if (.not. omp_in_parallel())
    do i = 1, Grid_nPoints
      src(i) = conjg(orbConjg(i)) * orb(i) * Grid_Fem_weights(i)
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Apply the nodal interaction potential: dOrb(i) = V(i) orb(i)
  !> (reentrant: parallel when called at top level, e.g. by the frozen
  !> mean-field reference of the Cayley step, serial inside the OpenMP loops
  !> of Method).
  subroutine MultiplyWithInteractionPotential(dOrb, interactionPotential, orb)
    !$ use omp_lib, only: omp_in_parallel

    complex(R64), intent(out), contiguous :: dOrb(:)
    complex(R64), intent(in), contiguous  :: interactionPotential(:)
    complex(R64), intent(in), contiguous  :: orb(:)

    integer(I32) :: i

    !$omp parallel do default(shared) private(i) schedule(static) if (.not. omp_in_parallel())
    do i = 1, size(dOrb)
      dOrb(i) = interactionPotential(i) * orb(i)
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Potential of the swapped pair: conjg(V) in nodal storage.
  subroutine ConjugateInteractionPotential(potOut, potIn)

    complex(R64), intent(out), contiguous :: potOut(:)
    complex(R64), intent(in), contiguous  :: potIn(:)

    potOut(:) = conjg(potIn(:))

  end subroutine

end submodule
