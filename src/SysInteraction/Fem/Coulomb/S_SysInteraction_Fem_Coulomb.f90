! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Implementation submodule for the Coulomb interaction on grid.fem.
submodule(M_SysInteraction_Fem_Coulomb) S_SysInteraction_Fem_Coulomb

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysInteraction_Fem_Coulomb_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysInteraction

    call Say_Fabricate("sysInteraction.fem.coulomb")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysInteraction_Fem_Coulomb_strength = Json_Get("sysInteraction.fem.coulomb.strength", 1.0_R64)

    SysInteraction_timeIndependentQ = .true.
    SysInteraction_bodyTypeIndependentQ = .true.

    ! Real kernel and real weights: the swapped source pair gives the
    ! conjugated potential
    SysInteraction_conjSymmetricQ = .true.

    SysInteraction_FillInteractionPotential => FillInteractionPotential

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief V = strength * (Galerkin Coulomb potential of the weighted source
  !> src = rho w); time, bt1_ and bt2_ are ignored (Method passes -1 for
  !> body-type-independent interactions).
  subroutine FillInteractionPotential(interactionPotential, src, time, bt1_, bt2_)
    !$ use omp_lib, only: omp_in_parallel
    use M_Utils_UnusedVariables
    use M_Grid
    use M_Grid_Fem

    complex(R64), intent(out), allocatable :: interactionPotential(:)
    complex(R64), intent(in), contiguous   :: src(:)
    real(R64), intent(in)                  :: time
    integer(I32), intent(in), optional     :: bt1_
    integer(I32), intent(in), optional     :: bt2_

    integer(I32) :: i

    if (.false.) call UnusedVariables_Mark(time, bt1_, bt2_)

    allocate (interactionPotential(Grid_nPoints))

    call Grid_Fem_SolvePoisson(interactionPotential, src)

    if (abs(SysInteraction_Fem_Coulomb_strength - 1.0_R64) > 0.0_R64) then
      !$omp parallel do default(shared) private(i) schedule(static) if (.not. omp_in_parallel())
      do i = 1, Grid_nPoints
        interactionPotential(i) = SysInteraction_Fem_Coulomb_strength * interactionPotential(i)
      end do
      !$omp end parallel do
    end if

  end subroutine

end submodule
