! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> FEM-grid fabrication submodule for external potentials.
!>
!> Dispatches to the potential model and provides the point-wise
!> multiplication shared by all nodal FEM potentials.
submodule(M_SysPotential_Fem) S_SysPotential_Fem

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Check the grid, bind the point-wise multiplication and wire the model.
  !>
  !> Available models:
  !> - `coulomb`: static nuclear attraction of the grid.fem centers (Galerkin
  !>   point nuclei or erf-smeared nuclei)
  module subroutine SysPotential_Fem_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysPotential, only: SysPotential_MultiplyWithExternalPotential
    use M_SysPotential_Fem_Coulomb

    call Say_Fabricate("sysPotential.fem")

    !------------------------------------
    ! validate grid requirements
    !------------------------------------

    if (.not. Json_GetExistence("grid.fem")) then
      error stop "grid.fem is required for sysPotential.fem"
    end if

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysPotential_MultiplyWithExternalPotential => MultiplyWithExternalPotential

    !------------------------------------
    ! branch by potential model
    !------------------------------------

    if (Json_GetExistence("sysPotential.fem.coulomb")) then
      call SysPotential_Fem_Coulomb_Fabricate

    else
      error stop "sysPotential.fem is missing one of: coulomb"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Apply the nodal potential: dOrb(i) = V(i) orb(i). Reentrant: serial
  !> inside the OpenMP loops of Method, parallel over the nodes when called
  !> outside a parallel region (e.g. by the Cayley exact part of Method).
  subroutine MultiplyWithExternalPotential(dOrb, externalPotential, orb)
    !$ use omp_lib, only: omp_in_parallel

    complex(R64), intent(out), contiguous :: dOrb(:)
    complex(R64), intent(in), contiguous  :: externalPotential(:)
    complex(R64), intent(in), contiguous  :: orb(:)

    integer(I32) :: i

    !$omp parallel do default(shared) private(i) schedule(static) if (.not. omp_in_parallel())
    do i = 1, size(dOrb)
      dOrb(i) = externalPotential(i) * orb(i)
    end do
    !$omp end parallel do

  end subroutine

end submodule
