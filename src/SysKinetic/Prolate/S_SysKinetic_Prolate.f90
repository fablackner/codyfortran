! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file S_SysKinetic_Prolate.f90
!> @brief Implementation submodule for prolate kinetic fabrication and operator.
submodule(M_SysKinetic_Prolate) S_SysKinetic_Prolate

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Bind the global kinetic operator and dispatch to the channel backend.
  module subroutine SysKinetic_Prolate_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysKinetic
    use M_SysKinetic_Prolate_Laplacian

    call Say_Fabricate("sysKinetic.prolate")

    !------------------------------------
    ! bind the global kinetic operator
    !------------------------------------

    SysKinetic_MultiplyWithKineticOp => MultiplyWithKineticOp

    !------------------------------------
    ! dispatch to channel backend
    !------------------------------------

    if (Json_GetExistence("sysKinetic.prolate.laplacian")) then
      call SysKinetic_Prolate_Laplacian_Fabricate

    else
      error stop "sysKinetic.prolate is missing one of: laplacian"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Apply T psi by looping over all azimuthal channels.
  subroutine MultiplyWithKineticOp(dOrb, orb, time, bt_)
    use M_Grid
    use M_Grid_Prolate

    complex(R64), intent(out), contiguous :: dOrb(:)
    complex(R64), intent(in), contiguous :: orb(:)
    real(R64), intent(in)                :: time
    integer(I32), intent(in), optional   :: bt_

    integer(I32) :: m
    complex(R64), allocatable :: orbM(:), dOrbM(:)

    dOrb = 0.0_R64

    allocate (orbM(Grid_Prolate_nSpatial), dOrbM(Grid_Prolate_nSpatial))

    ! Process each azimuthal channel separately
    do m = -Grid_Prolate_mmax, Grid_Prolate_mmax
      call Grid_Prolate_GetMComponent(orbM, m, orb)

      ! Skip negligible channels for efficiency
      if (all(abs(orbM) < 1.0e-14_R64)) cycle

      call SysKinetic_Prolate_MultiplyWithChannelKineticOp(dOrbM, orbM, m, time, bt_)

      call Grid_Prolate_SetMComponent(dOrb, m, dOrbM)
    end do

    deallocate (orbM, dOrbM)

  end subroutine

end submodule
