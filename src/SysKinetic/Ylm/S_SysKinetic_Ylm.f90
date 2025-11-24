! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysKinetic_Ylm) S_SysKinetic_Ylm

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysKinetic_Ylm_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysKinetic
    use M_SysKinetic_Ylm_Laplacian

    call Say_Fabricate("sysKinetic.ylm")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysKinetic_MultiplyWithKineticOp => MultiplyWithKineticOp

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("sysKinetic.ylm.laplacian")) then
      call SysKinetic_Ylm_Laplacian_Fabricate

    else
      error stop "sysKinetic.ylm is missing one of: laplacian"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine MultiplyWithKineticOp(dOrb, orb, time, bt_)
    use M_Grid
    use M_Grid_Ylm
    use M_Utils_Constants

    complex(R64), intent(out), contiguous :: dOrb(:)
    complex(R64), intent(in), contiguous :: orb(:)
    real(R64), intent(in)                :: time
    integer(I32), intent(in), optional   :: bt_

    integer(I32) :: l, m, lmax
    integer(I32) :: nRad
    complex(R64), allocatable :: orbLm(:), dOrbLm(:)

    nRad = Grid_Ylm_nRadial
    lmax = Grid_Ylm_lmax

    dOrb = 0.0_R64

    ! Allocate arrays for radial functions
    allocate (orbLm(nRad), dOrbLm(nRad))

    ! Process each (l,m) pair separately
    do l = 0, lmax
      do m = -l, l
        ! Extract all radial points for this (l,m) pair
        call Grid_Ylm_GetLmComponent(orbLm, l, m, orb)
        if (all(abs(orbLm) < 1.0e-14_R64)) cycle

        ! MultiplyWith radial part of the Laplacian
        call SysKinetic_Ylm_MultiplyWithRadialKineticOp(dOrbLm, orbLm, l, m, time, bt_)

        ! Add the results back to the output array
        call Grid_Ylm_SetLmComponent(dOrb, l, m, dOrbLm)
      end do
    end do

    deallocate (orbLm, dOrbLm)

  end subroutine

end submodule
