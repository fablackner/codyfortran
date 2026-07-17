! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysGauge_Ylm_Velocity_Fedvr) S_SysGauge_Ylm_Velocity_Fedvr

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysGauge_Ylm_Velocity_Fedvr_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysGauge
    use M_SysGauge_Ylm

    call Say_Fabricate("sysGauge.ylm.velocity.fedvr")

    !------------------------------------
    ! bind procedure pointers
    !------------------------------------

    SysGauge_Setup => Setup
    SysGauge_Ylm_ApplyRadialFirstDerivative => ApplyRadialFirstDerivative

    !------------------------------------
    ! validate grid requirements
    !------------------------------------

    if (.not. Json_GetExistence("grid.ylm.fedvr")) then
      error stop "grid.ylm.fedvr is required for sysGauge.ylm.velocity.fedvr"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Capture the radial coordinates entering the ±k/r coupling terms.
  subroutine Setup
    use M_Utils_Say
    use M_Grid_Ylm
    use M_SysGauge_Ylm

    call Say_Setup("sysGauge.ylm.velocity.fedvr")

    SysGauge_Ylm_radialCoordinates = cmplx(Grid_Ylm_radialPoints(:), 0.0_R64, R64)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Compute d fLm / dr using the FEDVR derivative context of the grid.
  subroutine ApplyRadialFirstDerivative(dfLm, fLm)
    use M_Utils_Fedvr
    use M_Utils_DerivativeFedvr
    use M_Grid_Ylm_Fedvr, fedvrCtx => Grid_Ylm_Fedvr_fedvrCtx, derivativeCtx => Grid_Ylm_Fedvr_derivativeCtx

    complex(R64), intent(out) :: dfLm(:)
    complex(R64), intent(in)  :: fLm(:)

    call DerivativeFedvr_Do1stDerivative(dfLm, fLm, derivativeCtx, fedvrCtx)

  end subroutine

end submodule
