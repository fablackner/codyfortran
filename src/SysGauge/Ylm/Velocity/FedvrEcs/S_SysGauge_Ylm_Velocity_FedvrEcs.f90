! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysGauge_Ylm_Velocity_FedvrEcs) S_SysGauge_Ylm_Velocity_FedvrEcs

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysGauge_Ylm_Velocity_FedvrEcs_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysGauge
    use M_SysGauge_Ylm

    call Say_Fabricate("sysGauge.ylm.velocity.fedvrEcs")

    !------------------------------------
    ! bind procedure pointers
    !------------------------------------

    SysGauge_Setup => Setup
    SysGauge_Ylm_ApplyRadialFirstDerivative => ApplyRadialFirstDerivative

    !------------------------------------
    ! validate grid requirements
    !------------------------------------

    if (.not. Json_GetExistence("grid.ylm.fedvrEcs")) then
      error stop "grid.ylm.fedvrEcs is required for sysGauge.ylm.velocity.fedvrEcs"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Capture the contour points entering the ±k/z coupling terms.
  subroutine Setup
    use M_Utils_Say
    use M_Grid_Ylm_FedvrEcs, contourPoints => Grid_Ylm_FedvrEcs_contourPoints
    use M_SysGauge_Ylm

    call Say_Setup("sysGauge.ylm.velocity.fedvrEcs")

    SysGauge_Ylm_radialCoordinates = contourPoints(:)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Compute d fLm / dz along the ECS contour using the grid derivative context.
  subroutine ApplyRadialFirstDerivative(dfLm, fLm)
    use M_Utils_FedvrEcs
    use M_Utils_DerivativeFedvrEcs
    use M_Grid_Ylm_FedvrEcs, fedvrEcsCtx => Grid_Ylm_FedvrEcs_fedvrEcsCtx, &
      derivativeCtx => Grid_Ylm_FedvrEcs_derivativeCtx

    complex(R64), intent(out) :: dfLm(:)
    complex(R64), intent(in)  :: fLm(:)

    call DerivativeFedvrEcs_Do1stDerivative(dfLm, fLm, derivativeCtx, fedvrEcsCtx)

  end subroutine

end submodule
