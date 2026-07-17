! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> FEDVR-ECS radial first derivative for the Ylm velocity-gauge coupling.
!>
!> Binds the radial building blocks of M_SysGauge_Ylm to the FEDVR-ECS
!> derivative context: derivatives are taken along the complex contour z(r)
!> and the ±k/r coupling terms use the complex contour points, matching the
!> analytic continuation used by the ECS kinetic operator (c-product metric).
module M_SysGauge_Ylm_Velocity_FedvrEcs
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Bind the FEDVR-ECS radial derivative and validate grid requirements.
    module subroutine SysGauge_Ylm_Velocity_FedvrEcs_Fabricate
    end subroutine
  end interface

end module
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
