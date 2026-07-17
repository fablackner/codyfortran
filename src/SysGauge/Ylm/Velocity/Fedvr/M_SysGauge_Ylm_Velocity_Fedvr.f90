! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> FEDVR radial first derivative for the Ylm velocity-gauge coupling.
!>
!> Binds the radial building blocks of M_SysGauge_Ylm to the FEDVR derivative
!> context owned by the Ylm FEDVR grid. Applied to g = r·f (see the gauge
!> model), the weak-form first derivative is anti-Hermitian in the
!> FEDVR-weighted metric up to the boundary term g_bra(rmax)·g_ket(rmax):
!> the radial grid includes the rmax endpoint, so A·p is Hermitian exactly
!> for wavefunctions that vanish at rmax and to |ψ(rmax)|² accuracy otherwise.
module M_SysGauge_Ylm_Velocity_Fedvr
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Bind the FEDVR radial derivative and validate grid requirements.
    module subroutine SysGauge_Ylm_Velocity_Fedvr_Fabricate
    end subroutine
  end interface

end module
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
