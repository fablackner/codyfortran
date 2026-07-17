! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Linearly polarized laser pulse in a Ylm representation (velocity gauge,
!> dipole approximation).
!>
!> The coupling term for an electron (charge -1, a.u.) is
!>   H_A = A(t) p_z / m + A(t)² / (2 m)
!> with the vector potential along z given by a sin²-envelope pulse
!>   A(t) = -(E0/omega) * sin²(pi (t - t0)/T) * sin(omega (t - t0 - T/2) + cep)
!> for t0 <= t <= t0 + T (zero outside), where T = 2 pi nCycles / omega. The
!> carrier is centered on the envelope maximum, matching the CoulombLaser
!> convention. The electric field E(t) = -dA/dt then equals the CoulombLaser
!> length-gauge field E0 sin² cos plus the envelope-derivative correction of
!> relative size 1/nCycles; A(t) vanishes identically at both pulse edges.
!>
!> The spatially uniform A²/(2m) term only contributes a global time-dependent
!> phase; it is included so that the operator is exactly (p + A)²/2m − p²/2m.
module M_SysGauge_Ylm_Velocity
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Set pulse parameters from config and wire the implementation.
    module subroutine SysGauge_Ylm_Velocity_Fabricate
    end subroutine

    !> Vector potential amplitude A(t) of the pulse along z (a.u.).
    module function SysGauge_Ylm_Velocity_VectorPotentialAmplitude(time) result(res)
      real(R64), intent(in) :: time
      real(R64) :: res
    end function
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Peak electric field strength E0 of the pulse (a.u.); A0 = E0 / omega.
  real(R64) :: SysGauge_Ylm_Velocity_fieldStrength
  !> Carrier angular frequency omega of the pulse (a.u.).
  real(R64) :: SysGauge_Ylm_Velocity_omega
  !> Number of optical cycles under the sin² envelope (T = 2 pi nCycles / omega).
  real(R64) :: SysGauge_Ylm_Velocity_nCycles
  !> Carrier-envelope phase (rad).
  real(R64) :: SysGauge_Ylm_Velocity_cep
  !> Start time t0 of the pulse (a.u.).
  real(R64) :: SysGauge_Ylm_Velocity_tStart
  !> Mass per body type entering A·p/m and A²/(2m).
  real(R64), allocatable :: SysGauge_Ylm_Velocity_bodyMass(:)

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
