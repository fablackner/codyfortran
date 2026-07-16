! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Central Coulomb potential plus a linearly polarized laser pulse in a Ylm
!> representation (length gauge, dipole approximation).
!>
!> The total one-body potential is
!>   V(r, Omega, t) = -Z/r + E(t) * z
!> with the electric field along z given by a sin^2-envelope pulse
!>   E(t) = E0 * sin^2(pi * (t - t0) / T) * cos(omega * (t - t0 - T/2) + cep)
!> for t0 <= t <= t0 + T (zero outside), where T = 2*pi*nCycles/omega. The
!> carrier is centered on the envelope maximum so that cep = 0 gives a
!> cos-like pulse peaking at the envelope peak.
!>
!> In the Ylm expansion only two components are nonzero:
!>   V_00(r) = -Z * sqrt(4 pi) / r
!>   V_10(r) = E(t) * r * sqrt(4 pi / 3)      (since z = r sqrt(4 pi/3) Y_10)
module M_SysPotential_Ylm_CoulombLaser
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> set pulse parameters from config and wire the implementation.
    module subroutine SysPotential_Ylm_CoulombLaser_Fabricate
    end subroutine

    !> Electric field amplitude E(t) of the pulse along z (a.u.).
    module function SysPotential_Ylm_CoulombLaser_FieldAmplitude(time) result(res)
      real(R64), intent(in) :: time
      real(R64) :: res
    end function
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Total charge Z of the central Coulomb potential (V(r) = -Z / r in a.u.).
  real(R64) :: SysPotential_Ylm_CoulombLaser_charge
  !> Peak electric field strength E0 of the pulse (a.u.).
  real(R64) :: SysPotential_Ylm_CoulombLaser_fieldStrength
  !> Carrier angular frequency omega of the pulse (a.u.).
  real(R64) :: SysPotential_Ylm_CoulombLaser_omega
  !> Number of optical cycles under the sin^2 envelope (T = 2 pi nCycles / omega).
  real(R64) :: SysPotential_Ylm_CoulombLaser_nCycles
  !> Carrier-envelope phase (rad).
  real(R64) :: SysPotential_Ylm_CoulombLaser_cep
  !> Start time t0 of the pulse (a.u.).
  real(R64) :: SysPotential_Ylm_CoulombLaser_tStart

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
