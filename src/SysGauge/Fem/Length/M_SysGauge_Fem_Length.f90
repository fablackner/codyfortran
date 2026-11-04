! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Linearly polarized laser pulse on grid.fem in length gauge (dipole
!> approximation).
!>
!> For an electron (charge -1, a.u.) in the electric field E(t) e_pol the
!> coupling is the diagonal operator
!>   H_E(t) = E(t) (e_pol . r),
!> applied point-wise at the nodes (e_pol . r is cached in Setup by
!> Grid_Fem_FillNodeCoordinate; r is measured from the coordinate origin).
!> The pulse is that of sysPotential.ylm.coulombLaser (Pulse_Sin2Field of
!> M_Utils_Pulse):
!>   E(t) = E0 sin^2(pi tRel/T) cos(omega (tRel - T/2) + cep),
!>   tRel = t - tStart,  T = 2 pi nCycles/omega,
!> and zero outside [tStart, tStart + T]; the carrier is centered on the
!> envelope maximum (cep = 0: cos-like pulse).
!>
!> Why a SysGauge backend although the operator is diagonal: SysPotential on
!> grid.fem stays static, so the exact field-free part
!> (method.mb.orbBased.exactFieldFreeQ) can propagate h0 = T + V_static exactly
!> while the integrator only sees this non-stiff time-dependent term.
!>
!> JSON Configuration
!> ------------------
!>   "sysGauge": { "fem": { "length": {
!>       "fieldStrength": 0.01, "omega": 0.057, "nCycles": 3.0, "cep": 0.0,
!>       "tStart": 0.0, "polarization": [0.0, 0.0, 1.0]
!>   } } }
!>
!> | Key             | Type    | Default     | Description                              |
!> |-----------------|---------|-------------|------------------------------------------|
!> | `fieldStrength` | real    | 0.01        | Peak field E0 (a.u.)                     |
!> | `omega`         | real    | 0.057       | Carrier frequency (a.u.), > 0            |
!> | `nCycles`       | real    | 3.0         | Cycles under the sin^2 envelope, > 0     |
!> | `cep`           | real    | 0.0         | Carrier-envelope phase (rad)             |
!> | `tStart`        | real    | 0.0         | Pulse start (a.u.)                       |
!> | `polarization`  | real(3) | [0, 0, 1]   | Polarization direction (normalized, ≠ 0) |
module M_SysGauge_Fem_Length
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Read the pulse, validate it, set the flags and bind Setup and the
    !> gauge operator.
    module subroutine SysGauge_Fem_Length_Fabricate
    end subroutine

    !> Electric field amplitude E(t) along the polarization (a.u.), zero
    !> outside the pulse (for printers and observables).
    module function SysGauge_Fem_Length_FieldAmplitude(time) result(res)
      real(R64), intent(in) :: time
      real(R64) :: res
    end function
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Peak electric field strength E0 of the pulse (a.u.).
  real(R64) :: SysGauge_Fem_Length_fieldStrength = 0.0_R64
  !> Carrier angular frequency omega of the pulse (a.u.).
  real(R64) :: SysGauge_Fem_Length_omega = 0.0_R64
  !> Number of optical cycles under the sin^2 envelope (T = 2 pi nCycles / omega).
  real(R64) :: SysGauge_Fem_Length_nCycles = 0.0_R64
  !> Carrier-envelope phase (rad).
  real(R64) :: SysGauge_Fem_Length_cep = 0.0_R64
  !> Start time of the pulse (a.u.).
  real(R64) :: SysGauge_Fem_Length_tStart = 0.0_R64
  !> Unit polarization vector e_pol of the field.
  real(R64) :: SysGauge_Fem_Length_polarization(3) = [0.0_R64, 0.0_R64, 1.0_R64]

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
