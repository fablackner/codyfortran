! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Linearly polarized laser pulse in length gauge (dipole approximation), the
!> model of sysGauge.fem.length and sysGauge.prolate.length.
!>
!> For an electron (charge -1, a.u.) in the electric field E(t) e_pol the
!> coupling is the diagonal operator
!>   H_E(t) = E(t) (e_pol . r),
!> applied point-wise at the grid nodes; e_pol . r is filled in Setup by the
!> coordinate filler of the grid family (Grid_Fem_FillNodeCoordinate,
!> Grid_Prolate_FillNodeCoordinate; r from the coordinate origin). The pulse
!> is that of sysPotential.ylm.coulombLaser (Pulse_Sin2Field of M_Utils_Pulse):
!>   E(t) = E0 sin^2(pi tRel/T) cos(omega (tRel - T/2) + cep),
!>   tRel = t - tStart,  T = 2 pi nCycles/omega,
!> and zero outside [tStart, tStart + T]; the carrier is centered on the
!> envelope maximum (cep = 0: cos-like pulse).
!>
!> Why a SysGauge backend although the operator is diagonal: SysPotential
!> stays static, so the exact field-free part (a fieldFree block) can
!> propagate h0 = T + V_static exactly while the integrator only sees this
!> non-stiff time-dependent term.
!>
!> JSON Configuration (<family>: fem or prolate)
!> ------------------
!>   "sysGauge": { "<family>": { "length": {
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
!> | `polarization`  | real(3) | [0, 0, 1]   | Direction (normalized, ≠ 0; prolate: z)  |
module M_SysGauge_Length
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  abstract interface
    !> coord = direction . r at the grid nodes (Grid_nPoints entries).
    subroutine I_SysGauge_Length_FillCoordinate(coord, direction)
      import :: R64
      real(R64), allocatable, intent(out) :: coord(:)
      real(R64), intent(in)               :: direction(3)
    end subroutine
  end interface

  interface
    !> Reads the pulse under path (e.g. "sysGauge.fem.length"), validates it,
    !> sets the flags and binds Setup and the gauge operator; FillCoordinate
    !> gives e_pol . r in Setup.
    module subroutine SysGauge_Length_Fabricate(path, FillCoordinate)
      character(len=*), intent(in) :: path
      procedure(I_SysGauge_Length_FillCoordinate) :: FillCoordinate
    end subroutine

    !> Electric field amplitude E(t) along the polarization (a.u.), zero
    !> outside the pulse (for printers and observables).
    module function SysGauge_Length_FieldAmplitude(time) result(res)
      real(R64), intent(in) :: time
      real(R64) :: res
    end function
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Peak electric field strength E0 of the pulse (a.u.).
  real(R64) :: SysGauge_Length_fieldStrength = 0.0_R64
  !> Carrier angular frequency omega of the pulse (a.u.).
  real(R64) :: SysGauge_Length_omega = 0.0_R64
  !> Number of optical cycles under the sin^2 envelope (T = 2 pi nCycles / omega).
  real(R64) :: SysGauge_Length_nCycles = 0.0_R64
  !> Carrier-envelope phase (rad).
  real(R64) :: SysGauge_Length_cep = 0.0_R64
  !> Start time of the pulse (a.u.).
  real(R64) :: SysGauge_Length_tStart = 0.0_R64
  !> Unit polarization vector e_pol of the field.
  real(R64) :: SysGauge_Length_polarization(3) = [0.0_R64, 0.0_R64, 1.0_R64]

end module
