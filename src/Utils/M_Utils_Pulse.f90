! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Laser-pulse shapes (stateless, pure).
!>
!> Pulse_Sin2Field is the electric field of the length-gauge laser couplings
!> sysPotential.ylm.coulombLaser and sysGauge.fem.length:
!>   E(t) = E0 sin^2(pi tRel/T) cos(omega (tRel - T/2) + cep),
!>   tRel = t - tStart,  T = 2 pi nCycles/omega,
!> zero outside [tStart, tStart + T] and for omega <= 0. The carrier is
!> centered on the envelope maximum (cep = 0: cos-like pulse).
!>
!> The velocity gauge sysGauge.ylm.velocity keeps its own shapes: there the
!> envelope sits on the vector potential,
!>   A(t) = -(E0/omega) sin^2(pi tRel/T) sin(omega (tRel - T/2) + cep),
!> or, with envelope "field", on E(t) with the uncentered carrier
!> cos(omega tRel + cep), so neither is this formula.
module M_Utils_Pulse
  use M_Utils_Types

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Electric field E(t) of the sin^2-envelope pulse (a.u.):
  !> amplitude sin^2(pi tRel/T) cos(omega (tRel - T/2) + cep), tRel = time -
  !> tStart, T = 2 pi nCycles/omega; zero outside [tStart, tStart + T] and for
  !> omega <= 0.
  pure function Pulse_Sin2Field(time, amplitude, omega, nCycles, cep, tStart) result(res)
    use M_Utils_Constants, only: PI

    !> Time t (a.u.).
    real(R64), intent(in) :: time
    !> Peak field E0 (a.u.).
    real(R64), intent(in) :: amplitude
    !> Carrier angular frequency omega (a.u.); no pulse for omega <= 0.
    real(R64), intent(in) :: omega
    !> Optical cycles under the envelope (T = 2 pi nCycles/omega).
    real(R64), intent(in) :: nCycles
    !> Carrier-envelope phase (rad).
    real(R64), intent(in) :: cep
    !> Start time of the pulse (a.u.).
    real(R64), intent(in) :: tStart
    !> E(t) (a.u.).
    real(R64) :: res

    real(R64) :: tRel, duration, envelope, carrier

    res = 0.0_R64

    if (omega <= 0.0_R64) return

    duration = 2.0_R64 * PI * nCycles / omega
    tRel = time - tStart

    if (tRel < 0.0_R64 .or. tRel > duration) return

    envelope = sin(PI * tRel / duration)**2
    carrier = cos(omega * (tRel - 0.5_R64 * duration) + cep)

    res = amplitude * envelope * carrier

  end function

end module
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
