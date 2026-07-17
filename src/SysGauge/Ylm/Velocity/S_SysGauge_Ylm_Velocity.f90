! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysGauge_Ylm_Velocity) S_SysGauge_Ylm_Velocity

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysGauge_Ylm_Velocity_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysGauge
    use M_SysGauge_Ylm_Velocity_Fedvr
    use M_SysGauge_Ylm_Velocity_FedvrEcs

    implicit none

    call Say_Fabricate("sysGauge.ylm.velocity")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysGauge_Ylm_Velocity_fieldStrength = &
      Json_Get("sysGauge.ylm.velocity.fieldStrength", 0.01_R64)
    SysGauge_Ylm_Velocity_omega = &
      Json_Get("sysGauge.ylm.velocity.omega", 0.057_R64)
    SysGauge_Ylm_Velocity_nCycles = &
      Json_Get("sysGauge.ylm.velocity.nCycles", 3.0_R64)
    SysGauge_Ylm_Velocity_cep = &
      Json_Get("sysGauge.ylm.velocity.cep", 0.0_R64)
    SysGauge_Ylm_Velocity_tStart = &
      Json_Get("sysGauge.ylm.velocity.tStart", 0.0_R64)
    SysGauge_Ylm_Velocity_bodyMass = &
      Json_Get("sysGauge.ylm.velocity.bodyMass", [1.0_R64])

    SysGauge_timeIndependentQ = .false.
    SysGauge_bodyTypeIndependentQ = .true.

    SysGauge_MultiplyWithGaugeOp => MultiplyWithGaugeOp

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("sysGauge.ylm.velocity.fedvr")) then
      call SysGauge_Ylm_Velocity_Fedvr_Fabricate

    else if (Json_GetExistence("sysGauge.ylm.velocity.fedvrEcs")) then
      call SysGauge_Ylm_Velocity_FedvrEcs_Fabricate

    else
      error stop "sysGauge.ylm.velocity is missing one of: fedvr, fedvrEcs"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Vector potential A(t) of the sin²-envelope pulse (zero outside the pulse).
  module function SysGauge_Ylm_Velocity_VectorPotentialAmplitude(time) result(res)
    use M_Utils_Constants

    real(R64), intent(in) :: time
    real(R64) :: res

    real(R64) :: tRel, duration, envelope, carrier

    res = 0.0_R64

    if (SysGauge_Ylm_Velocity_omega <= 0.0_R64) return

    duration = 2.0_R64 * PI * SysGauge_Ylm_Velocity_nCycles / SysGauge_Ylm_Velocity_omega
    tRel = time - SysGauge_Ylm_Velocity_tStart

    if (tRel < 0.0_R64 .or. tRel > duration) return

    envelope = sin(PI * tRel / duration)**2
    carrier = sin(SysGauge_Ylm_Velocity_omega * (tRel - 0.5_R64 * duration) &
                  + SysGauge_Ylm_Velocity_cep)

    res = -SysGauge_Ylm_Velocity_fieldStrength / SysGauge_Ylm_Velocity_omega * envelope * carrier

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Apply H_A = A(t) p_z / m + A(t)²/(2m) by looping over all (l,m) channels.
  !>
  !> Each source channel f_{lm}(r) Y_{lm} contributes to the two neighboring
  !> angular momenta of the same m (see M_SysGauge_Ylm):
  !>   dOrb_{l+1,m} += (A/m) (-i) c_{l,m}   (f' − l f/r)
  !>   dOrb_{l−1,m} += (A/m) (-i) c_{l−1,m} (f' + (l+1) f/r)
  !> plus the diagonal A²/(2m) f on the source channel itself.
  !>
  !> Implementation trick (same as the kinetic operator): with g(r) = r·f(r),
  !>   f' − l f/r      = (1/r) (g' − (l+1) g/r)
  !>   f' + (l+1) f/r  = (1/r) (g' + l g/r)
  !> The r² metric is absorbed by the transform, so the weak-form derivative
  !> acting on g is exactly antisymmetric under the plain radial weights and
  !> the ±k/r pieces become exactly Hermitian diagonal multiplications: the
  !> discrete operator is Hermitian in the weighted metric to machine
  !> precision. One radial derivative per channel serves both terms.
  subroutine MultiplyWithGaugeOp(dOrb, orb, time, bt_)
    use M_Utils_Constants
    use M_Grid_Ylm
    use M_SysGauge_Ylm

    complex(R64), intent(out), contiguous :: dOrb(:)
    complex(R64), intent(in), contiguous  :: orb(:)
    real(R64), intent(in)                 :: time
    integer(I32), intent(in), optional    :: bt_

    integer(I32) :: l, m, lmax, nRad, bt
    real(R64) :: aPot, mass
    complex(R64), allocatable :: orbLm(:), dOrbLm(:), gRadial(:), gRadialDeriv(:)

    dOrb = 0.0_R64

    aPot = SysGauge_Ylm_Velocity_VectorPotentialAmplitude(time)
    if (aPot .eq. 0.0_R64) return

    bt = 1
    if (present(bt_)) bt = bt_
    mass = SysGauge_Ylm_Velocity_bodyMass(bt)

    nRad = Grid_Ylm_nRadial
    lmax = Grid_Ylm_lmax

    allocate (orbLm(nRad), dOrbLm(nRad), gRadial(nRad), gRadialDeriv(nRad))

    ! Process each (l,m) source channel separately
    do l = 0, lmax
      do m = -l, l
        ! Extract radial points for this (l,m) channel
        call Grid_Ylm_GetLmComponent(orbLm, l, m, orb)

        ! Skip negligible channels for efficiency
        if (all(abs(orbLm) < 1.0e-14_R64)) cycle

        ! Transform: g(r) = r · f(r), then dg/dr
        associate (radial => SysGauge_Ylm_radialCoordinates)

          gRadial(:) = radial(:) * orbLm(:)
          call SysGauge_Ylm_ApplyRadialFirstDerivative(gRadialDeriv, gRadial)

          ! Diagonal A²/(2m) term on the source channel
          dOrbLm(:) = 0.5_R64 * aPot**2 / mass * orbLm(:)
          call Grid_Ylm_AddLmComponent(dOrb, l, m, dOrbLm)

          ! Raising term: couples into (l+1, m)
          if (l + 1 <= lmax) then
            dOrbLm(:) = -IU * aPot / mass * CosThetaCoupling(l, m) * &
                        (gRadialDeriv(:) - (l + 1) * gRadial(:) / radial(:)) / radial(:)
            call Grid_Ylm_AddLmComponent(dOrb, l + 1, m, dOrbLm)
          end if

          ! Lowering term: couples into (l−1, m)
          if (l - 1 >= abs(m)) then
            dOrbLm(:) = -IU * aPot / mass * CosThetaCoupling(l - 1, m) * &
                        (gRadialDeriv(:) + l * gRadial(:) / radial(:)) / radial(:)
            call Grid_Ylm_AddLmComponent(dOrb, l - 1, m, dOrbLm)
          end if

        end associate

      end do
    end do

    deallocate (orbLm, dOrbLm, gRadial, gRadialDeriv)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Coupling coefficient c_{l,m} of cosθ Y_{lm} = c_{l,m} Y_{l+1,m} + c_{l−1,m} Y_{l−1,m}.
  pure function CosThetaCoupling(l, m) result(res)
    integer(I32), intent(in) :: l
    integer(I32), intent(in) :: m
    real(R64) :: res

    res = sqrt(real((l + 1)**2 - m**2, R64) / real((2 * l + 1) * (2 * l + 3), R64))

  end function

end submodule
