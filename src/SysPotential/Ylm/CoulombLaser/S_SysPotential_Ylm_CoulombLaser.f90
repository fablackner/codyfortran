! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysPotential_Ylm_CoulombLaser) S_SysPotential_Ylm_CoulombLaser

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysPotential_Ylm_CoulombLaser_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysPotential
    use M_SysPotential_Ylm
    use M_SysPotential_Ylm_CoulombLaser_StdImpl

    implicit none

    call Say_Fabricate("sysPotential.ylm.coulombLaser")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysPotential_Ylm_CoulombLaser_charge = &
      Json_Get("sysPotential.ylm.coulombLaser.charge", 2.0_R64)
    SysPotential_Ylm_CoulombLaser_fieldStrength = &
      Json_Get("sysPotential.ylm.coulombLaser.fieldStrength", 0.01_R64)
    SysPotential_Ylm_CoulombLaser_omega = &
      Json_Get("sysPotential.ylm.coulombLaser.omega", 0.057_R64)
    SysPotential_Ylm_CoulombLaser_nCycles = &
      Json_Get("sysPotential.ylm.coulombLaser.nCycles", 3.0_R64)
    SysPotential_Ylm_CoulombLaser_cep = &
      Json_Get("sysPotential.ylm.coulombLaser.cep", 0.0_R64)
    SysPotential_Ylm_CoulombLaser_tStart = &
      Json_Get("sysPotential.ylm.coulombLaser.tStart", 0.0_R64)

    SysPotential_timeIndependentQ = .false.
    SysPotential_bodyTypeIndependentQ = .true.
    SysPotential_Ylm_lmax = 1
    SysPotential_Ylm_mIndependentQ = .true.

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("sysPotential.ylm.coulombLaser.stdImpl")) then
      call SysPotential_Ylm_CoulombLaser_StdImpl_Fabricate

    else
      error stop "sysPotential.ylm.coulombLaser is missing one of: stdImpl"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Electric field E(t) of the sin^2-envelope pulse (zero outside the pulse).
  module function SysPotential_Ylm_CoulombLaser_FieldAmplitude(time) result(res)
    use M_Utils_Constants

    real(R64), intent(in) :: time
    real(R64) :: res

    real(R64) :: tRel, duration, envelope, carrier

    res = 0.0_R64

    if (SysPotential_Ylm_CoulombLaser_omega <= 0.0_R64) return

    duration = 2.0_R64 * PI * SysPotential_Ylm_CoulombLaser_nCycles / SysPotential_Ylm_CoulombLaser_omega
    tRel = time - SysPotential_Ylm_CoulombLaser_tStart

    if (tRel < 0.0_R64 .or. tRel > duration) return

    envelope = sin(PI * tRel / duration)**2
    carrier = cos(SysPotential_Ylm_CoulombLaser_omega * (tRel - 0.5_R64 * duration) &
                  + SysPotential_Ylm_CoulombLaser_cep)

    res = SysPotential_Ylm_CoulombLaser_fieldStrength * envelope * carrier

  end function

end submodule
