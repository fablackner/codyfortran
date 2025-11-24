! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysPotential_Ylm_Coulomb) S_SysPotential_Ylm_Coulomb

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysPotential_Ylm_Coulomb_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysPotential
    use M_SysPotential_Ylm
    use M_SysPotential_Ylm_Coulomb_StdImpl

    implicit none

    real(R64) :: charge

    call Say_Fabricate("sysPotential.ylm.coulomb")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    charge = Json_Get("sysPotential.ylm.coulomb.charge", 2.0_R64)

    SysPotential_Ylm_Coulomb_charge = charge
    SysPotential_timeIndependentQ = .true.
    SysPotential_bodyTypeIndependentQ = .true.
    SysPotential_Ylm_lmax = 0
    SysPotential_Ylm_mIndependentQ = .true.

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("sysPotential.ylm.coulomb.stdImpl")) then
      call SysPotential_Ylm_Coulomb_StdImpl_Fabricate

    else
      error stop "sysPotential.ylm.coulomb is missing one of: stdImpl"
    end if

  end subroutine

end submodule
