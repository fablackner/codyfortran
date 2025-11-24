! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysPotential_Lattice_RandomUniform) S_SysPotential_Lattice_RandomUniform

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysPotential_Lattice_RandomUniform_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysPotential
    use M_SysPotential_Lattice_RandomUniform_StdImpl

    implicit none

    real(R64) :: minValue, maxValue
    integer(I32) :: seed

    call Say_Fabricate("sysPotential.lattice.randomUniform")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    minValue = Json_Get("sysPotential.lattice.randomUniform.minValue", -1.0_R64)
    maxValue = Json_Get("sysPotential.lattice.randomUniform.maxValue", 1.0_R64)
    seed = Json_Get("sysPotential.lattice.randomUniform.seed", -1_I32)

    SysPotential_Lattice_RandomUniform_minValue = minValue
    SysPotential_Lattice_RandomUniform_maxValue = maxValue
    SysPotential_Lattice_RandomUniform_seed = seed

    SysPotential_timeIndependentQ = .true.
    SysPotential_bodyTypeIndependentQ = .true.

    ! RNG seeding is handled in the implementation's Setup using the system clock

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("sysPotential.lattice.randomUniform.stdImpl")) then
      call SysPotential_Lattice_RandomUniform_StdImpl_Fabricate

    else
      error stop "sysPotential.lattice.randomUniform is missing one of: stdImpl"
    end if

  end subroutine

end submodule
