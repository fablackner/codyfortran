! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysPotential_Lattice_RandomGauss) S_SysPotential_Lattice_RandomGauss

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysPotential_Lattice_RandomGauss_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysPotential
    use M_SysPotential_Lattice_RandomGauss_StdImpl

    implicit none

    real(R64) :: meanValue, stdValue
    integer(I32) :: seed

    call Say_Fabricate("sysPotential.lattice.randomGauss")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    meanValue = Json_Get("sysPotential.lattice.randomGauss.meanValue", 0.0_R64)
    stdValue = Json_Get("sysPotential.lattice.randomGauss.stdImplValue", 1.0_R64)
    seed = Json_Get("sysPotential.lattice.randomGauss.seed", -1_I32)

    SysPotential_Lattice_RandomGauss_meanValue = meanValue
    SysPotential_Lattice_RandomGauss_stdValue = stdValue
    SysPotential_Lattice_RandomGauss_seed = seed

    SysPotential_timeIndependentQ = .true.
    SysPotential_bodyTypeIndependentQ = .true.

    ! RNG seeding is handled in the implementation's Setup using the system clock

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("sysPotential.lattice.randomGauss.stdImpl")) then
      call SysPotential_Lattice_RandomGauss_StdImpl_Fabricate

    else
      error stop "sysPotential.lattice.randomGauss is missing one of: stdImpl"
    end if

  end subroutine

end submodule
