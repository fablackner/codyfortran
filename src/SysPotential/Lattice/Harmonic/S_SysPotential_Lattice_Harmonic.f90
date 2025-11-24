! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysPotential_Lattice_Harmonic) S_SysPotential_Lattice_Harmonic

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysPotential_Lattice_Harmonic_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysPotential
    use M_SysPotential_Lattice_Harmonic_StdImpl

    implicit none

    real(R64) :: omegaX
    real(R64) :: omegaY
    real(R64) :: omegaZ

    real(R64) :: positionX
    real(R64) :: positionY
    real(R64) :: positionZ

    call Say_Fabricate("sysPotential.lattice.harmonic")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    omegaX = Json_Get("sysPotential.lattice.harmonic.omegaX", 0.0_R64)
    omegaY = Json_Get("sysPotential.lattice.harmonic.omegaY", 0.0_R64)
    omegaZ = Json_Get("sysPotential.lattice.harmonic.omegaZ", 0.0_R64)

    SysPotential_Lattice_Harmonic_omegaX = omegaX
    SysPotential_Lattice_Harmonic_omegaY = omegaY
    SysPotential_Lattice_Harmonic_omegaZ = omegaZ

    positionX = Json_Get("sysPotential.lattice.harmonic.positionX", 0.0_R64)
    positionY = Json_Get("sysPotential.lattice.harmonic.positionY", 0.0_R64)
    positionZ = Json_Get("sysPotential.lattice.harmonic.positionZ", 0.0_R64)

    SysPotential_Lattice_Harmonic_positionX = positionX
    SysPotential_Lattice_Harmonic_positionY = positionY
    SysPotential_Lattice_Harmonic_positionZ = positionZ

    SysPotential_timeIndependentQ = .true.
    SysPotential_bodyTypeIndependentQ = .true.

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("sysPotential.lattice.harmonic.stdImpl")) then
      call SysPotential_Lattice_Harmonic_StdImpl_Fabricate

    else
      error stop "sysPotential.lattice.harmonic is missing one of: stdImpl"
    end if

  end subroutine

end submodule
