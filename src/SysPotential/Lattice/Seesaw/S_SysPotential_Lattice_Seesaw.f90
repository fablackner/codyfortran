! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysPotential_Lattice_Seesaw) S_SysPotential_Lattice_Seesaw

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysPotential_Lattice_Seesaw_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysPotential
    use M_SysPotential_Lattice_Seesaw_StdImpl

    implicit none

    real(R64) :: slopeMinX, slopeMinY, slopeMinZ
    real(R64) :: slopeMaxX, slopeMaxY, slopeMaxZ
    real(R64) :: frequencyX, frequencyY, frequencyZ

    call Say_Fabricate("sysPotential.lattice.seesaw")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    slopeMinX = Json_Get("sysPotential.lattice.seesaw.slopeMinX", 0.0_R64)
    slopeMinY = Json_Get("sysPotential.lattice.seesaw.slopeMinY", 0.0_R64)
    slopeMinZ = Json_Get("sysPotential.lattice.seesaw.slopeMinZ", 0.0_R64)
    slopeMaxX = Json_Get("sysPotential.lattice.seesaw.slopeMaxX", 0.0_R64)
    slopeMaxY = Json_Get("sysPotential.lattice.seesaw.slopeMaxY", 0.0_R64)
    slopeMaxZ = Json_Get("sysPotential.lattice.seesaw.slopeMaxZ", 0.0_R64)
    frequencyX = Json_Get("sysPotential.lattice.seesaw.frequencyX", 0.0_R64)
    frequencyY = Json_Get("sysPotential.lattice.seesaw.frequencyY", 0.0_R64)
    frequencyZ = Json_Get("sysPotential.lattice.seesaw.frequencyZ", 0.0_R64)

    SysPotential_Lattice_Seesaw_slopeMinX = slopeMinX
    SysPotential_Lattice_Seesaw_slopeMinY = slopeMinY
    SysPotential_Lattice_Seesaw_slopeMinZ = slopeMinZ
    SysPotential_Lattice_Seesaw_slopeMaxX = slopeMaxX
    SysPotential_Lattice_Seesaw_slopeMaxY = slopeMaxY
    SysPotential_Lattice_Seesaw_slopeMaxZ = slopeMaxZ
    SysPotential_Lattice_Seesaw_frequencyX = frequencyX
    SysPotential_Lattice_Seesaw_frequencyY = frequencyY
    SysPotential_Lattice_Seesaw_frequencyZ = frequencyZ

    SysPotential_timeIndependentQ = .false.
    SysPotential_bodyTypeIndependentQ = .true.

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("sysPotential.lattice.seesaw.stdImpl")) then
      call SysPotential_Lattice_Seesaw_StdImpl_Fabricate

    else
      error stop "sysPotential.lattice.seesaw is missing one of: stdImpl"
    end if

  end subroutine

end submodule
