! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysKinetic_Lattice) S_SysKinetic_Lattice

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysKinetic_Lattice_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysKinetic_Lattice_NearestNeighbor

    call Say_Fabricate("sysKinetic.lattice")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("sysKinetic.lattice.nearestNeighbor")) then
      call SysKinetic_Lattice_NearestNeighbor_Fabricate

    else
      error stop "sysKinetic.lattice is missing one of: nearestNeighbor"
    end if

  end subroutine

end submodule
