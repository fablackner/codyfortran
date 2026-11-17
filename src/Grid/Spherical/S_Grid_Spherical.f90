! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Grid_Spherical) S_Grid_Spherical

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Grid_Spherical_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Grid_Spherical_Const

    call Say_Fabricate("grid.spherical")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    ! Read radial grid parameters from JSON
    Grid_Spherical_rmin = Json_Get("grid.spherical.rmin", 0.0_R64)
    Grid_Spherical_rmax = Json_Get("grid.spherical.rmax", 20.0_R64)

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("grid.spherical.const")) then
      call Grid_Spherical_Const_Fabricate
    else
      error stop "grid.spherical is missing one of: const"
    end if

  end subroutine

end submodule
