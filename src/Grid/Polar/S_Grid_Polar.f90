! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Grid_Polar) S_Grid_Polar

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Grid_Polar_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Grid_Polar_Const

    call Say_Fabricate("grid.polar")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    ! Read radial grid parameters from JSON
    Grid_Polar_rmin = Json_Get("grid.polar.rmin", 0.0_R64)
    Grid_Polar_rmax = Json_Get("grid.polar.rmax", 20.0_R64)

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("grid.polar.const")) then
      call Grid_Polar_Const_Fabricate
    else
      error stop "grid.polar is missing one of: const"
    end if

  end subroutine

end submodule
