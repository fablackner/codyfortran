! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Grid_Square) S_Grid_Square

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Grid_Square_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Grid_Square_Const

    call Say_Fabricate("grid.square")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    ! Set domain boundaries for both dimensions
    Grid_Square_xmin = Json_Get("grid.square.xmin", -20.0_R64)
    Grid_Square_xmax = Json_Get("grid.square.xmax", 20.0_R64)
    Grid_Square_ymin = Json_Get("grid.square.ymin", -20.0_R64)
    Grid_Square_ymax = Json_Get("grid.square.ymax", 20.0_R64)

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("grid.square.const")) then
      call Grid_Square_Const_Fabricate

    else
      error stop "grid.square is missing one of: const"
    end if

  end subroutine

end submodule
