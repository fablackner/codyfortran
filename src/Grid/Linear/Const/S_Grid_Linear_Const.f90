! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Grid_Linear_Const) S_Grid_Linear_Const

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Grid_Linear_Const_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Grid

    call Say_Fabricate("grid.linear.const")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    Grid_nPoints = Json_Get("grid.linear.const.nPoints", 100)

    Grid_Setup => Setup

    !------------------------------------
    ! branch
    !------------------------------------

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Setup
    use M_Utils_Json
    use M_Utils_Say
    use M_Grid
    use M_Grid_Linear

    integer(I32) :: iGrid

    call Say_Setup("grid.linear.const")

    allocate (Grid_Linear_weights(Grid_nPoints))
    Grid_Linear_weights(:) = (Grid_Linear_xmax - Grid_Linear_xmin) / (Grid_nPoints - 1)

    ! Allocate coordinates for 1D
    allocate (Grid_Linear_xCoord(Grid_nPoints))
    do iGrid = 1, Grid_nPoints
      Grid_Linear_xCoord(iGrid) = Grid_Linear_xmin + Grid_Linear_weights(iGrid) * (iGrid - 1)
    end do
  end subroutine

end submodule
