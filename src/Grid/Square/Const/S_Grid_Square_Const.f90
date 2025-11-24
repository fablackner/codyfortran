! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Grid_Square_Const) S_Grid_Square_Const

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Grid_Square_Const_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Grid
    use M_Grid_Square

    call Say_Fabricate("grid.square.const")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    ! Read grid dimensions for both x and y directions
    Grid_Square_nPointsX = Json_Get("grid.square.const.nPointsX", 100)
    Grid_Square_nPointsY = Json_Get("grid.square.const.nPointsY", 100)

    ! Total grid points is the product of points in each dimension
    Grid_nPoints = Grid_Square_nPointsX * Grid_Square_nPointsY

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
    use M_Grid_Square

    integer(I32) :: iX, iY, iGrid
    real(R64) :: dx, dy

    call Say_Setup("grid.square.const")

    ! Calculate step sizes
    dx = (Grid_Square_xmax - Grid_Square_xmin) / (Grid_Square_nPointsX - 1)
    dy = (Grid_Square_ymax - Grid_Square_ymin) / (Grid_Square_nPointsY - 1)

    ! Allocate and calculate weights (for 2D integration)
    allocate (Grid_Square_weights(Grid_nPoints))
    Grid_Square_weights(:) = dx * dy

    ! Allocate arrays for coordinate storage and indexing
    allocate (Grid_Square_xCoord(Grid_nPoints)) ! Store x, y coordinates
    allocate (Grid_Square_yCoord(Grid_nPoints)) ! Store x, y coordinates

    ! Set up the grid coordinates and indexing
    iGrid = 0
    do iY = 1, Grid_Square_nPointsY
      do iX = 1, Grid_Square_nPointsX
        iGrid = iGrid + 1

        ! Store coordinates
        Grid_Square_xCoord(iGrid) = Grid_Square_xmin + (iX - 1) * dx
        Grid_Square_yCoord(iGrid) = Grid_Square_ymin + (iY - 1) * dy

      end do
    end do

  end subroutine

end submodule
