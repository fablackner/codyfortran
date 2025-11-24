! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Grid_Polar_Const) S_Grid_Polar_Const

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Grid_Polar_Const_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Grid
    use M_Grid_Polar

    call Say_Fabricate("grid.polar.constRadial2d")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    ! Read radial grid parameters from JSON
    Grid_Polar_nRadial = Json_Get("grid.polar.constRadial2d.nRadial", 100)
    Grid_Polar_rmin = Json_Get("grid.polar.constRadial2d.rmin", 0.0_R64)
    Grid_Polar_rmax = Json_Get("grid.polar.constRadial2d.rmax", 20.0_R64)

    ! For a standard polar grid, we just need the number of phi points
    ! No need for mmax parameter
    Grid_Polar_nPhi = Json_Get("grid.polar.constRadial2d.nPhi", 21)

    ! Set total grid size - product of radial points and phi points
    Grid_nPoints = Grid_Polar_nRadial * Grid_Polar_nPhi

    ! Set the setup procedure
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
    use M_Grid_Polar
    use M_Utils_Constants

    integer(I32) :: iRad, iPhi, iGrid
    real(R64) :: radialPos, phi, dr

    call Say_Setup("grid.polar.const")

    ! Allocate arrays for the grid
    allocate (Grid_Polar_xCoord(Grid_nPoints))  ! 2D coordinates
    allocate (Grid_Polar_yCoord(Grid_nPoints))  ! 2D coordinates
    allocate (Grid_Polar_weights(Grid_nPoints))      ! Integration weights
    allocate (Grid_Polar_rCoord(Grid_nPoints))  ! Polar coordinates (r, phi)
    allocate (Grid_Polar_phiCoord(Grid_nPoints))  ! Polar coordinates (r, phi)

    ! Allocate module arrays for radial and phi grid
    allocate (Grid_Polar_radialPoints(Grid_Polar_nRadial))
    allocate (Grid_Polar_radialWeights(Grid_Polar_nRadial))
    allocate (Grid_Polar_phiPoints(Grid_Polar_nPhi))
    allocate (Grid_Polar_phiWeights(Grid_Polar_nPhi))

    ! Set up the radial grid spacing
    dr = (Grid_Polar_rmax - Grid_Polar_rmin) / Grid_Polar_nRadial

    ! Set up the radial grid (assuming origin is not included)
    do iRad = 1, Grid_Polar_nRadial
      ! Calculate radial points starting from rmin + dr
      ! Ensure points are strictly positive if rmin is 0
      Grid_Polar_radialPoints(iRad) = Grid_Polar_rmin + iRad * dr

      ! For 2D polar coordinates, the integration weight includes r (not r²)
      Grid_Polar_radialWeights(iRad) = dr * Grid_Polar_radialPoints(iRad)
    end do

    ! Set up phi grid with uniform spacing
    do iPhi = 1, Grid_Polar_nPhi
      Grid_Polar_phiPoints(iPhi) = 2.0_R64 * PI * (iPhi - 1) / Grid_Polar_nPhi
      Grid_Polar_phiWeights(iPhi) = 2.0_R64 * PI / Grid_Polar_nPhi
    end do

    ! Now combine radial and phi parts to form the complete 2D grid
    ! Changed loop order: outer loop is phi, inner loop is radial
    ! This makes radial index vary fastest in memory
    iGrid = 0
    do iPhi = 1, Grid_Polar_nPhi
      phi = Grid_Polar_phiPoints(iPhi)

      do iRad = 1, Grid_Polar_nRadial
        radialPos = Grid_Polar_radialPoints(iRad)

        iGrid = iGrid + 1

        ! Convert from polar (r,θ) to Cartesian (x,y,z) coordinates
        Grid_Polar_xCoord(iGrid) = radialPos * cos(phi)  ! x
        Grid_Polar_yCoord(iGrid) = radialPos * sin(phi)  ! y

        ! Store polar coordinates
        Grid_Polar_rCoord(iGrid) = radialPos  ! r
        Grid_Polar_phiCoord(iGrid) = phi        ! phi

        ! Set integration weight as product of radial and phi weights
        Grid_Polar_weights(iGrid) = Grid_Polar_radialWeights(iRad) * Grid_Polar_phiWeights(iPhi)
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end submodule
