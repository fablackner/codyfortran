! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Grid_Spherical_Const) S_Grid_Spherical_Const

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Grid_Spherical_Const_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Grid
    use M_Grid_Spherical

    call Say_Fabricate("grid.spherical.constRadial3d")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    ! Read radial grid parameters from JSON
    Grid_Spherical_nRadial = Json_Get("grid.spherical.constRadial3d.nRadial", 100)
    Grid_Spherical_nTheta = Json_Get("grid.spherical.constRadial3d.nTheta", 21)
    Grid_Spherical_nPhi = Json_Get("grid.spherical.constRadial3d.nPhi", 42)

    ! Set total grid size - product of radial, theta and phi points
    Grid_nPoints = Grid_Spherical_nRadial * Grid_Spherical_nTheta * Grid_Spherical_nPhi

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
    use M_Utils_Constants
    use M_Grid
    use M_Grid_Spherical

    integer(I32) :: iRad, iTheta, iPhi, iGrid
    real(R64) :: radialPos, theta, phi, dr

    call Say_Setup("grid.spherical.const")

    ! Allocate arrays for the grid
    allocate (Grid_Spherical_xCoord(Grid_nPoints))  ! 3D coordinates (x,y,z)
    allocate (Grid_Spherical_yCoord(Grid_nPoints))  ! 3D coordinates (x,y,z)
    allocate (Grid_Spherical_zCoord(Grid_nPoints))  ! 3D coordinates (x,y,z)
    allocate (Grid_Spherical_rCoord(Grid_nPoints))     ! Spherical coordinates (r,θ,φ)
    allocate (Grid_Spherical_thetaCoord(Grid_nPoints))     ! Spherical coordinates (r,θ,φ)
    allocate (Grid_Spherical_phiCoord(Grid_nPoints))     ! Spherical coordinates (r,θ,φ)
    allocate (Grid_Spherical_weights(Grid_nPoints))      ! Integration weights

    ! Allocate module arrays for radial and phi grid
    allocate (Grid_Spherical_radialPoints(Grid_Spherical_nRadial))
    allocate (Grid_Spherical_radialWeights(Grid_Spherical_nRadial))
    allocate (Grid_Spherical_thetaPoints(Grid_Spherical_nTheta))
    allocate (Grid_Spherical_thetaWeights(Grid_Spherical_nTheta))
    allocate (Grid_Spherical_phiPoints(Grid_Spherical_nPhi))
    allocate (Grid_Spherical_phiWeights(Grid_Spherical_nPhi))

    ! Set up the radial grid spacing
    dr = (Grid_Spherical_rmax - Grid_Spherical_rmin) / real(Grid_Spherical_nRadial, R64)

    ! Set up the radial grid
    do iRad = 1, Grid_Spherical_nRadial
      ! Calculate radial points starting from rmin + dr
      Grid_Spherical_radialPoints(iRad) = Grid_Spherical_rmin + iRad * dr

      ! For 3D spherical coordinates, radial weight includes r²
      Grid_Spherical_radialWeights(iRad) = dr * Grid_Spherical_radialPoints(iRad)**2
    end do

    ! Set up theta grid (polar angle, 0 to π)
    do iTheta = 1, Grid_Spherical_nTheta
      Grid_Spherical_thetaPoints(iTheta) = PI * (iTheta) / real(Grid_Spherical_nTheta + 1, R64)
      ! Weight includes sin(theta) for spherical coordinates
      Grid_Spherical_thetaWeights(iTheta) = PI / Grid_Spherical_nTheta * &
                                            sin(Grid_Spherical_thetaPoints(iTheta))
    end do

    ! Set up phi grid (azimuthal angle, 0 to 2π)
    do iPhi = 1, Grid_Spherical_nPhi
      Grid_Spherical_phiPoints(iPhi) = 2.0_R64 * PI * (iPhi - 1) / real(Grid_Spherical_nPhi, R64)
      Grid_Spherical_phiWeights(iPhi) = 2.0_R64 * PI / Grid_Spherical_nPhi
    end do

    ! Now combine radial, theta and phi to form the complete 3D grid
    ! Loop order: phi (outer), theta (middle), radial (inner)
    ! This makes radial index vary fastest in memory
    iGrid = 0
    do iPhi = 1, Grid_Spherical_nPhi
      phi = Grid_Spherical_phiPoints(iPhi)

      do iTheta = 1, Grid_Spherical_nTheta
        theta = Grid_Spherical_thetaPoints(iTheta)

        do iRad = 1, Grid_Spherical_nRadial
          radialPos = Grid_Spherical_radialPoints(iRad)

          iGrid = iGrid + 1

          ! Store spherical coordinates
          Grid_Spherical_rCoord(iGrid) = radialPos   ! r
          Grid_Spherical_thetaCoord(iGrid) = theta       ! θ
          Grid_Spherical_phiCoord(iGrid) = phi         ! φ

          ! Convert from spherical (r,θ,φ) to Cartesian (x,y,z) coordinates
          Grid_Spherical_xCoord(iGrid) = radialPos * sin(theta) * cos(phi)  ! x
          Grid_Spherical_yCoord(iGrid) = radialPos * sin(theta) * sin(phi)  ! y
          Grid_Spherical_zCoord(iGrid) = radialPos * cos(theta)             ! z

          ! Set integration weight as product of radial, theta and phi weights
          Grid_Spherical_weights(iGrid) = Grid_Spherical_radialWeights(iRad) * &
                                          Grid_Spherical_thetaWeights(iTheta) * &
                                          Grid_Spherical_phiWeights(iPhi)
        end do
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end submodule
