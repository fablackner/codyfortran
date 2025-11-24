! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> 2D polar grid (r, phi) parameters and runtime wiring.
!>
!> Provides arrays for coordinates and quadrature weights for a two-dimensional
!> polar grid where the radial coordinate is typically uniform or piecewise
!> uniform. The Fabricate routine initializes the discretization and wires
!> polar-specific operators as needed.
module M_Grid_Polar
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Instantiate and configure a 2D polar grid.
    !>
    !> Initializes radial and angular points and their integration weights,
    !> constructs Cartesian coordinate views, and associates polar operators.
    module subroutine Grid_Polar_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Minimum radial coordinate value rmin (> 0 recommended for kinetic operator).
  real(R64) :: Grid_Polar_rmin

  !> Maximum radial coordinate value rmax.
  real(R64) :: Grid_Polar_rmax

  !> Total integration weights for each (r,phi) grid point in flattened order.
  !> Combines radial and angular quadrature as used by the inner product.
  real(R64), allocatable :: Grid_Polar_weights(:)

  !> Polar coordinates (component-wise views) for each (r,phi) grid point.
  real(R64), allocatable :: Grid_Polar_rCoord(:)
  real(R64), allocatable :: Grid_Polar_phiCoord(:)

  !> Cartesian coordinates (component-wise views) for each (r,phi) grid point.
  real(R64), allocatable :: Grid_Polar_xCoord(:)
  real(R64), allocatable :: Grid_Polar_yCoord(:)

  !> Number of grid points along the radial direction (size of radialPoints).
  integer(I32) :: Grid_Polar_nRadial

  !> Number of grid points along the azimuthal angle phi.
  integer(I32) :: Grid_Polar_nPhi

  !> Radial quadrature nodes r_i.
  real(R64), allocatable :: Grid_Polar_radialPoints(:)

  !> Radial quadrature weights w_i^r.
  real(R64), allocatable :: Grid_Polar_radialWeights(:)

  !> Angular quadrature nodes phi_j.
  real(R64), allocatable :: Grid_Polar_phiPoints(:)

  !> Angular quadrature weights w_j^phi.
  real(R64), allocatable :: Grid_Polar_phiWeights(:)

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
