! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> 2D Cartesian grid (x,y) parameters and runtime wiring.
!>
!> Stores coordinates and weights for a two-dimensional Cartesian grid with
!> uniform or piecewise-uniform spacing. A Fabricate routine initializes the
!> discretization and wires 2D operators (e.g., finite differences, FFT).
module M_Grid_Square
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Instantiate and configure a 2D Cartesian grid back-end.
    module subroutine Grid_Square_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Minimum x-coordinate (left boundary of the domain).
  real(R64)  ::   Grid_Square_xmin

  !> Maximum x-coordinate (right boundary of the domain).
  real(R64)  ::   Grid_Square_xmax

  !> Minimum y-coordinate (bottom boundary of the domain).
  real(R64)  ::   Grid_Square_ymin

  !> Maximum y-coordinate (top boundary of the domain).
  real(R64)  ::   Grid_Square_ymax

  !> Total integration weights for each (x,y) grid point in flattened order.
  !> Combines per-dimension quadrature as used by the inner product.
  real(R64), allocatable :: Grid_Square_weights(:)

  !> Component-wise Cartesian coordinates for each grid point (flattened order).
  real(R64), allocatable :: Grid_Square_xCoord(:)
  real(R64), allocatable :: Grid_Square_yCoord(:)

  !> Number of grid points in the x-dimension.
  integer(I32) :: Grid_Square_nPointsX

  !> Number of grid points in the y-dimension.
  integer(I32) :: Grid_Square_nPointsY

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
