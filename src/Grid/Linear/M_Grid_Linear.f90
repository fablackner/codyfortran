! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Common state for 1D linear (x) grids and runtime wiring.
!>
!> M_Grid_Linear collects generic 1D grid parameters used by multiple
!> implementations (constant spacing, FEDVR, …). A Fabricate routine selects
!> and initializes a concrete back-end and populates the arrays below.
module M_Grid_Linear
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Instantiate and configure a 1D linear grid back-end.
    !>
    !> Selects the concrete representation (e.g., constant spacing, FEDVR),
    !> initializes coordinates and weights, and wires any 1D-specific operators.
    module subroutine Grid_Linear_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Minimum x-coordinate of the 1D grid (left boundary of the domain).
  real(R64) :: Grid_Linear_xmin

  !> Maximum x-coordinate of the 1D grid (right boundary of the domain).
  real(R64) :: Grid_Linear_xmax

  !> Quadrature/metric weights for the 1D grid points used in inner products.
  !> Grid_Linear_weights(i) corresponds to Grid_Linear_xCoord(i).
  real(R64), allocatable :: Grid_Linear_weights(:)

  !> x-coordinates of the grid points ordered as used by the active back-end.
  !> Grid_Linear_xCoord(i) gives the position of point i.
  real(R64), allocatable :: Grid_Linear_xCoord(:)

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
