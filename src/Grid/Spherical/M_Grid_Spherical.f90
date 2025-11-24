! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> 3D spherical grid (r, theta, phi) parameters and runtime wiring.
!>
!> Holds coordinates and quadrature weights for a spherical grid and exposes
!> sizes for each angular and radial dimension. A Fabricate routine initializes
!> the discretization and wires spherical-specific operators.
module M_Grid_Spherical
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Instantiate and configure a 3D spherical grid back-end.
    module subroutine Grid_Spherical_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Minimum radius rmin (inner boundary of the domain).
  real(R64) :: Grid_Spherical_rmin

  !> Maximum radius rmax (outer boundary of the domain).
  real(R64) :: Grid_Spherical_rmax

  !> Total integration weights for each (r,theta,phi) grid point in flattened order.
  !> Combines radial and angular quadrature per the active spherical back-end.
  real(R64), allocatable :: Grid_Spherical_weights(:)

  !> Spherical coordinate components (component-wise views) for each grid point.
  real(R64), allocatable :: Grid_Spherical_rCoord(:)
  real(R64), allocatable :: Grid_Spherical_thetaCoord(:)
  real(R64), allocatable :: Grid_Spherical_phiCoord(:)

  !> Cartesian coordinates (component-wise views) for each spherical grid point.
  real(R64), allocatable :: Grid_Spherical_xCoord(:)
  real(R64), allocatable :: Grid_Spherical_yCoord(:)
  real(R64), allocatable :: Grid_Spherical_zCoord(:)

  !> Number of radial grid points.
  integer(I32) :: Grid_Spherical_nRadial

  !> Number of theta (polar angle) grid points.
  integer(I32) :: Grid_Spherical_nTheta

  !> Number of phi (azimuthal angle) grid points.
  integer(I32) :: Grid_Spherical_nPhi

  !> Radial quadrature nodes r_i.
  real(R64), allocatable :: Grid_Spherical_radialPoints(:)

  !> Radial quadrature weights w_i^r.
  real(R64), allocatable :: Grid_Spherical_radialWeights(:)

  !> Polar angle quadrature nodes theta_j.
  real(R64), allocatable :: Grid_Spherical_thetaPoints(:)

  !> Polar angle quadrature weights w_j^theta.
  real(R64), allocatable :: Grid_Spherical_thetaWeights(:)

  !> Azimuthal angle quadrature nodes phi_k.
  real(R64), allocatable :: Grid_Spherical_phiPoints(:)

  !> Azimuthal angle quadrature weights w_k^phi.
  real(R64), allocatable :: Grid_Spherical_phiWeights(:)

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
