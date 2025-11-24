! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Regular 3D lattice grid parameters and runtime wiring.
!>
!> M_Grid_Lattice holds the shape and boundary configuration for a regular
!> three-dimensional lattice along with an optional index map that encodes the
!> linearization used by higher-level code. A corresponding Fabricate routine
!> sets sizes, periodicity flags, and constructs any mapping arrays used by
!> stencil or hopping operators.
module M_Grid_Lattice
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Instantiate and configure a 3D lattice grid back-end.
    !>
    !> Wires lattice-specific procedures (e.g., neighbor traversal, hopping
    !> patterns) and initializes module data below based on input. May branch
    !> to specialized implementations (nearest-neighbor, longer-range, …).
    module subroutine Grid_Lattice_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Mapping from 3D site coordinates to a flattened linear index.
  !>
  !> Grid_Lattice_code(i1,i2,i3) gives the 1-based linear index of the site at
  !> (i1,i2,i3). The exact linearization order is defined by the active
  !> back-end but is guaranteed to be consistent across the lattice modules.
  integer(I32), allocatable :: Grid_Lattice_code(:, :, :)

  !> Number of grid points along the x-dimension of the lattice.
  integer(I32)          :: Grid_Lattice_xSize

  !> Number of grid points along the y-dimension of the lattice.
  integer(I32)          :: Grid_Lattice_ySize

  !> Number of grid points along the z-dimension of the lattice.
  integer(I32)          :: Grid_Lattice_zSize

  !> Periodic boundary condition flag in the x-direction.
  !> When true, sites wrap between the first and last x-layers.
  logical               :: Grid_Lattice_xPeriodicQ

  !> Periodic boundary condition flag in the y-direction.
  !> When true, sites wrap between the first and last y-layers.
  logical               :: Grid_Lattice_yPeriodicQ

  !> Periodic boundary condition flag in the z-direction.
  !> When true, sites wrap between the first and last z-layers.
  logical               :: Grid_Lattice_zPeriodicQ

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
