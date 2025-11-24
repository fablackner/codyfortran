! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Nearest-neighbor lattice kinetic operator on a linear lattice.
!>
!> Implements a tight-binding style kinetic term with direction-dependent
!> hopping amplitudes. This module provides configuration storage and a
!> fabrication entry point; concrete multiply/apply routines are assigned by
!> the higher-level lattice fabric.
module M_SysKinetic_Lattice_NearestNeighbor
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Assign the nearest-neighbor lattice backend and read configuration.
    !>
    !> Expected configuration keys (names may vary with the calling fabric):
    !> - hoppX, hoppY, hoppZ: real hopping amplitudes per lattice direction.
    !> The number of active directions depends on the problem dimensionality.
    module subroutine SysKinetic_Lattice_NearestNeighbor_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================
  !> Hopping amplitude along x-direction (default: 1.0)
  real(R64) :: SysKinetic_Lattice_NearestNeighbor_hoppX = 1.0_R64
  !> Hopping amplitude along y-direction (default: 1.0)
  real(R64) :: SysKinetic_Lattice_NearestNeighbor_hoppY = 1.0_R64
  !> Hopping amplitude along z-direction (default: 1.0)
  real(R64) :: SysKinetic_Lattice_NearestNeighbor_hoppZ = 1.0_R64

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
