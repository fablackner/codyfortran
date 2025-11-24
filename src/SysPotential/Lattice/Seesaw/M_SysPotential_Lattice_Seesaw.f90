! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Seesaw potential on lattice grids.
!>
!> Implements a separable 3D linear potential on a lattice grid.
!> The slope of the potential oscillates in time between a minimum and maximum
!> value with a given frequency. The concrete procedures are wired by
!> `SysPotential_Lattice_Seesaw_Fabricate`.
module M_SysPotential_Lattice_Seesaw
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Parse parameters and wire the lattice harmonic implementation.
    !>
    !> Reads the configuration (typically from JSON), initializes the module
    !> data (trap center and frequencies), and connects the implementation into
    !> the lattice backend and grid-agnostic pointers.
    module subroutine SysPotential_Lattice_Seesaw_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Minimum slope of the seesaw potential along X.
  real(R64) :: SysPotential_Lattice_Seesaw_slopeMinX
  !> Minimum slope of the seesaw potential along Y.
  real(R64) :: SysPotential_Lattice_Seesaw_slopeMinY
  !> Minimum slope of the seesaw potential along Z.
  real(R64) :: SysPotential_Lattice_Seesaw_slopeMinZ

  !> Maximum slope of the seesaw potential along X.
  real(R64) :: SysPotential_Lattice_Seesaw_slopeMaxX
  !> Maximum slope of the seesaw potential along Y.
  real(R64) :: SysPotential_Lattice_Seesaw_slopeMaxY
  !> Maximum slope of the seesaw potential along Z.
  real(R64) :: SysPotential_Lattice_Seesaw_slopeMaxZ

  !> Frequency of the seesaw potential along X.
  real(R64) :: SysPotential_Lattice_Seesaw_frequencyX
  !> Frequency of the seesaw potential along Y.
  real(R64) :: SysPotential_Lattice_Seesaw_frequencyY
  !> Frequency of the seesaw potential along Z.
  real(R64) :: SysPotential_Lattice_Seesaw_frequencyZ

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
