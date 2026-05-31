! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Seesaw (time-dependent linear tilt) potential on lattice grids.
!>
!> Implements a separable 3D linear potential whose slope oscillates in time:
!>   V(r,t) = sₓ(t)(x-cₓ) + sᵧ(t)(y-cᵧ) + s_z(t)(z-c_z)
!>
!> where the slope along each axis varies sinusoidally:
!>   s(t) = (s_max + s_min)/2 + (s_max - s_min)/2 × sin(2πft)
!>
!> This potential is useful for studying transport phenomena, Bloch oscillations,
!> and non-equilibrium dynamics in optical lattice systems.
module M_SysPotential_Lattice_Seesaw
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Parse parameters and wire the lattice seesaw implementation.
    !>
    !> Reads the configuration (typically from JSON), initializes the module
    !> data (slope bounds and frequencies), and connects the implementation into
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

  !> Oscillation frequency of the seesaw potential along X (Hz).
  real(R64) :: SysPotential_Lattice_Seesaw_frequencyX
  !> Oscillation frequency of the seesaw potential along Y (Hz).
  real(R64) :: SysPotential_Lattice_Seesaw_frequencyY
  !> Oscillation frequency of the seesaw potential along Z (Hz).
  real(R64) :: SysPotential_Lattice_Seesaw_frequencyZ

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
