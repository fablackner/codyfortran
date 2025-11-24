! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Harmonic potential on lattice grids.
!>
!> Implements a separable 3D harmonic trapping potential on a lattice grid,
!> centered at `(positionX, positionY, positionZ)` with angular frequencies
!> `(omegaX, omegaY, omegaZ)`. The concrete procedures are wired by
!> `SysPotential_Lattice_Harmonic_Fabricate`.
module M_SysPotential_Lattice_Harmonic
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
    module subroutine SysPotential_Lattice_Harmonic_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Harmonic angular frequency along X.
  real(R64) :: SysPotential_Lattice_Harmonic_omegaX
  !> Harmonic angular frequency along Y.
  real(R64) :: SysPotential_Lattice_Harmonic_omegaY
  !> Harmonic angular frequency along Z.
  real(R64) :: SysPotential_Lattice_Harmonic_omegaZ

  !> Trap center position along X.
  real(R64) :: SysPotential_Lattice_Harmonic_positionX
  !> Trap center position along Y.
  real(R64) :: SysPotential_Lattice_Harmonic_positionY
  !> Trap center position along Z.
  real(R64) :: SysPotential_Lattice_Harmonic_positionZ

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
