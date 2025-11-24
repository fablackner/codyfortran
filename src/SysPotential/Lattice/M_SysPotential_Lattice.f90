! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Lattice-grid external potential backend.
!>
!> This module exposes the lattice-specific setup hook and a factory that wires
!> lattice external potential implementations (e.g., harmonic traps) at runtime.
!> All concrete algorithms are selected and connected by
!> `SysPotential_Lattice_Fabricate` based on the JSON configuration.
module M_SysPotential_Lattice
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Lattice factory: parse config and wire concrete implementations.
    !>
    !> Selects the lattice potential model and connects its procedures to the
    !> pointers exported by the grid-agnostic interface. It also parses and
    !> stores lattice-specific parameters as needed by the chosen model.
    module subroutine SysPotential_Lattice_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> Lattice-specific setup hook (allocations, masks, caches for lattice grids).
  procedure(I_SysPotential_Lattice_Setup), pointer :: SysPotential_Setup => NoOpProcedures_Setup
  abstract interface
    !> Initialize the lattice potential backend.
    subroutine I_SysPotential_Lattice_Setup
    end subroutine
  end interface

end module
