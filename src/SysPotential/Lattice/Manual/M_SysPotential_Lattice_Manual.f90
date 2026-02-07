! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Manual potential on lattice grids.
!>
!> Behaviour:
!> - Time-independent: V(i) = fixed manual values drawn once from the chosen
!>   distribution during setup and cached for all subsequent calls
module M_SysPotential_Lattice_Manual
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Parse parameters and wire the lattice manual implementation.
    !>
    !> Reads the configuration (typically from JSON), initializes the module
    !> data (distribution parameters), and connects the implementation into
    !> the lattice backend and grid-agnostic pointers.
    module subroutine SysPotential_Lattice_Manual_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> List of sites to apply values to. If empty, apply to all sites in order.
  integer(I32), allocatable :: SysPotential_Lattice_Manual_sites(:)
  !> List of manual values to apply. Must match the size of sites if provided, or the grid size otherwise.
  real(R64), allocatable :: SysPotential_Lattice_Manual_values(:)
  !> Cached site values generated in Setup
  real(R64), allocatable :: SysPotential_Lattice_Manual_manualValues(:)

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
