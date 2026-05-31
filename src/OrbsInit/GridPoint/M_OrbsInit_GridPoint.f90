! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Grid-point (Kronecker δ) orbital initialization backend.
!>
!> @details
!> Initializes orbitals as Kronecker δ functions on a generic discrete grid:
!>
!>   φ_i(j) = δ_{ij}
!>
!> Each orbital is localized at exactly one grid point, with the orbital index
!> directly mapping to the grid point index.
!>
!> This backend is useful for:
!> - Debugging and testing grid infrastructure
!> - Building custom localized basis sets
!> - Direct grid-point representation without physical interpretation
!>
!> JSON Configuration
!> ------------------
!>   {"orbsInit": {"gridPoint": {}}}
!>
!> @note Unlike Lattice/OnSite, this backend does not impose any spatial structure.
module M_OrbsInit_GridPoint
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Select and wire the grid-point initialization implementation.
    !> Expected to parse relevant configuration (e.g., grid size/spacing,
    !> masks, seed functions) and assign the `M_OrbsInit` procedure pointers.
    module subroutine OrbsInit_GridPoint_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
