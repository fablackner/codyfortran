! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief On-site (Hubbard-like) interaction model for lattice grids.
!>
!> @details Implements a local density-density interaction of the form:
!>    Ŵ = U Σᵢ nᵢ↑ nᵢ↓
!> where U is the coupling strength and nᵢσ is the occupation of spin σ at site i.
!>
!> This is the standard Hubbard interaction used in condensed matter simulations
!> of strongly correlated electrons on tight-binding lattices.
!>
!> **JSON configuration:**
!> ```json
!> "sysInteraction": {
!>   "lattice": {
!>     "onSite": {
!>       "strength": 1.0,
!>       "stdImpl": {}
!>     }
!>   }
!> }
!> ```
!>
!> **Parameters:**
!>   - `strength`: Coupling constant U (positive = repulsive, negative = attractive)
module M_SysInteraction_Lattice_OnSite
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Register the on-site lattice interaction with the core system. Reads and
    !> validates parameters, then wires the concrete implementation.
    module subroutine SysInteraction_Lattice_OnSite_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Coupling constant of the local on-site interaction. Positive values
  !> indicate repulsion, negative values attraction (convention depends on
  !> the chosen Hamiltonian).
  real(R64) :: SysInteraction_Lattice_OnSite_Strength

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
