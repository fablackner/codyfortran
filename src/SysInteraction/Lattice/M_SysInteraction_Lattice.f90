! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Lattice-domain interaction interface.
!>
!> @details This module registers interaction implementations that operate on
!> discrete lattice grids, typical for tight-binding or Hubbard-type models.
!>
!> **Physics context:**
!> On a lattice, particle-particle interactions are typically local (on-site)
!> or short-ranged. The canonical example is the Hubbard-U term:
!>    Ŵ = U Σᵢ n↑ᵢ n↓ᵢ
!> which penalizes double occupancy at each site i.
!>
!> **Available implementations:**
!>   - `OnSite`: Local density-density interaction (Hubbard U)
!>
!> The lattice back-end provides `FillInteractionSrc` and
!> `MultiplyWithInteractionPotential` common to all lattice models; the
!> specific physics (e.g., strength U) is set by the leaf implementation.
module M_SysInteraction_Lattice
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Bind lattice-specific interaction procedures according to the selected
    !> lattice model and user parameters (read during global fabrication).
    module subroutine SysInteraction_Lattice_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
