! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Singlet Valdemoro reconstruction scheme.
!>
!> @details
!> Implements the Valdemoro approximation for spin-singlet states. Works with
!> the α-β block of the 2-RDM and reconstructs the corresponding rdm3aab block.
!>
!> **JSON path**: `reconstruction.singlet.valdemoro`
!>
!> **Options**:
!> - `ccQ` (bool, default=false): Apply contraction-consistent trace constraints
!>   to enforce that partial traces of D³ match the known D² exactly.
!>
!> @see M_Utils_RdmEnforceTraceM3aab for trace constraint enforcement
module M_Reconstruction_Singlet_Valdemoro
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Bind the singlet Valdemoro implementation at runtime.
    !>
    !> Wires the singlet-specific procedure pointers (e.g., alpha–beta block
    !> builders) to Valdemoro implementations and may read/apply variant options
    !> as required by the codebase’s configuration system.
    module subroutine Reconstruction_Singlet_Valdemoro_Fabricate
    end subroutine
  end interface

end module
