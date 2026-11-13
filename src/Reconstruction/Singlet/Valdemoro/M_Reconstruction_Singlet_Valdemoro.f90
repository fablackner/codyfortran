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
!> **Options** (read by the family, S_Reconstruction_Singlet):
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
    !> Binds the closure Reconstruction_Singlet_FillRdm3aab.
    module subroutine Reconstruction_Singlet_Valdemoro_Fabricate
    end subroutine
  end interface

end module
