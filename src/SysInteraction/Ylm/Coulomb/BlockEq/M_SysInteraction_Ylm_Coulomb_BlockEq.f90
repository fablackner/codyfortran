! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Block-factorized FEDVR solver with precomputed Schur complement.
!>
!> @details Exploits the block-tridiagonal structure of the FEDVR discretization
!> to achieve efficient repeated solves. During setup, precomputes:
!>   - LU factorizations for each element block
!>   - Unit response vectors for interface coupling
!>   - Schur complement structure for interface system
!>
!> At solve time, only a small interface system needs to be solved, with the
!> interior recovered via back-substitution using precomputed responses.
!>
!> **Complexity:**
!>   - Setup: O(nE × nLoc³) + O(lmax × nE × nLoc³) per l
!>   - Solve: O(nE² + nE × nLoc) per (l,m) channel
!>
!> **Use case:** Many repeated Poisson solves (time propagation), large grids.
module M_SysInteraction_Ylm_Coulomb_BlockEq
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Bind the block-equation Coulomb solver to the core Ylm interaction hooks.
    module subroutine SysInteraction_Ylm_Coulomb_BlockEq_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
