! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief File-based orbital initialization backend.
!>
!> @details
!> Loads pre-computed orbital data from binary files rather than constructing
!> orbitals analytically. Useful for:
!> - Restarting calculations from a previous run
!> - Loading orbitals computed by external codes
!> - Using optimized/converged orbitals as initial guesses
!>
!> File Naming Convention
!> ----------------------
!> Orbitals are loaded from files named: orb{bt}_{ind}.in
!> where bt is the 2-digit body type and ind is the 2-digit orbital index.
!> Example: orb01_03.in contains orbital 3 of body type 1.
!>
!> File Format
!> -----------
!> Binary format compatible with M_Utils_DataStorage (LoadData routine).
!> Files must contain complex(R64) arrays of length nGrid.
!>
!> JSON Configuration
!> ------------------
!>   {"orbsInit": {"load": {}}}
!>
!> @note Files must exist in the working directory at runtime.
module M_OrbsInit_Load
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Select and wire the load-from-disk initialization implementation.
    !> Expected to parse data source specification (paths, formats, indexing) and
    !> to assign the `M_OrbsInit` procedure pointers appropriately.
    module subroutine OrbsInit_Load_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
