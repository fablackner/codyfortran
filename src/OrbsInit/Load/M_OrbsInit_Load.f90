! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> File-backed orbital initialization backend.
!>
!> Provides an initialization mode that loads orbitals from external sources
!> (e.g., files on disk) rather than constructing them analytically. The
!> fabricate routine configures file paths/format and wires the generic
!> `M_OrbsInit` pointers to load-based implementations.
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
