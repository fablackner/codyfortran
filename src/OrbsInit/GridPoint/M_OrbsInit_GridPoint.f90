! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Grid-point orbital initialization backend.
!>
!> This module fabricates a backend that initializes orbitals defined directly
!> on a discrete grid of points (no special structure assumed). Its fabricate
!> routine configures and connects the generic `M_OrbsInit` procedure pointers
!> to grid-point specific implementations.
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
