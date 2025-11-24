! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> 2D Cartesian grid with constant spacing in x and y.
!>
!> Provides the Fabricate entry point for a uniform Cartesian grid and wires
!> constant-spacing specific operators (e.g., finite differences, FFT) as
!> supplied by specialized submodules.
module M_Grid_Square_Const
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Instantiate and configure a constant-spacing 2D Cartesian grid.
    module subroutine Grid_Square_Const_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
