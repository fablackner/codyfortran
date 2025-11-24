! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> 1D linear grid with constant spacing (uniform Cartesian x-grid).
!>
!> This module hosts the Fabricate entry for constant-spacing 1D grids and
!> wires the corresponding operators (e.g., finite differences, FFT-based
!> kinetic operators) provided by specialized submodules.
module M_Grid_Linear_Const
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Instantiate and configure a constant-spacing 1D linear grid.
    !>
    !> Initializes coordinates and weights for a uniform grid on [xmin,xmax]
    !> with a chosen number of points and associates constant-grid operators.
    module subroutine Grid_Linear_Const_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
