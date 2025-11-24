! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> 2D polar grid with constant radial spacing.
!>
!> Provides the Fabricate entry point for a uniform-in-r polar grid. The
!> Fabricate routine initializes (r,phi) nodes and weights and wires any
!> polar-constant specific operators.
module M_Grid_Polar_Const
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Instantiate and configure a constant-r 2D polar grid back-end.
    module subroutine Grid_Polar_Const_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
