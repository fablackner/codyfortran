! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Orbital initialization backends for the prolate-spheroidal grid.
!>
!> @details
!> Dispatches to the concrete initializer:
!> - `lcao`: gerade/ungerade combinations of hydrogen-like atomic orbitals
!>   centered on the two nuclei (LCAO guess for diatomic molecules).
module M_OrbsInit_Prolate
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Fabricate a prolate-grid orbital initializer.
    module subroutine OrbsInit_Prolate_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
