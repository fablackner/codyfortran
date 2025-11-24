! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Lattice-domain interaction registration.
!>
!> This module is responsible for fabricating interaction implementations that
!> operate on discrete lattice grids (e.g., tight-binding models). It wires the
!> concrete lattice back-end into the public pointers exposed by `M_SysInteraction`.
module M_SysInteraction_Lattice
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Bind lattice-specific interaction procedures according to the selected
    !> lattice model and user parameters (read during global fabrication).
    module subroutine SysInteraction_Lattice_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
