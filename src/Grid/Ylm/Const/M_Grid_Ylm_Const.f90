! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Ylm radial grid with constant spacing.
!>
!> Provides the Fabricate entry for a Ylm grid whose radial coordinate uses a
!> uniform partition. The routine initializes radial nodes/weights and wires
!> Ylm-specific operators compatible with the constant radial layout.
module M_Grid_Ylm_Const
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Instantiate and configure a constant-r Ylm radial grid back-end.
    module subroutine Grid_Ylm_Const_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
