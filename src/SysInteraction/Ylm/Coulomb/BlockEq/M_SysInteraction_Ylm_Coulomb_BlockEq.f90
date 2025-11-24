! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Coulomb solver variant: block-equation method in Ylm representation.
!>
!> Provides a memory-efficient scheme that solves radial equations in blocks
!> per (l, m) channel. Useful trade-off between complexity and performance.
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
