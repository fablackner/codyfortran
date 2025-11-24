! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Default/reference Coulomb solver in Ylm representation.
!>
!> A robust, easy-to-understand baseline implementation used when no specialized
!> method is requested. Prioritizes clarity over peak performance.
module M_SysInteraction_Ylm_Coulomb_StdImpl
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Bind the default Ylm Coulomb solver to the core interaction API.
    module subroutine SysInteraction_Ylm_Coulomb_StdImpl_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
