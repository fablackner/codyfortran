! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Reference implementation for the lattice on-site interaction.
!>
!> Provides a straightforward, portable realization of the on-site model. More
!> optimized or specialized variants can coexist alongside this default.
module M_SysInteraction_Lattice_OnSite_StdImpl
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Bind the default on-site lattice implementation to the core interaction
    !> procedure pointers.
    module subroutine SysInteraction_Lattice_OnSite_StdImpl_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
