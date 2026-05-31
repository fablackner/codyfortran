! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Reference implementation for on-site lattice interaction.
!>
!> @details Provides the straightforward realization: V(i) = U * ρ(i).
!> The potential at each site is simply the coupling strength times the
!> local density, implementing the Hubbard-U term directly.
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
