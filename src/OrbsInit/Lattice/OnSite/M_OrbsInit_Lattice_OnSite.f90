! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
module M_OrbsInit_Lattice_OnSite
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Fabricate an on-site lattice initializer.
    !> Specializes the lattice backend to on-site shapes (e.g., localized orbitals
    !> centered on lattice points). Assigns any necessary callouts in the lattice
    !> module and caches on-site parameters from configuration.
    module subroutine OrbsInit_Lattice_OnSite_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
