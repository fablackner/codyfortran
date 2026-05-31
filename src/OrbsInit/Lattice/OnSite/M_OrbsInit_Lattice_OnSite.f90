! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief On-site δ-orbital initializer for lattice models.
!>
!> @details
!> Implements orbital initialization as Kronecker δ functions localized at
!> individual lattice sites. Each orbital is nonzero only at one site:
!>
!>   φ_i(j) = δ_{ij}
!>
!> The orbital index maps to site coordinates via row-major ordering:
!>   index → (ix, iy, iz) where x varies fastest, then y, then z.
!>
!> This is the standard basis for Hubbard model calculations where orbitals
!> represent localized atomic states on each lattice site.
!>
!> JSON Configuration
!> ------------------
!>   "orbsInit": {
!>     "lattice": {
!>       "onSite": { }  // no additional parameters needed
!>     }
!>   }
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
