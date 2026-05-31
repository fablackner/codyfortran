! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Lattice-based (3D discrete) orbital initialization backend.
!>
!> @details
!> Provides an initialization function interface for orbitals defined on
!> discrete lattice sites (ix, iy, iz). Designed for tight-binding and Hubbard
!> model simulations where orbitals are localized on lattice points.
!>
!> Available Sub-Backends
!> ----------------------
!> - **OnSite** : Kronecker δ orbitals localized at individual lattice sites.
!>                Orbital index maps to site via row-major ordering (x fastest).
!>
!> Procedure Pointer Contract
!> --------------------------
!> `OrbsInit_Lattice_InitFunction(ix, iy, iz, index, bt_)` returns the amplitude
!> at site (ix, iy, iz) for orbital `index`. Normalization applied after sampling.
!>
!> JSON Configuration
!> ------------------
!>   {"orbsInit": {"lattice": {"onSite": {}}}}
!>
!> @note Lattice dimensions come from Grid_Lattice_xSize/ySize/zSize.
!> @see M_OrbsInit_Lattice_OnSite for the on-site δ-orbital implementation.
module M_OrbsInit_Lattice
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Select and wire the lattice initialization implementation.
    !> Reads any lattice-related configuration (dimensions, boundary conditions,
    !> sublattices, etc.) and assigns the `M_OrbsInit` procedure pointers.
    module subroutine OrbsInit_Lattice_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> Lattice-site initializer used by lattice backends.
  procedure(I_OrbsInit_Lattice_InitFunction), pointer :: OrbsInit_Lattice_InitFunction
  abstract interface
    !> Compute the real-valued amplitude for a lattice orbital at site (ix,iy,iz).
    function I_OrbsInit_Lattice_InitFunction(ix, iy, iz, index, bt_) result(res)
      import :: I32, R64
      !> Resulting orbital amplitude/value at the site.
      real(R64)                 :: res
      !> Lattice index along x (1-based unless stated otherwise by the backend).
      integer(I32), intent(in) :: ix
      !> Lattice index along y.
      integer(I32), intent(in) :: iy
      !> Lattice index along z.
      integer(I32), intent(in) :: iz
      !> Orbital index within the current species/body type.
      integer(I32), intent(in) :: index
      !> Optional body type/species/bath index.
      integer(I32), intent(in), optional :: bt_
    end function
  end interface

end module
