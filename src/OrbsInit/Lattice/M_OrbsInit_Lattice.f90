! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Lattice-based orbital initialization backend.
!>
!> Provides an initialization function interface for orbitals defined on
!> discrete lattice indices (ix, iy, iz). The fabricate routine in this module
!> connects lattice-specific implementations to the generic pointers in
!> `M_OrbsInit`.
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
