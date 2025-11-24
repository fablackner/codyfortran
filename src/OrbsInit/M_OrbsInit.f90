! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Central API for orbital initialization.
!>
!> This module exposes procedure pointers that higher-level components use to
!> initialize single orbitals or full orbital sets on a chosen grid/representation
!> (linear, lattice, grid-point, Ylm, file-load, ...). The concrete procedures
!> are selected and wired at runtime by `OrbsInit_Fabricate`, typically based on
!> user configuration (e.g., JSON input) and problem geometry.
!>
!> Contract
!> - `OrbsInit_Setup()` prepares module-local state (parameters, masks, buffers).
!> - `OrbsInit_Initialize(orbs)` fills a 2D array of orbitals in-place.
!> - `OrbsInit_InitializeOrb(orb, ind, bt_)` initializes one orbital vector by index
!>   and optional body-type/species.
!>
!> Notes
!> - All arrays are assumed contiguous and column-major (Fortran default).
!> - Normalization and phase conventions are backend-specific but should be
!>   consistent across a single run once fabricated.
module M_OrbsInit
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Selects and assigns the concrete initialization backend.
    !>
    !> This routine sets the procedure pointers exported by this module based on
    !> runtime input (e.g., JSON configuration) and available submodules. It may
    !> also parse and cache parameters that are needed during setup/initialization.
    module subroutine OrbsInit_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> Pointer to the setup procedure for initializing the OrbsInit system.
  procedure(I_OrbsInit_Setup), pointer :: OrbsInit_Setup => NoOpProcedures_Setup
  abstract interface
    !> Prepare backend-specific state required before any orbital is created.
    !> Typical tasks: allocate scratch buffers, precompute masks/basis functions,
    !> validate configuration, and preload constants needed by initialization.
    subroutine I_OrbsInit_Setup
    end subroutine
  end interface

  !> Pointer to the procedure for initializing orbitals.
  procedure(I_OrbsInit_Initialize), pointer :: OrbsInit_Initialize
  abstract interface
    !> Initialize a set of orbitals in one call.
    subroutine I_OrbsInit_Initialize(orbs)
      import :: I32, R64
      !> Output 2D array of orbitals with shape (nGrid, nOrbitals).
      !> Each column corresponds to one orbital sampled on the active grid.
      complex(R64), intent(out), contiguous :: orbs(:, :)
    end subroutine
  end interface

  !> Pointer to the procedure for initializing a single orbital of specific body type.
  procedure(I_OrbsInit_InitializeOrb), pointer :: OrbsInit_InitializeOrb
  abstract interface
    !> Initialize one orbital vector by index and optional body type/species.
    !> The active backend defines the basis/grid and any normalization.
    subroutine I_OrbsInit_InitializeOrb(orb, ind, bt_)
      import :: I32, R64
      !> Output vector containing the initialized orbital sampled on the grid.
      complex(R64), intent(out), contiguous :: orb(:)
      !> Orbital index within the body type/species.
      integer(I32), intent(in)              :: ind
      !> Optional body type/species index (e.g., electron spin species, bath id).
      integer(I32), intent(in), optional    :: bt_
    end subroutine
  end interface

end module
