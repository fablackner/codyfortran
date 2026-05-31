! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Central API for orbital initialization in the CodyFortranRDM framework.
!>
!> @details
!> This module exposes procedure pointers that higher-level components use to
!> initialize single orbitals or full orbital sets on a chosen grid/representation
!> (linear, lattice, grid-point, Ylm, file-load, ...). The concrete procedures
!> are selected and wired at runtime by `OrbsInit_Fabricate`, typically based on
!> user configuration (e.g., JSON input) and problem geometry.
!>
!> Architecture
!> ------------
!> The module follows the Interface/Implementation split pattern used throughout
!> CodyFortranRDM. This file (M_OrbsInit.f90) defines the public interface via
!> abstract interfaces and procedure pointers. The companion submodule
!> (S_OrbsInit.f90) implements the branching logic in `OrbsInit_Fabricate`.
!>
!> Available Backends
!> ------------------
!> - **Linear**    : 1D continuous grids (harmonic oscillator eigenstates, etc.)
!> - **Lattice**   : Discrete 3D lattice sites (on-site δ-orbitals for Hubbard models)
!> - **GridPoint** : Generic discrete grid (Kronecker δ basis)
!> - **Ylm**       : Spherical harmonics basis (hydrogen-like radial functions)
!> - **Load**      : File-based initialization from binary orbital data
!>
!> Contract (Public API)
!> ---------------------
!> - `OrbsInit_Fabricate()` : Reads JSON config, selects backend, binds pointers.
!> - `OrbsInit_Setup()`     : Prepares module-local state (parameters, buffers).
!> - `OrbsInit_Initialize(orbs)` : Fills a 2D array of orbitals (nGrid × nOrbs).
!> - `OrbsInit_InitializeOrb(orb, ind, bt_)` : Initializes a single orbital vector.
!>
!> Typical Usage Sequence
!> ----------------------
!> 1. call OrbsInit_Fabricate   ! JSON-driven backend selection
!> 2. call OrbsInit_Setup       ! Allocate buffers, precompute constants
!> 3. call OrbsInit_Initialize(orbs)  ! Fill orbital array
!>
!> JSON Configuration Examples
!> ---------------------------
!> Linear/Harmonic:    {"orbsInit": {"linear": {"harmonic": {"omega": 1.0}}}}
!> Lattice/OnSite:     {"orbsInit": {"lattice": {"onSite": {}}}}
!> Ylm/HydrogenLike:   {"orbsInit": {"ylm": {"hydrogenLike": {"charge": 1.0,
!>                                    "n": [1,2], "l": [0,0], "m": [0,0]}}}}
!> Load from file:     {"orbsInit": {"load": {}}}
!>
!> Notes
!> -----
!> - All arrays are assumed contiguous and column-major (Fortran default).
!> - Normalization is performed via Grid_InnerProduct after raw initialization.
!> - Phase conventions are backend-specific but consistent within a single run.
!> - Body types (bt_) enable species-dependent initialization (e.g., spin-up/down).
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
