! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Ground-state solver interface module.
!>
!> @details
!> Provides the top-level contract for self-consistent-field (SCF) ground-state
!> solvers. This module does not contain numerical algorithms itself. Instead,
!> it declares procedure pointers that are bound at runtime by a corresponding
!> `GroundSolver_Fabricate` routine (provided by a backend/implementation).
!>
!> The bound procedures implement a single SCF/relaxation iteration that drives
!> an input quantum state towards a stationary (ground) state. Backends may use
!> different physical models or discretizations (e.g., SCF Hartree–Fock, Ylm
!> radial grids), but they all conform to the same small interface declared here.
!>
!> ## Architecture
!>
!> The module follows the CodyFortranRDM fabrication pattern:
!>   - **Interface Module** (`M_GroundSolver.f90`): Declares procedure pointers
!>   - **Dispatch Submodule** (`S_GroundSolver.f90`): Routes JSON to backends
!>   - **Implementation Submodules**: Provide concrete algorithms
!>
!> ## Usage
!>
!> 1. Call `GroundSolver_Fabricate()` to bind JSON-selected backend
!> 2. Call `GroundSolver_Setup()` to allocate backend resources
!> 3. Loop: Call `GroundSolver_Approach(state, time)` until convergence
!>
!> ## Available Backends
!>
!> - **scf.stdImpl**: Standard SCF Hartree–Fock (general grids)
!> - **scf.ylmOpt**: SCF optimized for spherical-harmonics (per-l channels)
!>
!> @see M_GroundSolver_Scf, M_GroundSolver_Mcscf, M_DiagonalizerList
module M_GroundSolver
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface GroundSolver_Fabricate
    !> @brief Bind concrete implementations to the module's procedure pointers.
    !>
    !> @details
    !> Reads the JSON configuration key `groundSolver.*` and dispatches to the
    !> appropriate backend's `_Fabricate` routine. This assigns the procedure
    !> pointers `GroundSolver_Setup` and `GroundSolver_Approach` to backend-
    !> specific procedures (for example an SCF or MCSCF variant).
    !>
    !> No allocations or state creation occur here—only pointer assignment.
    !>
    !> @pre JSON configuration loaded via `Json_LoadJsonFile`
    !> @post `GroundSolver_Setup` and `GroundSolver_Approach` are bound
    module subroutine GroundSolver_Fabricate()
    end subroutine
  end interface

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> @brief Pointer to the setup routine for the ground-state solver backend.
  !>
  !> Initialized to a no-op; bound to a real implementation by `_Fabricate`.
  procedure(I_GroundSolver_Setup), pointer :: GroundSolver_Setup => NoOpProcedures_Setup
  abstract interface
    !> @brief Initialize the ground-state solver backend.
    !>
    !> @details
    !> Allocates working buffers (potentials, orbital temporaries), reads any
    !> backend-specific JSON configuration, and validates input dimensions
    !> against the configured grid. This routine performs no SCF iterations.
    !>
    !> @pre All required modules (Grid, SysKinetic, etc.) must be set up.
    !> @post Backend is ready for calls to `GroundSolver_Approach`.
    subroutine I_GroundSolver_Setup
    end subroutine
  end interface

  !> @brief Pointer to one ground-state approach/relaxation step.
  !>
  !> This is the main workhorse—call repeatedly until energy convergence.
  procedure(I_GroundSolver_Approach), pointer :: GroundSolver_Approach
  abstract interface
    !> @brief Perform a single SCF approach/relaxation iteration on the state.
    !>
    !> @details
    !> Implements one iteration of the self-consistent-field cycle:
    !>   1. Build mean-field potentials (Hartree + external) from current orbitals
    !>   2. Diagonalize the effective Fock operator
    !>   3. Gauge-align the eigenvectors with the old orbitals
    !>      (`Orbs_AlignOnReference`), then mix via `Mixing_Mix`
    !>   4. Orthonormalize the mixed orbitals
    !>
    !> The caller is responsible for the convergence loop (checking energy change).
    !>
    !> @param[inout] state   Packed state vector (orbitals). Updated in place.
    !> @param[in]    time    Backend-defined parameter. Typically 0 for ground state;
    !>                       some backends may use it for time-dependent potentials.
    !>
    !> @note For spin-restricted calculations, the spin-up block is copied to
    !>       spin-down after mixing, ensuring identical spatial orbitals.
    subroutine I_GroundSolver_Approach(state, time)
      import :: R64
      !> Packed state vector containing orbital coefficients. Updated in place.
      complex(R64), intent(inout), contiguous, target :: state(:)
      !> Backend-defined step parameter (e.g., physical or imaginary time).
      real(R64), intent(in) :: time
    end subroutine
  end interface

end module
