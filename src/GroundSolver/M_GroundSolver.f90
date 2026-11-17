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
!> The module follows the CodyFortran fabrication pattern:
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
!> - **mcscf.stdImpl**, **mcscf.ylmOpt**: multiconfiguration SCF
!>
!> ## Orbital preconditioner (all backends)
!>
!> `GroundSolver_PreconditionAction` is bound by `GroundSolver_Fabricate` for
!> every backend. It is the `ApplyPrecOnVec` callback of a preconditioned
!> diagonalizer (lobpcg) for any orbital eigenproblem of the solvers (the SCF
!> Fock operator, the MCSCF RDM-weighted Fock operator):
!>   dOrb = (T̂ + max(−shift, σ_min))⁻¹·orb
!> through the optional kinetic hook `SysKinetic_SolveShifted` (exact on
!> grid.fem), the identity when the kinetic backend does not provide it;
!> σ_min = `GroundSolver_precShiftMin` (JSON `groundSolver.precShiftMin`,
!> default 0.5).
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
    !> For every backend it also reads `groundSolver.precShiftMin` and binds
    !> `GroundSolver_PreconditionAction`.
    !>
    !> @pre JSON configuration loaded via `Json_LoadJsonFile`
    !> @post `GroundSolver_Setup`, `GroundSolver_Approach` and
    !>       `GroundSolver_PreconditionAction` are bound
    module subroutine GroundSolver_Fabricate()
    end subroutine
  end interface

  interface GroundSolver_NoteUnconverged
    !> @brief One-line note if diagonalizer iDiag did not converge.
    !>
    !> @details
    !> Called by every backend after each `DiagonalizerList(iDiag)` solve:
    !> prints a note when its `convergedQ` is false (lobpcg at maxIter or
    !> after a stagnation exit) and returns; the solver continues with the
    !> best approximations the diagonalizer returned (the next Approach
    !> warm-starts from them), so an SCF can still converge.
    module subroutine GroundSolver_NoteUnconverged(iDiag)
      integer(I32), intent(in) :: iDiag
    end subroutine
  end interface

  interface GroundSolver_NumberOfSpinUpOrbs
    !> @brief Number of spin-up (spatial) orbitals in the packed state.
    !>
    !> @details
    !> The ground-state solvers implement spin-restricted physics: one spatial
    !> orbital per doubly occupied shell is relaxed. For `orbs.restrictedQ`
    !> the state stores only that shared spatial set, so all stored orbitals
    !> are spin-up; otherwise the state stores both spin blocks and the
    !> spin-up half comes first.
    module function GroundSolver_NumberOfSpinUpOrbs() result(nUp)
      integer(I32) :: nUp
    end function
  end interface

  interface GroundSolver_AlignOrbitals
    !> @brief Prepare new spin-up orbitals for acceptance.
    !>
    !> @details
    !> Orthonormalizes the new orbitals (e.g. Fock eigenvectors, normalized
    !> in the Euclidean metric) in the grid metric, which alignment and mixing
    !> require, and gauge-aligns them on the old ones (`Orbs_AlignOnReference`:
    !> removes the arbitrary per-vector phase and the rotation within
    !> degenerate shells; with blockOf_ only within blocks of equal declared m).
    module subroutine GroundSolver_AlignOrbitals(orbsRaw, oldOrbs, blockOf_)
      !> New spin-up orbitals (nPoints x nUp), orthonormalized and aligned in place.
      complex(R64), intent(inout), contiguous :: orbsRaw(:, :)
      !> Old spin-up orbitals (nPoints x nUp).
      complex(R64), intent(in), contiguous :: oldOrbs(:, :)
      !> Symmetry block of every orbital (Orbs_blockOfOrb(1:nUp)).
      integer(I32), intent(in), contiguous, optional :: blockOf_(:)
    end subroutine
  end interface

  interface GroundSolver_StoreOrbitals
    !> @brief Accept new spin-up orbitals into the state.
    !>
    !> @details
    !> Mixes them into the spin-up orbitals of the state (`Mixing_Mix`, the
    !> mixTarget "orbitals") or replaces those, copies spin up to spin down
    !> unless the orbitals are stored restricted, and orthonormalizes all
    !> stored orbitals (`Orbs_Orthonormalize`).
    module subroutine GroundSolver_StoreOrbitals(orbsState, newOrbs, mixQ)
      !> The orbital part of the state (nPoints x Orbs_nOrbsInState, flat).
      complex(R64), intent(inout), contiguous, target :: orbsState(:)
      !> New spin-up orbitals (nPoints x nUp), aligned (GroundSolver_AlignOrbitals).
      complex(R64), intent(in), contiguous, target :: newOrbs(:, :)
      !> Mix (true) or replace (false).
      logical, intent(in) :: mixQ
    end subroutine
  end interface

  interface GroundSolver_MixPotential
    !> @brief Accept the next mean-field potential of mixTarget "potential".
    !>
    !> @details
    !> raw is the potential of the current orbitals. The first call (mixed
    !> not allocated) takes it as is: the initial orbitals are not the
    !> eigenvectors of any accepted potential, so their potential is the start
    !> iterate, not a fixed-point step (a zero start costs up to two SCF
    !> steps and puts a pair of no map into the DIIS history). Every later
    !> call mixes raw, the image of mixed, into mixed (`Mixing_Mix`).
    module subroutine GroundSolver_MixPotential(mixed, raw)
      complex(R64), intent(inout), allocatable :: mixed(:)
      complex(R64), intent(in), contiguous :: raw(:)
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Smallest shift σ_min of the orbital preconditioner (T̂ + σ)⁻¹
  !> (`GroundSolver_PreconditionAction`); read by `GroundSolver_Fabricate`
  !> from `groundSolver.precShiftMin` (default 0.5, must be positive).
  real(R64) :: GroundSolver_precShiftMin = 0.5_R64

  !> Options that only some backends implement. The fabricate of a backend
  !> sets the flags of those it implements; GroundSolver_Fabricate resets them
  !> before it dispatches and refuses a configured option afterwards if its
  !> flag is unset:
  !>   - declared orbital m (method.mb.orbBased.mOfOrbs): the backend keeps
  !>     every orbital in its m channel (radial Fock operators per l,
  !>     alignment per m block, or an m-conserving update);
  !>   - Kohn-Sham DFT (sysXc): the backend builds the mean field v_KS - a K.
  logical :: GroundSolver_implementsDeclaredMQ = .false.
  logical :: GroundSolver_implementsKohnShamQ = .false.

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
  procedure(I_GroundSolver_Approach), pointer :: GroundSolver_Approach => null()
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
    !> @note The solvers implement spin-restricted physics. Without
    !>       `orbs.restrictedQ` the state stores both spin blocks and the
    !>       spin-up block is copied to spin-down after mixing; with
    !>       `orbs.restrictedQ` only the shared spatial set is stored and no
    !>       copy is needed (see `GroundSolver_NumberOfSpinUpOrbs`).
    subroutine I_GroundSolver_Approach(state, time)
      import :: R64
      !> Packed state vector containing orbital coefficients. Updated in place.
      complex(R64), intent(inout), contiguous, target :: state(:)
      !> Backend-defined step parameter (e.g., physical or imaginary time).
      real(R64), intent(in) :: time
    end subroutine
  end interface

  !> @brief Pointer to the orbital preconditioner (≈ (F̂ − shift)⁻¹).
  !>
  !> Bound by `GroundSolver_Fabricate` for every backend (null before).
  !> Matches the `ApplyPrecOnVec` callback of DiagonalizerList, so a program
  !> wires it into a preconditioned diagonalizer (lobpcg) of any orbital
  !> eigenproblem (SCF Fock operator, MCSCF orbital step):
  !> `input(i) % ApplyPrecOnVec => GroundSolver_PreconditionAction`
  !> (after `GroundSolver_Fabricate`; the assignment captures the target).
  procedure(I_GroundSolver_PreconditionAction), pointer :: GroundSolver_PreconditionAction => null()
  abstract interface
    !> @brief Apply the orbital preconditioner to a residual.
    !>
    !> @details
    !> dOrb = (T̂ + max(−shift, σ_min))⁻¹·orb if the kinetic backend provides
    !> `SysKinetic_SolveShifted` (body type 1, i.e. the mass of the spin-
    !> restricted orbitals), else dOrb = orb (identity);
    !> σ_min = `GroundSolver_precShiftMin`. For occupied orbitals (Ritz value
    !> θ < 0) the shifted kinetic operator T̂ − θ approximates F̂ − θ in its
    !> stiff high-momentum part, which dominates the spectral radius on fine
    !> grids; it is Hermitian positive definite in the grid metric for every
    !> shift, as LOBPCG requires. Defined for orbital vectors
    !> (Grid_nPoints entries), not for CI vectors.
    !>
    !> @param[out] dOrb   Preconditioned vector
    !> @param[in]  orb    Input vector (residual of a Ritz pair)
    !> @param[in]  shift  Current Ritz value (orbital energy estimate)
    !> @param[in]  time   Time parameter (unused)
    subroutine I_GroundSolver_PreconditionAction(dOrb, orb, shift, time)
      import :: R64
      !> Output: preconditioned vector
      complex(R64), intent(out), contiguous, target :: dOrb(:)
      !> Input: vector to precondition
      complex(R64), intent(in), contiguous, target :: orb(:)
      !> Ritz value the preconditioner adapts to
      real(R64), intent(in) :: shift
      !> Time parameter
      real(R64), intent(in) :: time
    end subroutine
  end interface

end module
