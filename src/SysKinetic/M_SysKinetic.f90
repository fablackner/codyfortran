! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file M_SysKinetic.f90
!> @brief Central facade for kinetic-energy operators T̂ used throughout the codebase.
!>
!> @details
!> This module declares procedure pointers that are bound at runtime by
!> `SysKinetic_Fabricate`, based on the simulation setup (grid type, numerical
!> scheme, masses, etc.). Callers use these pointers without depending on a
!> specific backend implementation.
!>
!> Physics Background
!> ------------------
!> The kinetic energy operator in quantum mechanics is:
!>
!>     T̂ = −(ℏ²/2m) ∇²
!>
!> In atomic units (ℏ = m_e = 1), this simplifies to T̂ = −(1/2m)∇² for a
!> particle of mass m. This module provides pluggable discretizations of ∇²
!> on various spatial grids:
!>
!>   - **Linear**: 1D Cartesian grid using finite-difference or Fourier methods
!>   - **Lattice**: Tight-binding hopping on discrete 3D lattice sites
!>   - **Ylm**: Radial Laplacian in spherical harmonics: −(1/2m)[d²/dr² + (2/r)d/dr − l(l+1)/r²]
!>   - **Prolate**: Two-center prolate-spheroidal Laplacian per azimuthal channel
!>   - **Fem**: 3D Cartesian spectral-element Laplacian (nodal, lumped mass)
!>
!> Contract
!> --------
!> 1. Fabrication selects and assigns concrete procedures before first use.
!> 2. Setup is idempotent and prepares any masks, work arrays, or constants.
!> 3. Kinetic application operates on contiguous complex arrays representing a
!>    single orbital; an optional body-type index `bt_` allows mass scaling.
!> 4. Optional hooks (`SysKinetic_SolveShifted`) stay null unless the backend
!>    provides them (`SysKinetic_Fabricate` nullifies them before it
!>    dispatches); callers test them with `associated(...)`.
!>
!> Initialization Sequence
!> -----------------------
!> @code{.f90}
!> call SysKinetic_Fabricate    ! Reads JSON, binds procedure pointers
!> call SysKinetic_Setup        ! Allocates precomputed data (stencils, k-vectors, etc.)
!> @endcode
!>
!> @see M_SysKinetic_Linear, M_SysKinetic_Lattice, M_SysKinetic_Ylm, M_SysKinetic_Prolate,
!>      M_SysKinetic_Fem
module M_SysKinetic
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> @brief Runtime factory for kinetic operators.
    !>
    !> @details
    !> Reads user configuration (JSON input) and assigns the procedure
    !> pointers exported by this module to concrete implementations:
    !>
    !>   - `sysKinetic.linear`  → Cartesian grid Laplacian (FinDiff/Fourier)
    !>   - `sysKinetic.lattice` → Tight-binding hopping (nearest-neighbor)
    !>   - `sysKinetic.ylm`     → Spherical radial Laplacian (FinDiff/FEDVR)
    !>   - `sysKinetic.prolate` → Prolate-spheroidal Laplacian (FEDVR)
    !>   - `sysKinetic.fem`     → 3D spectral-element Laplacian (binds SolveShifted)
    !>
    !> Must be invoked exactly once per simulation context before calling
    !> any other procedure in this module. Subsequent calls will overwrite
    !> the previously assigned pointers.
    !>
    !> JSON Configuration Example
    !> --------------------------
    !> @code{.json}
    !> {
    !>   "sysKinetic": {
    !>     "linear": {
    !>       "laplacian": {
    !>         "bodyMass": [1.0],
    !>         "finDiff": {}
    !>       }
    !>     }
    !>   }
    !> }
    !> @endcode
    module subroutine SysKinetic_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data (set during Fabricate)
  !=============================================================================

  !> @brief True if the kinetic operator does not depend on time.
  !>
  !> Currently all backends (Linear, Lattice, Ylm) set this to `.true.` since
  !> kinetic operators are typically time-independent. Time-dependent mass
  !> scaling or absorbing boundaries would require setting this to `.false.`.
  logical :: SysKinetic_timeIndependentQ = .false.

  !> @brief True if the kinetic operator is the same for all body types.
  !>
  !> When `.true.`, the `bt_` argument to `MultiplyWithKineticOp` can be omitted
  !> (all particles share the same mass). When `.false.`, different masses apply
  !> per body type, e.g., electrons vs. nuclei or spin-up vs. spin-down.
  logical :: SysKinetic_bodyTypeIndependentQ = .false.

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> @brief Pointer to the setup procedure for initializing the kinetic subsystem.
  !>
  !> @details
  !> Responsibilities:
  !>   - Allocate/capture grid-dependent constants and masks
  !>   - Precompute operator coefficients (e.g., k² arrays for Fourier, FD stencils)
  !>   - Validate dimensional compatibility with the chosen backend
  !>
  !> Must be called after `Grid_Setup` since discretization depends on grid spacing.
  !> Default points to `NoOpProcedures_Setup` (no-op) until `Fabricate` assigns it.
  procedure(I_SysKinetic_Setup), pointer :: SysKinetic_Setup => NoOpProcedures_Setup
  abstract interface
    !> @brief Initializes the kinetic subsystem.
    !>
    !> Prepares mask functions, numerical stencils, and any other runtime data
    !> needed by the selected kinetic implementation.
    subroutine I_SysKinetic_Setup
    end subroutine
  end interface

  !> @brief Pointer to the procedure that applies T̂ψ = −(1/2m)∇²ψ to a single orbital.
  !>
  !> @details
  !> The concrete operator depends on the selected backend:
  !>
  !>   - **Linear**: −(1/2m)∇² via finite-difference or Fourier on Cartesian grid
  !>   - **Lattice**: Tight-binding hopping −t Σ_{⟨i,j⟩} (|i⟩⟨j| + h.c.)
  !>   - **Ylm**: Radial kinetic −(1/2m)[d²/dr² + (2/r)d/dr − l(l+1)/r²] per channel
  !>
  !> Notes:
  !>   - Arrays must be contiguous; shapes and BCs are backend-specific
  !>   - If `bt_` is provided, selects body-type specific mass for 1/(2 m_bt)
  !>   - Currently all backends are time-independent (time argument unused)
  procedure(I_SysKinetic_MultiplyWithKineticOp), pointer :: SysKinetic_MultiplyWithKineticOp
  abstract interface
    !> @brief Applies the kinetic-energy operator: dOrb = T̂ · orb
    !>
    !> @param[out] dOrb  Result of applying T̂ to the orbital (length = nGridPoints)
    !> @param[in]  orb   Input orbital wavefunction ψ(r) on the spatial grid
    !> @param[in]  time  Simulation time (reserved for time-dependent extensions)
    !> @param[in]  bt_   Optional body-type index for mass-scaling (default = 1)
    subroutine I_SysKinetic_MultiplyWithKineticOp(dOrb, orb, time, bt_)
      import :: I32, R64
      !> Output orbital after applying the kinetic operator: dOrb = T̂ · orb
      complex(R64), intent(out), contiguous :: dOrb(:)
      !> Input orbital wavefunction ψ(r)
      complex(R64), intent(in), contiguous :: orb(:)
      !> Current simulation time (unused by current backends)
      real(R64), intent(in) :: time
      !> Body type index (1-based) selecting mass from bodyMass array
      integer(I32), intent(in), optional :: bt_
    end subroutine
  end interface

  !> @brief Optional hook: solve (T̂ + shift) dOrb = orb exactly.
  !>
  !> @details
  !> Optional, only grids with a fast exact kinetic inverse bind it (Fem: fast
  !> diagonalization of the separable Laplacian); it stays null otherwise, so
  !> callers must test `associated(SysKinetic_SolveShifted)` and fall back
  !> (e.g. to no preconditioning). The shift is complex: real for
  !> preconditioners of eigenproblems, (T̂ + σ)⁻¹, and 2/(i dt) = −2i/dt for
  !> the Crank–Nicolson factor 1 + i dt T̂/2 = (i dt/2) (T̂ − 2i/dt).
  procedure(I_SysKinetic_SolveShifted), pointer :: SysKinetic_SolveShifted => null()
  abstract interface
    !> @brief Exact shifted solve: dOrb = (T̂ + shift)⁻¹ · orb
    !>
    !> @param[out] dOrb   Solution of (T̂ + shift) dOrb = orb (length = nGridPoints)
    !> @param[in]  orb    Right-hand side on the spatial grid
    !> @param[in]  shift  Complex shift; T̂ + shift must be invertible (for a
    !>                    real shift: shift > −(smallest eigenvalue of T̂))
    !> @param[in]  bt_    Optional body-type index selecting the mass (default = 1)
    subroutine I_SysKinetic_SolveShifted(dOrb, orb, shift, bt_)
      import :: I32, R64
      !> Output: solution of (T̂ + shift) dOrb = orb
      complex(R64), intent(out), contiguous :: dOrb(:)
      !> Input: right-hand side
      complex(R64), intent(in), contiguous  :: orb(:)
      !> Complex shift (real for preconditioners, −2i/dt for Crank–Nicolson)
      complex(R64), intent(in)              :: shift
      !> Body type index (1-based) selecting the mass from the bodyMass array
      integer(I32), intent(in), optional    :: bt_
    end subroutine
  end interface

end module
