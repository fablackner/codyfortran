! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file M_DiagonalizerList_Lobpcg.f90
!> @brief Preconditioned block eigensolver backend (LOBPCG).
!>
!> @details
!> Locally Optimal Block Preconditioned Conjugate Gradient (Knyazev 2001) for
!> the lowest `nEvals` eigenpairs of an operator A that is self-adjoint in the
!> metric ⟨·,·⟩ given by the optional diagonal `metricWeights` (Euclidean if
!> unallocated; Gram matrices with BLAS-3), e.g. A = W⁻¹H on a grid with
!> diagonal quadrature weights W.
!> Unlike ARPACK it accepts a preconditioner (`ApplyPrecOnVec`, ≈ (A − θ)⁻¹)
!> and converges in a number of iterations that does not grow with the
!> spectral radius of A when the preconditioner captures the stiff part
!> (e.g. (T + σ)⁻¹ for a Fock operator on a finite-element grid). It is
!> warm-started from the eigenvectors of the previous call, which makes
!> repeated solves of slowly changing operators (SCF iterations) cheap.
!>
!> **Algorithm** (robust variant, Hetmaniuk & Lehoucq 2006; Duersch et al.
!> 2018): block X of m = nEvals + nGuard columns, metric-orthonormal, with
!> A·X stored. Per iteration
!> 1. residuals R = AX − XΘ, norms in the metric; column j is converged when
!>    ‖r_j‖ ≤ tol·max(1, |θ_j|); stop when all j ≤ nEvals are converged;
!>    converged columns are soft-locked (no new directions, still in the
!>    Rayleigh–Ritz step)
!> 2. W = ApplyPrecOnVec(r_j, θ_j) for the active columns (W = R without
!>    preconditioner)
!> 3. W metric-orthonormalized against [X, P] and within (block classical
!>    Gram–Schmidt twice, Cholesky-QR with a fallback to an eigen-
!>    decomposition of the Gram matrix that drops near-dependent directions);
!>    a preconditioned residual without a (numerically) nonzero component
!>    outside [X, P] (e.g. a nearly exact shift-and-invert preconditioner,
!>    (A − θ_j)⁻¹ r_j ≈ x_j) is replaced by the raw residual r_j, which is
!>    orthogonal to X; only if that vanishes too is the direction dropped
!> 4. A·W (the only matrix–vector products of the iteration)
!> 5. Rayleigh–Ritz on the metric-orthonormal basis [X, P, W] (ZHEEV of the
!>    projected Hermitian matrix; the [X, P] block is known from the last
!>    step and needs no inner products)
!> 6. X ← lowest m Ritz vectors; P ← metric-orthonormal basis of the
!>    directions in which the active columns moved, orthogonal to the new X
!>    (computed in coefficient space); A·X, A·P updated linearly
!>
!> **Early exits** (both return `convergedQ = .false.` and `stalledQ =
!> .true.` with a one-line warning, instead of running to maxIter):
!> - stagnation: over the last STALL_WINDOW = 10 iterations the largest
!>   relative residual of the wanted columns moved by at most
!>   1 − STALL_FACTOR = 10 % of its smallest value there, summing all its
!>   changes (decreases and increases). Such a residual decreased by less
!>   than 10 % over the window and does not fluctuate: a floor. Floors come
!>   from operators that are self-adjoint in the metric only approximately,
!>   e.g. the exact exchange on grid.fem, Hermitian only up to the multipole
!>   truncation of the Poisson boundary values: the Ritz pairs reach a fixed
!>   point, the residual stays constant, and further iterations only cost
!>   matrix–vector products (H₂ RHF with exact exchange: floors of 7e-8 …
!>   3e-7 relative, detected after 19–49 iterations instead of maxIter).
!>   Summing the fluctuations keeps erratic but converging histories going:
!>   without a preconditioner a stiff operator needs hundreds of iterations
!>   with plateaus over which the residual drops by less than 10 % while it
!>   fluctuates by 40 % or more (T_Ne3d_17, identity preconditioner on FEDVR).
!>   A floor that fluctuates by more than 10 % (e.g. round-off for a tol
!>   below the attainable accuracy) is not detected and ends at maxIter;
!> - no new search direction is left (every new direction vanished).
!>
!> **JSON Configuration:**
!> @code{.json}
!> "diagonalizerList": {
!>   "lobpcg": { "nEvals": 5, "nGuard": 2, "maxIter": 200, "tol": 1e-8,
!>               "warmStartQ": true, "printLevel": 0 }
!> }
!> @endcode
!>
!> | Parameter    | Type | Default | Description                                    |
!> |--------------|------|---------|------------------------------------------------|
!> | `nEvals`     | int  | 1       | Number of wanted (lowest) eigenpairs           |
!> | `nGuard`     | int  | 2       | Extra block columns that speed up convergence  |
!> | `maxIter`    | int  | 200     | Iteration limit (warning, no stop, if reached) |
!> | `tol`        | real | 1e-8    | ‖r_j‖ ≤ tol·max(1, abs(θ_j)), metric norms      |
!> | `warmStartQ` | bool | true    | Start from the previous evecs if allocated     |
!> | `printLevel` | int  | 0       | 1: iterations and residuals                    |
!>
!> **Outputs:** `evals(1:nEvals)` ascending and `evecs(:, 1:nEvals)`
!> orthonormal in the metric: the nEvals best approximations, also after an
!> unconverged call (maxIter or an early exit). `nFound` is the number of
!> converged wanted columns (‖r_j‖ ≤ tol·max(1, |θ_j|)), nEvals on success;
!> `convergedQ` (component of T_DiagonalizerList_E) is `nFound .eq. nEvals`.
!> After an unconverged call the converged columns need not be the leading
!> ones; `resNorms` tells them apart. `evecs` always holds the whole block
!> (nEvals + nGuard columns, also for `evecsQ = .false.`), so that the next
!> call reuses the guard vectors. `SetStartVectors` seeds the start block of
!> the next call instead (`dim` rows, any number of columns, possibly
!> linearly dependent; e.g. before the first call); missing or dependent
!> columns are filled with deterministic pseudo-random vectors.
!>
!> @see M_DiagonalizerList for the callbacks, A_DiagonalizerList.md
module M_DiagonalizerList_Lobpcg
  use M_Utils_Types
  use M_DiagonalizerList, only: T_DiagonalizerList_E

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> @brief Allocate a LOBPCG diagonalizer instance.
    !> @param[out] e     Polymorphic pointer to allocated instance
    !> @param[in]  path  JSON configuration path
    module subroutine DiagonalizerList_Lobpcg_Allocate(e, path)
      class(T_DiagonalizerList_E), allocatable, intent(out) :: e
      character(len=*), intent(in) :: path
    end subroutine
  end interface

  !=============================================================================
  ! module types
  !=============================================================================

  !-----------------------------------------------------------------------------
  !> @brief LOBPCG eigensolver implementation.
  !>
  !> @details
  !> Extends `T_DiagonalizerList_E`; uses `ApplyMatOnVec`, the optional
  !> `ApplyPrecOnVec` callback and the optional `metricWeights`. Diagnostics
  !> of the last `Diagonalize` call are kept in `nIter`, `nMatVec`,
  !> `stalledQ`, `nRawResiduals` and `resNorms`, besides the generic
  !> `convergedQ` and `nFound`.
  !-----------------------------------------------------------------------------
  type, extends(T_DiagonalizerList_E) :: T_DiagonalizerList_E_Lobpcg
    !> Extra block columns beyond nEvals (not required to converge).
    integer(I32) :: nGuard = 2
    !> Maximal number of iterations per call.
    integer(I32) :: maxIter = 200
    !> Relative residual tolerance: ‖r_j‖ ≤ tol·max(1, |θ_j|).
    real(R64) :: tol = 1e-8_R64
    !> Start from the eigenvectors of the previous call if available.
    logical :: warmStartQ = .true.
    !> Iterations of the last call (0 if the start block was converged).
    integer(I32) :: nIter = 0
    !> Matrix–vector products of the last call.
    integer(I32) :: nMatVec = 0
    !> Whether the last call stopped early for lack of progress (stagnation
    !> exit or no new search direction left; convergedQ is then .false.).
    logical :: stalledQ = .false.
    !> Preconditioned residuals of the last call that had no component
    !> outside [X, P] and were replaced by the raw residuals.
    integer(I32) :: nRawResiduals = 0
    !> Metric residual norms ‖A x_j − θ_j x_j‖ of the block after the last call.
    real(R64), allocatable :: resNorms(:)
  contains
    procedure :: Setup      !< Prepare for diagonalization
    procedure :: Fabricate  !< Read JSON parameters
    procedure :: Diagonalize !< Run the LOBPCG iteration
    procedure :: SetStartVectors !< Seed the start block of the next call
  end type

  interface
    !> @brief Read LOBPCG parameters from JSON configuration.
    module subroutine Fabricate(this)
      class(T_DiagonalizerList_E_Lobpcg), intent(inout) :: this
    end subroutine
  end interface

  interface
    !> @brief Setup for the LOBPCG backend (workspace is allocated per call).
    module subroutine Setup(this)
      class(T_DiagonalizerList_E_Lobpcg), intent(inout) :: this
    end subroutine
  end interface

  interface
    !> @brief Compute the lowest `nEvals` eigenpairs with LOBPCG.
    !>
    !> @param[inout] this    Diagonalizer instance (results in %evals, %evecs)
    !> @param[in]    time    Time for time-dependent operators
    !> @param[in]    evecsQ  Unused: the eigenvectors are always returned
    module subroutine Diagonalize(this, time, evecsQ)
      class(T_DiagonalizerList_E_Lobpcg), intent(inout) :: this
      real(R64), intent(in) :: time
      logical, intent(in) :: evecsQ
    end subroutine
  end interface

  interface
    !> @brief Seed the start block of the next `Diagonalize` (overrides the
    !>        no-op default of T_DiagonalizerList_E).
    !>
    !> @param[inout] this  Diagonalizer instance (vecs stored in %evecs,
    !>                     %nFound = 0: no converged pairs yet)
    !> @param[in]    vecs  Start vectors (dim rows, any number of columns)
    module subroutine SetStartVectors(this, vecs)
      class(T_DiagonalizerList_E_Lobpcg), intent(inout) :: this
      complex(R64), intent(in), contiguous :: vecs(:, :)
    end subroutine
  end interface

end module
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
