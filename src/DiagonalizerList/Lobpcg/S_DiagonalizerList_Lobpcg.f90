! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file S_DiagonalizerList_Lobpcg.f90
!> @brief Implementation of the LOBPCG eigensolver backend.
!>
!> @details
!> **Storage** (n = dim, m = nEvals + nGuard block columns):
!> - `xp(n, 2m)` = [X | P]: X = xp(:, 1:m) metric-orthonormal Ritz vectors,
!>   P = xp(:, m+1:m+kP) metric-orthonormal and orthogonal to X
!> - `w(n, m)`: residuals → preconditioned residuals → new directions W
!> - `axp`, `aw`: A applied to the columns of `xp` and `w`
!> - `hKnown`: projected matrix of [X | P], diag(θ) ⊕ Uᴴ diag(θ_rest) U
!>
!> Every column of W gets exactly one matrix–vector product; A·X and A·P are
!> updated linearly with the Ritz coefficients. Since the basis [X, P, W] is
!> metric-orthonormal, the Rayleigh–Ritz step is a standard Hermitian
!> eigenproblem of size m + kP + kW ≤ 3m, and only the rows of W of the
!> projected matrix need inner products. P is chosen as in Hetmaniuk &
!> Lehoucq (2006): with the unitary Ritz coefficient matrix Y = [Y1 Y2]
!> (Y1 the m selected columns) the active columns J moved, within
!> span[X_old, X_new] and orthogonal to X_new, along Y2·(Y(1:m, m+1:)ᴴ Y(1:m, J));
!> an SVD of that small matrix gives a metric-orthonormal P without
!> inner products.
!>
!> **Inner products** use one set of diagonal metric weights `wts`, set up
!> by Diagonalize once per call: `metricWeights`, or unit weights for the
!> Euclidean metric (as SIL of IntegratorList does). Gram blocks are ZGEMMs
!> on weight-scaled row chunks (`WeightedGram`), norms weighted sums; the
!> block kernels take `wts`, not the solver object. In-place block updates
!> V ← V·C run over row chunks.
!>
!> **Structure.** `Diagonalize` holds the iteration state and the main loop
!> with every exit; its contained steps share that state by host association
!> (as the steps of SIL in IntegratorList do) and follow items 1–6 of the
!> algorithm in M_DiagonalizerList_Lobpcg: `BuildStartBlock`,
!> `UpdateResiduals` (1), `BuildSearchDirections` (2–4) and
!> `RayleighRitzStep` (5–6). The module-level procedures are stateless: the
!> stagnation rule `StagnatedQ`, the projected matrix of [X | P]
!> (`ResetKnownProjection`), the orthonormalization and the metric and block
!> kernels.
submodule(M_DiagonalizerList_Lobpcg) S_DiagonalizerList_Lobpcg

  implicit none

  !> Rows per chunk of the in-place block updates V ← V·C.
  integer(I32), parameter :: CHUNK_ROWS = 1024
  !> Minimal work (vector length × number of vector products) for threading.
  integer(I64), parameter :: MIN_PARALLEL_WORK = 65536_I64
  !> Smallest admissible squared pivot of the Cholesky factor of the unit-
  !> diagonal Gram matrix; below it the eigen-decomposition based
  !> orthonormalization (with dropping) is used instead of Cholesky-QR.
  real(R64), parameter :: CHOL_PIVOT_TOL = 1e-10_R64
  !> Directions whose eigenvalue of the unit-diagonal Gram matrix is below
  !> DROP_TOL times the largest one are dropped as linearly dependent.
  real(R64), parameter :: DROP_TOL = 1e-12_R64
  !> Columns whose norm the projection against [X, P] reduces by more than
  !> this factor carry no new direction: a preconditioned residual is then
  !> replaced by the raw residual, which is dropped only if it vanishes too.
  real(R64), parameter :: VANISH_TOL = 1e-10_R64
  !> Stagnation exit: the iteration stops (convergedQ = .false., stalledQ =
  !> .true., one-line warning) when, over the last STALL_WINDOW iterations,
  !> the largest relative residual of the wanted columns moved by at most
  !> 1 − STALL_FACTOR = 10 % of its smallest value in that window, counting
  !> every change (sum of the absolute changes between consecutive
  !> iterations). Such a history decreased by less than 10 % over the
  !> window and is stationary: the floor of an operator that is self-adjoint
  !> in the metric only approximately (the Ritz pairs reach a fixed point
  !> and the residual stays constant), e.g. the exact exchange on grid.fem.
  !> Counting the fluctuations keeps erratic but converging histories going:
  !> without a preconditioner a stiff operator converges in hundreds of
  !> iterations with plateaus of 10–30 iterations over which the residual
  !> drops by less than 10 % but fluctuates by 40 % or more (T_Ne3d_17).
  integer(I32), parameter :: STALL_WINDOW = 10
  real(R64), parameter :: STALL_FACTOR = 0.9_R64
  !> Deviation from orthonormality after two passes that triggers a third.
  real(R64), parameter :: REORTHO_TOL = 1e-8_R64
  !> Singular values of the coefficient matrix of P below this are dropped.
  real(R64), parameter :: P_DROP_TOL = 1e-14_R64
  !> Attempts to complete a linearly independent start block.
  integer(I32), parameter :: MAX_START_ATTEMPTS = 5
  !> Fixed seed and constants of the Park–Miller (MINSTD) generator used for
  !> the deterministic pseudo-random start vectors.
  integer(I64), parameter :: RANDOM_SEED = 20250707_I64
  integer(I64), parameter :: MINSTD_A = 48271_I64
  integer(I64), parameter :: MINSTD_M = 2147483647_I64

  complex(R64), parameter :: ONE = (1.0_R64, 0.0_R64)
  complex(R64), parameter :: ZERO = (0.0_R64, 0.0_R64)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Allocate a new LOBPCG diagonalizer instance.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine DiagonalizerList_Lobpcg_Allocate(e, path)
    use M_Utils_UnusedVariables
    class(T_DiagonalizerList_E), allocatable, intent(out) :: e
    character(len=*), intent(in) :: path

    if (.false.) call UnusedVariables_Mark(path)

    allocate (T_DiagonalizerList_E_Lobpcg :: e)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Read configuration parameters for the LOBPCG backend.
!>
!> @details
!> Reads from JSON at `this%path`: `nEvals` (1), `nGuard` (2), `maxIter`
!> (200), `tol` (1e-8), `warmStartQ` (true), `printLevel` (0).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Fabricate(this)
    use M_Utils_Json
    use M_Utils_Say

    class(T_DiagonalizerList_E_Lobpcg), intent(inout) :: this

    call Say_Fabricate(this % path)

    !------------------------------------
    ! Read configuration parameters
    !------------------------------------

    this % nEvals = Json_Get("nEvals", 1, path_=this % path)
    this % nGuard = Json_Get("nGuard", 2, path_=this % path)
    this % maxIter = Json_Get("maxIter", 200, path_=this % path)
    this % tol = Json_Get("tol", 1e-8_R64, path_=this % path)
    this % warmStartQ = Json_Get("warmStartQ", .true., path_=this % path)
    this % printLevel = Json_Get("printLevel", 0, path_=this % path)

    !------------------------------------
    ! Validate
    !------------------------------------

    if (this % nEvals < 1) error stop "diagonalizerList lobpcg: nEvals must be at least 1"
    if (this % nEvals > this % dim) error stop "diagonalizerList lobpcg: nEvals exceeds the dimension"
    if (this % nGuard < 0) error stop "diagonalizerList lobpcg: nGuard must not be negative"
    if (this % maxIter < 0) error stop "diagonalizerList lobpcg: maxIter must not be negative"
    if (.not. (this % tol > 0.0_R64)) error stop "diagonalizerList lobpcg: tol must be positive"

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Setup for the LOBPCG backend (workspace is allocated per call).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Setup(this)
    use M_Utils_Say

    class(T_DiagonalizerList_E_Lobpcg), intent(inout) :: this

    call Say_Setup(this % path)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Seed the start block of the next `Diagonalize`.
!>
!> @details
!> The vectors take the place of the previous eigenvectors (`evecs`, with
!> `nFound = 0`: none of them is a converged pair), so the next call starts
!> from them as from a warm start (with `warmStartQ`); see Diagonalize.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SetStartVectors(this, vecs)
    class(T_DiagonalizerList_E_Lobpcg), intent(inout) :: this
    complex(R64), intent(in), contiguous :: vecs(:, :)

    if (size(vecs, 1) .ne. this % dim) then
      error stop "diagonalizerList lobpcg: start vectors must have the size dim of the operator"
    end if

    this % evecs = vecs
    this % nFound = 0

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Run the LOBPCG iteration for the lowest `nEvals` eigenpairs.
!>
!> @details
!> Start block (`BuildStartBlock`): the previous `evecs` (if `warmStartQ` and
!> their row count matches; or the seed of `SetStartVectors`), completed by
!> deterministic pseudo-random columns, metric-orthonormalized; then one
!> Rayleigh–Ritz step on it. Iterations as described in
!> M_DiagonalizerList_Lobpcg, one contained step per phase: `UpdateResiduals`
!> (item 1), `BuildSearchDirections` (items 2–4) and `RayleighRitzStep`
!> (items 5–6). The main loop stops when all wanted columns satisfy
!> ‖r_j‖ ≤ tol·max(1, |θ_j|), when `maxIter` is reached (one-line warning),
!> when the largest wanted residual stagnates (`StagnatedQ`; one-line
!> warning) or when no new direction is left (one-line warning). A
!> preconditioned residual without a component outside [X, P] is replaced by
!> its raw residual (staged in `aw`, which is free until A·W is applied). On
!> output `nFound` counts the converged wanted columns; `evals`/`evecs` hold
!> the nEvals (evecs: whole block) best approximations in any case.
!>
!> The steps share the iteration state declared at the top of Diagonalize by
!> host association, as the steps of `Integrate` in SIL (IntegratorList) do;
!> all their other variables are local.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Diagonalize(this, time, evecsQ)
    use M_Utils_LapackLib, only: LapackLib_DiagonalizeGeneric, LapackLib_Svd
    use M_Utils_UnusedVariables

    class(T_DiagonalizerList_E_Lobpcg), intent(inout) :: this
    real(R64), intent(in) :: time
    logical, intent(in) :: evecsQ

    ! Iteration state (shared with the contained steps):
    ! [X | P] (X = xp(:, 1:m), P = xp(:, m+1:m+kP)) and A·[X | P]; w holds
    ! the residuals, then the preconditioned residuals, then the new
    ! directions W; aw holds A·W (before: the raw residuals as substitutes)
    complex(R64), allocatable :: xp(:, :), axp(:, :), w(:, :), aw(:, :)
    ! projected matrix of [X | P], known from the last Rayleigh–Ritz step
    complex(R64), allocatable :: hKnown(:, :)
    ! preconditioner output
    complex(R64), allocatable :: tmp(:)
    ! Ritz values of X, metric residual norms, and per iteration the largest
    ! relative residual of the wanted columns (stagnation test)
    real(R64), allocatable :: theta(:), resNorm(:), maxRelResHistory(:)
    ! metric weights of all inner products
    real(R64), allocatable :: wts(:)
    ! the nAct active (not converged) block columns
    integer(I32), allocatable :: active(:)
    ! converged block columns (soft-locked)
    logical, allocatable :: convQ(:)
    ! dim, block size, wanted columns, columns of P, of [X | P] and of W,
    ! active columns, iteration
    integer(I32) :: n, m, nE, kP, q, kW, nAct, iter
    ! largest relative residual of the wanted columns
    real(R64) :: maxRelRes
    logical :: precQ

    if (.false.) call UnusedVariables_Mark(evecsQ)

    n = this % dim
    nE = this % nEvals
    m = min(nE + this % nGuard, n)
    precQ = associated(this % ApplyPrecOnVec)

    ! The metric: metricWeights, or unit weights (Euclidean) if unallocated
    if (allocated(this % metricWeights)) then
      if (size(this % metricWeights) .ne. n) then
        error stop "diagonalizerList lobpcg: metricWeights must have the size dim of the operator"
      end if
      wts = this % metricWeights
    else
      allocate (wts(n), source=1.0_R64)
    end if

    allocate (xp(n, 2 * m), axp(n, 2 * m), w(n, m), aw(n, m), tmp(n))
    allocate (hKnown(2 * m, 2 * m), theta(m), resNorm(m), convQ(m), active(m))
    allocate (maxRelResHistory(0:this % maxIter))

    this % nIter = 0
    this % nMatVec = 0
    this % nRawResiduals = 0
    this % convergedQ = .false.
    this % stalledQ = .false.

    if (this % printLevel > 0) then
      write (*, *)
      write (*, '(1X, 3A, I0, A, I0, A, I0)') "LOBPCG (", this % path, "): dim ", n, ", block ", m, &
        ", wanted ", nE
    end if

    call BuildStartBlock

    !------------------------------------
    ! Iterate
    !------------------------------------

    iter = 0
    maxRelRes = huge(1.0_R64)

    do
      q = m + kP
      call UpdateResiduals

      if (all(convQ(1:nE))) then
        this % convergedQ = .true.
        exit
      end if
      if (iter >= this % maxIter) then
        write (*, '(1X, 3A, I0, A, ES10.3)') "Warning: LOBPCG (", this % path, ") not converged after ", iter, &
          " iterations, max relative residual ", maxRelRes
        exit
      end if
      if (StagnatedQ(maxRelResHistory(0:iter))) then
        write (*, '(1X, 3A, I0, A, I0, A, ES10.3)') "Warning: LOBPCG (", this % path, ") stagnates: no progress over ", &
          STALL_WINDOW, " iterations, stopped after ", iter, " at max relative residual ", maxRelRes
        this % stalledQ = .true.
        exit
      end if

      iter = iter + 1

      call BuildSearchDirections
      if (kW .eq. 0) then
        write (*, '(1X, 3A, I0, A, ES10.3)') "Warning: LOBPCG (", this % path, &
          "): no new search direction left, stopped after ", iter - 1, " iterations at max relative residual ", maxRelRes
        iter = iter - 1
        this % stalledQ = .true.
        exit
      end if

      call RayleighRitzStep
    end do

    if (this % printLevel > 0 .and. this % convergedQ) then
      write (*, '(1X, A, I0, A, I0, A)') "LOBPCG converged after ", iter, " iterations (", this % nMatVec, &
        " matrix-vector products)"
      write (*, *)
    end if
    if (this % printLevel > 0 .and. this % nRawResiduals > 0) then
      write (*, '(1X, A, I0, A)') "LOBPCG: ", this % nRawResiduals, &
        " preconditioned residuals without a component outside [X, P] replaced by raw residuals"
    end if

    !------------------------------------
    ! Outputs: all nEvals best
    ! approximations (the whole block is
    ! kept in evecs for the next warm
    ! start); nFound counts the converged
    ! wanted columns
    !------------------------------------

    this % evals = theta(1:nE)
    this % evecs = xp(:, 1:m)
    this % nFound = int(count(convQ(1:nE)), I32)
    this % nIter = iter
    this % resNorms = resNorm

  contains

    !> Start block X: the previous eigenvectors (or the SetStartVectors seed)
    !> and pseudo-random columns, metric-orthonormalized (columns that are
    !> linearly dependent are replaced by new pseudo-random ones); A·X; one
    !> Rayleigh–Ritz step on X. Leaves hKnown = diag(θ) and no P (kP = 0).
    subroutine BuildStartBlock
      complex(R64), allocatable :: h(:, :), y(:, :)
      real(R64), allocatable :: thetaAll(:)
      integer(I32) :: nStart, kX, kNew, attempt, j
      integer(I64) :: seed

      nStart = 0
      if (this % warmStartQ .and. allocated(this % evecs)) then
        if (size(this % evecs, 1) .eq. n) then
          nStart = min(int(size(this % evecs, 2), I32), m)
          w(:, 1:nStart) = this % evecs(:, 1:nStart)
        end if
      end if
      seed = RANDOM_SEED
      call FillRandom(w(:, nStart + 1:m), seed)

      kX = 0
      do attempt = 1, MAX_START_ATTEMPTS
        call OrthonormalizeBlock(kNew, w(:, 1:m - kX), wts, xp(:, 1:kX))
        xp(:, kX + 1:kX + kNew) = w(:, 1:kNew)
        kX = kX + kNew
        if (kX .eq. m) exit
        call FillRandom(w(:, 1:m - kX), seed)
      end do
      if (kX < m) error stop "diagonalizerList lobpcg: no linearly independent start block found"

      do j = 1, m
        call this % ApplyMatOnVec(axp(:, j), xp(:, j), time)
      end do
      this % nMatVec = m

      ! Rayleigh–Ritz on the start block
      allocate (h(m, m))
      call CalcProjected(h, xp(:, 1:m), axp(:, 1:m), wts)
      call LapackLib_DiagonalizeGeneric(thetaAll, y, h, .true.)
      call RightMultiplyInPlace(xp, y, n, m, m)
      call RightMultiplyInPlace(axp, y, n, m, m)
      theta(:) = thetaAll(1:m)

      call ResetKnownProjection(hKnown, theta)
      kP = 0

    end subroutine

    !> Item 1: residuals R = AX − XΘ of all block columns (staged in w), their
    !> metric norms and convergence (converged columns are soft-locked); the
    !> largest relative residual of the wanted columns enters the history of
    !> the stagnation test.
    subroutine UpdateResiduals

      call CalcResiduals(w, xp, axp, theta, n, m)
      call CalcNorms(resNorm, w(:, 1:m), wts)
      convQ(:) = resNorm(:) <= this % tol * max(1.0_R64, abs(theta(:)))
      maxRelRes = maxval(resNorm(1:nE) / max(1.0_R64, abs(theta(1:nE))))
      maxRelResHistory(iter) = maxRelRes

      if (this % printLevel > 0) then
        write (*, '(1X, A, I5, A, I4, A, ES10.3)') "LOBPCG iteration", iter, ": active ", count(.not. convQ), &
          ", max relative residual (wanted) ", maxRelRes
        write (*, '(1X, A, *(ES10.3, 1X))') "  relative residuals: ", resNorm(:) / max(1.0_R64, abs(theta(:)))
      end if

    end subroutine

    !> Items 2–4: the preconditioned residuals of the active columns,
    !> compacted into w, metric-orthonormalized against [X, P] and among
    !> themselves into the kW new directions W, and A·W (skipped if kW = 0).
    !> The raw residuals of the active columns are staged in aw (free until
    !> A·W is applied) as substitutes for preconditioned ones without a
    !> component outside [X, P].
    subroutine BuildSearchDirections
      integer(I32) :: nRaw, j

      nAct = 0
      do j = 1, m
        if (convQ(j)) cycle
        nAct = nAct + 1
        active(nAct) = j
        if (precQ) then
          call this % ApplyPrecOnVec(tmp, w(:, j), theta(j), time)
          aw(:, nAct) = w(:, j)
          w(:, nAct) = tmp(:)
        else if (nAct .ne. j) then
          w(:, nAct) = w(:, j)
        end if
      end do

      ! New directions: metric-orthonormal and orthogonal to [X, P]
      if (precQ) then
        call OrthonormalizeBlock(kW, w(:, 1:nAct), wts, xp(:, 1:q), nSubstituted_=nRaw, substitutes_=aw(:, 1:nAct))
        this % nRawResiduals = this % nRawResiduals + nRaw
      else
        call OrthonormalizeBlock(kW, w(:, 1:nAct), wts, xp(:, 1:q))
      end if
      if (kW .eq. 0) return

      do j = 1, kW
        call this % ApplyMatOnVec(aw(:, j), w(:, j), time)
      end do
      this % nMatVec = this % nMatVec + kW

    end subroutine

    !> Items 5–6: Rayleigh–Ritz on the metric-orthonormal basis [X, P, W]
    !> (the [X, P] block of the projected matrix is hKnown; only the rows of
    !> W need inner products); X ← the lowest m Ritz vectors, P ← the
    !> metric-orthonormal directions in which the active columns moved,
    !> orthogonal to the new X (computed in coefficient space); A·X and A·P
    !> updated linearly; hKnown of the new [X, P].
    subroutine RayleighRitzStep
      complex(R64), allocatable :: h(:, :), hc(:, :), hw(:, :), y(:, :), mP(:, :), uP(:, :), vtP(:, :), coef(:, :)
      real(R64), allocatable :: thetaAll(:), svP(:)
      integer(I32) :: s, nKeepP

      ! Projected matrix on the metric-orthonormal basis [X, P, W]; the
      ! [X, P] block is known from the previous Rayleigh–Ritz step
      s = q + kW
      allocate (h(s, s), hc(q, kW), hw(kW, kW))
      h(1:q, 1:q) = hKnown(1:q, 1:q)
      call CalcGram(hc, xp(:, 1:q), aw(:, 1:kW), wts)
      call CalcProjected(hw, w(:, 1:kW), aw(:, 1:kW), wts)
      h(1:q, q + 1:s) = hc(:, :)
      h(q + 1:s, 1:q) = conjg(transpose(hc(:, :)))
      h(q + 1:s, q + 1:s) = hw(:, :)
      deallocate (hc, hw)

      call LapackLib_DiagonalizeGeneric(thetaAll, y, h, .true.)

      ! P: directions in which the active columns moved, orthogonal to the
      ! new X, as a metric-orthonormal combination of the other Ritz vectors
      allocate (mP(s - m, nAct))
      mP(:, :) = matmul(conjg(transpose(y(1:m, m + 1:s))), y(1:m, active(1:nAct)))
      call LapackLib_Svd(uP, svP, vtP, mP)
      nKeepP = count(svP(:) > P_DROP_TOL)

      ! [X, P] ← [X, P, W]·[Y1, Y2·U]  (A·[X, P] alike)
      allocate (coef(s, m + nKeepP))
      coef(:, 1:m) = y(:, 1:m)
      if (nKeepP > 0) coef(:, m + 1:m + nKeepP) = matmul(y(:, m + 1:s), uP(:, 1:nKeepP))
      call CombineInPlace(xp, w, coef(1:q, :), coef(q + 1:s, :), n, q, kW, m + nKeepP)
      call CombineInPlace(axp, aw, coef(1:q, :), coef(q + 1:s, :), n, q, kW, m + nKeepP)

      ! Projected matrix of the new [X, P]: diag(θ) ⊕ Uᴴ diag(θ_rest) U
      theta(:) = thetaAll(1:m)
      call ResetKnownProjection(hKnown, theta, uP_=uP(:, 1:nKeepP), thetaRest_=thetaAll(m + 1:s))
      kP = nKeepP

    end subroutine

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Stagnation test on the history of the largest wanted residual.
!>
!> @details
!> True when over the last STALL_WINDOW iterations the largest relative
!> residual of the wanted columns moved by at most 1 − STALL_FACTOR of its
!> smallest value in that window, counting every change (the sum of the
!> absolute changes between consecutive iterations): it decreased by less
!> than 10 % without fluctuating, a floor. False before STALL_WINDOW
!> iterations.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure function StagnatedQ(history) result(res)
    logical :: res
    !> Largest relative residual of the wanted columns of iterations 0 … iter.
    real(R64), intent(in), contiguous :: history(0:)

    integer(I32) :: iter
    real(R64) :: variation

    iter = ubound(history, 1, kind=I32)
    res = .false.
    if (iter < STALL_WINDOW) return

    variation = sum(abs(history(iter - STALL_WINDOW + 1:iter) &
                        - history(iter - STALL_WINDOW:iter - 1)))
    res = variation <= (1.0_R64 - STALL_FACTOR) * minval(history(iter - STALL_WINDOW:iter))

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Projected matrix of [X | P] after a Rayleigh–Ritz step.
!>
!> @details
!> hKnown = diag(θ) ⊕ Uᴴ diag(θ_rest) U, zero elsewhere: X holds the Ritz
!> vectors with the Ritz values θ, and P = (other Ritz vectors)·U, whose
!> block is Hermitized. Without `uP_` and `thetaRest_` (the start block has
!> no P): hKnown = diag(θ).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ResetKnownProjection(hKnown, theta, uP_, thetaRest_)
    !> Projected matrix of [X | P] (at least m + size(uP_, 2) rows/columns).
    complex(R64), intent(out), contiguous :: hKnown(:, :)
    !> Ritz values θ of the m columns of X.
    real(R64), intent(in), contiguous :: theta(:)
    !> Coefficients U of the P columns in the other Ritz vectors.
    complex(R64), intent(in), contiguous, optional :: uP_(:, :)
    !> Ritz values θ_rest of the other Ritz vectors (one per row of uP_).
    real(R64), intent(in), contiguous, optional :: thetaRest_(:)

    complex(R64), allocatable :: uTheta(:, :)
    integer(I32) :: m, nRest, nKeepP, i, j

    m = int(size(theta), I32)

    hKnown(:, :) = ZERO
    do j = 1, m
      hKnown(j, j) = cmplx(theta(j), 0.0_R64, kind=R64)
    end do

    if (.not. (present(uP_) .and. present(thetaRest_))) return
    nRest = int(size(uP_, 1), I32)
    nKeepP = int(size(uP_, 2), I32)
    if (nKeepP .eq. 0) return

    allocate (uTheta(nRest, nKeepP))
    do j = 1, nKeepP
      do i = 1, nRest
        uTheta(i, j) = thetaRest_(i) * uP_(i, j)
      end do
    end do
    hKnown(m + 1:m + nKeepP, m + 1:m + nKeepP) = matmul(conjg(transpose(uP_)), uTheta)
    hKnown(m + 1:m + nKeepP, m + 1:m + nKeepP) = 0.5_R64 * (hKnown(m + 1:m + nKeepP, m + 1:m + nKeepP) &
                                                   + conjg(transpose(hKnown(m + 1:m + nKeepP, m + 1:m + nKeepP))))

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Metric-orthonormalize a block against an orthonormal set.
!>
!> @details
!> Orthonormalizes the columns of `v` in the metric against the metric-
!> orthonormal columns of `q` and among themselves: block classical Gram–
!> Schmidt twice ("twice is enough"; a third pass if the second one still
!> sees a deviation above REORTHO_TOL), each pass followed by
!> `OrthonormalizeWithGram`. A column annihilated by the first projection
!> (relative norm below VANISH_TOL) carries no new direction: it is replaced
!> by the same column of `substitutes_` (if present; projected against q
!> and tested alike) and dropped only if that vanishes too. Near-dependent
!> directions are dropped; the `kOut` kept columns are compacted into
!> v(:, 1:kOut).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine OrthonormalizeBlock(kOut, v, wts, q, nSubstituted_, substitutes_)
    use M_Utils_BlasLib, only: ZGEMM

    !> Number of kept (orthonormal) columns, stored in v(:, 1:kOut).
    integer(I32), intent(out) :: kOut
    !> Columns to orthonormalize (in), kept orthonormal columns (out).
    complex(R64), intent(inout), contiguous :: v(:, :)
    !> Metric weights.
    real(R64), intent(in), contiguous :: wts(:)
    !> Metric-orthonormal set to orthogonalize against (may have 0 columns).
    complex(R64), intent(in), contiguous :: q(:, :)
    !> Optional output: number of columns replaced by their substitute.
    integer(I32), intent(out), optional :: nSubstituted_
    !> Optional substitutes, one per column of v (LOBPCG: the raw residuals
    !> of preconditioned ones), for columns without a component outside q.
    complex(R64), intent(in), contiguous, optional :: substitutes_(:, :)

    complex(R64), allocatable :: c(:, :), g(:, :), sub(:, :), cSub(:, :)
    real(R64), allocatable :: normSub(:)
    integer(I32), allocatable :: keep(:), iSubOf(:)
    logical, allocatable :: vanishedQ(:)
    integer(I32) :: nRows, nQ, pass, j, kNew, nSub, iSub
    real(R64) :: devC, devG, normAfterSq, normBeforeSq

    nRows = int(size(v, 1), I32)
    nQ = int(size(q, 2), I32)
    kOut = int(size(v, 2), I32)
    if (present(nSubstituted_)) nSubstituted_ = 0
    if (kOut .eq. 0) return

    do pass = 1, 3
      devC = 0.0_R64
      devG = 0.0_R64

      ! Projection against q: v ← v − q (qᴴ B v)
      if (nQ > 0) then
        allocate (c(nQ, kOut))
        call CalcGram(c, q, v(:, 1:kOut), wts)
        call ZGEMM('N', 'N', nRows, kOut, nQ, -ONE, q, nRows, c, nQ, ONE, v, nRows)
        devC = maxval(abs(c))
      end if

      ! Gram matrix of the projected columns
      allocate (g(kOut, kOut))
      call CalcProjected(g, v(:, 1:kOut), v(:, 1:kOut), wts)

      if (pass .eq. 1) then
        ! Columns without a new direction: the projection removed all but a
        ! fraction VANISH_TOL of their norm (‖v_before‖² = ‖v_after‖² +
        ! ‖qᴴ B v‖² for metric-orthonormal q), or they vanish entirely
        allocate (vanishedQ(kOut))
        do j = 1, kOut
          normAfterSq = real(g(j, j), kind=R64)
          normBeforeSq = normAfterSq
          if (nQ > 0) normBeforeSq = normBeforeSq + sum(abs(c(:, j))**2)
          vanishedQ(j) = .not. (normAfterSq > 0.0_R64 .and. normAfterSq > VANISH_TOL**2 * normBeforeSq)
        end do

        ! Replace them by their substitutes, projected and tested alike
        if (present(substitutes_) .and. any(vanishedQ)) then
          nSub = int(count(vanishedQ), I32)
          allocate (sub(nRows, nSub), normSub(nSub), iSubOf(nSub))
          iSub = 0
          do j = 1, kOut
            if (.not. vanishedQ(j)) cycle
            iSub = iSub + 1
            iSubOf(iSub) = j
            sub(:, iSub) = substitutes_(:, j)
          end do
          if (nQ > 0) then
            allocate (cSub(nQ, nSub))
            call CalcGram(cSub, q, sub, wts)
            call ZGEMM('N', 'N', nRows, nSub, nQ, -ONE, q, nRows, cSub, nQ, ONE, sub, nRows)
          end if
          call CalcNorms(normSub, sub, wts)
          do iSub = 1, nSub
            normAfterSq = normSub(iSub)**2
            normBeforeSq = normAfterSq
            if (nQ > 0) normBeforeSq = normBeforeSq + sum(abs(cSub(:, iSub))**2)
            if (normAfterSq > 0.0_R64 .and. normAfterSq > VANISH_TOL**2 * normBeforeSq) then
              j = iSubOf(iSub)
              v(:, j) = sub(:, iSub)
              vanishedQ(j) = .false.
              if (present(nSubstituted_)) nSubstituted_ = nSubstituted_ + 1
            end if
          end do
          deallocate (sub, normSub, iSubOf)
          if (allocated(cSub)) deallocate (cSub)
          ! Gram matrix including the substituted columns
          call CalcProjected(g, v(:, 1:kOut), v(:, 1:kOut), wts)
        end if

        ! Drop the remaining ones
        allocate (keep(kOut))
        kNew = 0
        do j = 1, kOut
          if (vanishedQ(j)) cycle
          kNew = kNew + 1
          keep(kNew) = j
          if (kNew .ne. j) v(:, kNew) = v(:, j)
        end do
        if (kNew < kOut) then
          g = g(keep(1:kNew), keep(1:kNew))
          kOut = kNew
        end if
        deallocate (keep, vanishedQ)
      else
        do j = 1, kOut
          g(j, j) = g(j, j) - ONE
        end do
        devG = maxval(abs(g))
        do j = 1, kOut
          g(j, j) = g(j, j) + ONE
        end do
      end if
      if (allocated(c)) deallocate (c)

      if (kOut .eq. 0) return

      call OrthonormalizeWithGram(kNew, v(:, 1:kOut), g)
      deallocate (g)
      kOut = kNew
      if (kOut .eq. 0) return

      if (pass .eq. 2 .and. devC <= REORTHO_TOL .and. devG <= REORTHO_TOL) exit
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Metric-orthonormalize a block given its Gram matrix.
!>
!> @details
!> With D = diag(g)^(1/2) and the unit-diagonal Gram matrix G_s = D⁻¹ g D⁻¹:
!> Cholesky-QR v ← v D⁻¹ R⁻¹ (G_s = Rᴴ R) if the factorization succeeds with
!> all squared pivots above CHOL_PIVOT_TOL; otherwise the eigen-decomposition
!> G_s = U Λ Uᴴ gives v ← v D⁻¹ U Λ^(-1/2) restricted to λ > DROP_TOL·λ_max,
!> which drops near-dependent directions. The kOut resulting columns are
!> stored in v(:, 1:kOut).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine OrthonormalizeWithGram(kOut, v, g)
    use M_Utils_LapackLib, only: ZPOTRF, LapackLib_DiagonalizeGeneric

    !> Number of resulting orthonormal columns.
    integer(I32), intent(out) :: kOut
    !> Columns (in), orthonormal columns in v(:, 1:kOut) (out).
    complex(R64), intent(inout), contiguous :: v(:, :)
    !> Gram matrix vᴴ B v (Hermitian).
    complex(R64), intent(in), contiguous :: g(:, :)

    complex(R64), allocatable :: gs(:, :), r(:, :), t(:, :), u(:, :)
    real(R64), allocatable :: d(:), lambda(:)
    integer(I32) :: k, nRows, i, j, info, first
    logical :: choleskyQ

    k = int(size(g, 1), I32)
    nRows = int(size(v, 1), I32)

    allocate (d(k), gs(k, k), t(k, k))
    do j = 1, k
      d(j) = sqrt(max(real(g(j, j), kind=R64), tiny(1.0_R64)))
    end do
    do j = 1, k
      do i = 1, k
        gs(i, j) = g(i, j) / (d(i) * d(j))
      end do
    end do

    ! Cholesky-QR
    r = gs
    info = 0
    call ZPOTRF('U', k, r, k, info)
    choleskyQ = info .eq. 0
    if (choleskyQ) then
      do j = 1, k
        if (.not. (real(r(j, j), kind=R64)**2 > CHOL_PIVOT_TOL)) choleskyQ = .false.
      end do
    end if

    if (choleskyQ) then
      ! t = D⁻¹ R⁻¹ (R upper triangular, back substitution per column)
      t(:, :) = ZERO
      do j = 1, k
        t(j, j) = ONE / r(j, j)
        do i = j - 1, 1, -1
          t(i, j) = -sum(r(i, i + 1:j) * t(i + 1:j, j)) / r(i, i)
        end do
      end do
      do j = 1, k
        t(:, j) = t(:, j) / d(:)
      end do
      kOut = k
    else
      ! Eigen-decomposition based orthonormalization with dropping
      call LapackLib_DiagonalizeGeneric(lambda, u, gs, .true.)
      kOut = int(count(lambda(:) > DROP_TOL * lambda(k)), I32)
      first = k - kOut
      do j = 1, kOut
        t(:, j) = u(:, first + j) / (d(:) * sqrt(lambda(first + j)))
      end do
    end if

    if (kOut > 0) call RightMultiplyInPlace(v, t(:, 1:kOut), nRows, k, kOut)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Metric inner product ⟨a|b⟩ = Σ_i conj(a_i) w_i b_i.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure function MetricProduct(a, b, wts) result(res)
    complex(R64) :: res
    !> Bra vector.
    complex(R64), intent(in), contiguous :: a(:)
    !> Ket vector.
    complex(R64), intent(in), contiguous :: b(:)
    !> Metric weights.
    real(R64), intent(in), contiguous :: wts(:)

    res = sum(conjg(a(:)) * wts(:) * b(:))

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Metric norms of the columns of v.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine CalcNorms(res, v, wts)
    !> Output: res(j) = ‖v(:, j)‖ in the metric.
    real(R64), intent(out), contiguous :: res(:)
    !> Input columns.
    complex(R64), intent(in), contiguous :: v(:, :)
    !> Metric weights.
    real(R64), intent(in), contiguous :: wts(:)

    integer(I32) :: j, nCols
    logical :: threadsQ

    nCols = int(size(v, 2), I32)
    threadsQ = ParallelWorkQ(int(size(v, 1), I32), nCols)

    !$omp parallel do default(shared) private(j) schedule(dynamic) if (threadsQ)
    do j = 1, nCols
      res(j) = sqrt(max(0.0_R64, real(MetricProduct(v(:, j), v(:, j), wts), kind=R64)))
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Metric Gram block g(i, j) = ⟨a(:, i)|b(:, j)⟩.
!>
!> @details
!> ZGEMM on weight-scaled row chunks (`WeightedGram`).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine CalcGram(g, a, b, wts)
    !> Output block, size (size(a, 2), size(b, 2)).
    complex(R64), intent(out), contiguous :: g(:, :)
    !> Bra columns.
    complex(R64), intent(in), contiguous :: a(:, :)
    !> Ket columns.
    complex(R64), intent(in), contiguous :: b(:, :)
    !> Metric weights.
    real(R64), intent(in), contiguous :: wts(:)

    integer(I32) :: nRows, nA, nB

    nRows = int(size(a, 1), I32)
    nA = int(size(a, 2), I32)
    nB = int(size(b, 2), I32)
    if (nA .eq. 0 .or. nB .eq. 0) return

    call WeightedGram(g, a, b, wts, nRows, nA, nB)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Hermitian projected matrix h(i, j) = ⟨v(:, i)|av(:, j)⟩.
!>
!> @details
!> The block from `WeightedGram`, Hermitized: only its upper triangle is
!> kept, the lower one is set to the conjugate and the diagonal to its real
!> part (exact for an operator self-adjoint in the metric). With av = v this
!> is the Gram matrix of v.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine CalcProjected(h, v, av, wts)
    !> Output Hermitian matrix, size (size(v, 2), size(v, 2)).
    complex(R64), intent(out), contiguous :: h(:, :)
    !> Basis columns.
    complex(R64), intent(in), contiguous :: v(:, :)
    !> Operator applied to the basis columns.
    complex(R64), intent(in), contiguous :: av(:, :)
    !> Metric weights.
    real(R64), intent(in), contiguous :: wts(:)

    integer(I32) :: nRows, k, i, j

    nRows = int(size(v, 1), I32)
    k = int(size(v, 2), I32)
    if (k .eq. 0) return

    call WeightedGram(h, v, av, wts, nRows, k, k)

    do j = 1, k
      h(j, j) = cmplx(real(h(j, j), kind=R64), 0.0_R64, kind=R64)
      do i = j + 1, k
        h(i, j) = conjg(h(j, i))
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Weighted Gram block g = aᴴ diag(w) b with BLAS-3.
!>
!> @details
!> Row chunks of b are scaled by the weights into a small buffer and
!> accumulated with ZGEMM into one partial matrix per thread (static chunk
!> schedule); the partials are summed in thread order, so the result is
!> deterministic for a given number of threads.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine WeightedGram(g, a, b, wts, nRows, nA, nB)
    !$ use omp_lib, only: omp_get_thread_num, omp_get_max_threads
    use M_Utils_BlasLib, only: ZGEMM

    integer(I32), intent(in) :: nRows, nA, nB
    !> Output block aᴴ diag(w) b.
    complex(R64), intent(out) :: g(nA, nB)
    !> Bra columns (leading dimension nRows).
    complex(R64), intent(in) :: a(nRows, *)
    !> Ket columns (leading dimension nRows).
    complex(R64), intent(in) :: b(nRows, *)
    !> Metric weights.
    real(R64), intent(in) :: wts(nRows)

    complex(R64), allocatable :: buf(:, :), gParts(:, :, :)
    integer(I32) :: iChunk, nChunks, r0, nr, i, j, iThread, nThreads
    logical :: threadsQ

    nChunks = (nRows + CHUNK_ROWS - 1) / CHUNK_ROWS
    threadsQ = ParallelWorkQ(nRows, nA * nB)
    nThreads = 1
    !$ if (threadsQ) nThreads = int(omp_get_max_threads(), I32)

    allocate (gParts(nA, nB, nThreads))
    gParts(:, :, :) = ZERO

    !$omp parallel default(shared) private(buf, iChunk, r0, nr, i, j, iThread) num_threads(nThreads) if (threadsQ)
    iThread = 1
    !$ iThread = int(omp_get_thread_num(), I32) + 1
    allocate (buf(CHUNK_ROWS, nB))
    !$omp do schedule(static)
    do iChunk = 1, nChunks
      r0 = (iChunk - 1) * CHUNK_ROWS + 1
      nr = min(CHUNK_ROWS, nRows - r0 + 1)
      do j = 1, nB
        do i = 1, nr
          buf(i, j) = wts(r0 + i - 1) * b(r0 + i - 1, j)
        end do
      end do
      call ZGEMM('C', 'N', nA, nB, nr, ONE, a(r0, 1), nRows, buf, CHUNK_ROWS, ONE, gParts(1, 1, iThread), nA)
    end do
    !$omp end do
    deallocate (buf)
    !$omp end parallel

    g(:, :) = gParts(:, :, 1)
    do iThread = 2, nThreads
      g(:, :) = g(:, :) + gParts(:, :, iThread)
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Residuals r(:, j) = ax(:, j) − θ_j x(:, j) of the m block columns.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine CalcResiduals(r, x, ax, theta, nRows, m)
    integer(I32), intent(in) :: nRows, m
    !> Output residual columns r(:, 1:m).
    complex(R64), intent(out) :: r(nRows, *)
    !> Ritz vectors x(:, 1:m).
    complex(R64), intent(in) :: x(nRows, *)
    !> Operator applied to the Ritz vectors.
    complex(R64), intent(in) :: ax(nRows, *)
    !> Ritz values.
    real(R64), intent(in) :: theta(m)

    integer(I32) :: i, j
    logical :: threadsQ

    threadsQ = ParallelWorkQ(nRows, m)

    !$omp parallel do default(shared) private(i, j) schedule(static) if (threadsQ)
    do i = 1, nRows
      do j = 1, m
        r(i, j) = ax(i, j) - theta(j) * x(i, j)
      end do
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief In-place right multiplication v(:, 1:kOut) ← v(:, 1:kIn)·t.
!>
!> @details
!> kOut ≤ kIn. Processed in independent row chunks (threaded), each with one
!> ZGEMM into a small buffer that is copied back.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RightMultiplyInPlace(v, t, nRows, kIn, kOut)
    use M_Utils_BlasLib, only: ZGEMM

    integer(I32), intent(in) :: nRows, kIn, kOut
    !> Columns to transform (leading dimension nRows).
    complex(R64), intent(inout) :: v(nRows, *)
    !> Transformation matrix.
    complex(R64), intent(in) :: t(kIn, kOut)

    complex(R64), allocatable :: buf(:, :)
    integer(I32) :: iChunk, nChunks, r0, nr
    logical :: threadsQ

    nChunks = (nRows + CHUNK_ROWS - 1) / CHUNK_ROWS
    threadsQ = ParallelWorkQ(nRows, kIn * kOut)

    !$omp parallel default(shared) private(buf, iChunk, r0, nr) if (threadsQ)
    allocate (buf(CHUNK_ROWS, kOut))
    !$omp do schedule(static)
    do iChunk = 1, nChunks
      r0 = (iChunk - 1) * CHUNK_ROWS + 1
      nr = min(CHUNK_ROWS, nRows - r0 + 1)
      call ZGEMM('N', 'N', nr, kOut, kIn, ONE, v(r0, 1), nRows, t, kIn, ZERO, buf, CHUNK_ROWS)
      v(r0:r0 + nr - 1, 1:kOut) = buf(1:nr, 1:kOut)
    end do
    !$omp end do
    deallocate (buf)
    !$omp end parallel

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief In-place basis update v(:, 1:kOut) ← v(:, 1:kV)·c1 + u(:, 1:kU)·c2.
!>
!> @details
!> kOut ≤ kV + kU, kU ≥ 1. Processed in independent row chunks (threaded).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine CombineInPlace(v, u, c1, c2, nRows, kV, kU, kOut)
    use M_Utils_BlasLib, only: ZGEMM

    integer(I32), intent(in) :: nRows, kV, kU, kOut
    !> First block (leading dimension nRows), overwritten by the result.
    complex(R64), intent(inout) :: v(nRows, *)
    !> Second block (leading dimension nRows).
    complex(R64), intent(in) :: u(nRows, *)
    !> Coefficients of the first block.
    complex(R64), intent(in) :: c1(kV, kOut)
    !> Coefficients of the second block.
    complex(R64), intent(in) :: c2(kU, kOut)

    complex(R64), allocatable :: buf(:, :)
    integer(I32) :: iChunk, nChunks, r0, nr
    logical :: threadsQ

    nChunks = (nRows + CHUNK_ROWS - 1) / CHUNK_ROWS
    threadsQ = ParallelWorkQ(nRows, (kV + kU) * kOut)

    !$omp parallel default(shared) private(buf, iChunk, r0, nr) if (threadsQ)
    allocate (buf(CHUNK_ROWS, kOut))
    !$omp do schedule(static)
    do iChunk = 1, nChunks
      r0 = (iChunk - 1) * CHUNK_ROWS + 1
      nr = min(CHUNK_ROWS, nRows - r0 + 1)
      call ZGEMM('N', 'N', nr, kOut, kV, ONE, v(r0, 1), nRows, c1, kV, ZERO, buf, CHUNK_ROWS)
      call ZGEMM('N', 'N', nr, kOut, kU, ONE, u(r0, 1), nRows, c2, kU, ONE, buf, CHUNK_ROWS)
      v(r0:r0 + nr - 1, 1:kOut) = buf(1:nr, 1:kOut)
    end do
    !$omp end do
    deallocate (buf)
    !$omp end parallel

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Deterministic pseudo-random real columns with entries in (−1, 1).
!>
!> @details
!> Park–Miller minimal standard generator (MINSTD); independent of the
!> global state of `random_number`, so results are reproducible.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillRandom(v, seed)
    !> Output columns.
    complex(R64), intent(out), contiguous :: v(:, :)
    !> Generator state (in [1, MINSTD_M − 1]), advanced.
    integer(I64), intent(inout) :: seed

    integer(I64) :: i, j

    do j = 1, size(v, 2, kind=I64)
      do i = 1, size(v, 1, kind=I64)
        seed = mod(MINSTD_A * seed, MINSTD_M)
        v(i, j) = cmplx(2.0_R64 * real(seed, kind=R64) / real(MINSTD_M, kind=R64) - 1.0_R64, 0.0_R64, kind=R64)
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Whether a vector pass of the given size should be threaded.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function ParallelWorkQ(nRows, nProducts) result(res)
    !$ use omp_lib, only: omp_in_parallel
    use M_Utils_UnusedVariables

    logical :: res
    !> Vector length.
    integer(I32), intent(in) :: nRows
    !> Number of vector products (pairs or columns).
    integer(I32), intent(in) :: nProducts

    if (.false.) call UnusedVariables_Mark(nRows, nProducts)

    res = .false.
    !$ res = (.not. omp_in_parallel()) .and. int(nRows, I64) * int(nProducts, I64) >= MIN_PARALLEL_WORK

  end function

end submodule
