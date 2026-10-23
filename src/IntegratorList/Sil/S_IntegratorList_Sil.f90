! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file S_IntegratorList_Sil.f90
!> @brief Short Iterative Lanczos (SIL) propagator implementation.
!>
!> @details
!> Approximates exp(−i Δt Ĥ) |ψ⟩ using a short Krylov subspace built via
!> the symmetric Lanczos process. The method is memory-efficient and well-
!> suited for Hermitian generators.
!>
!> ## Algorithm Outline
!>
!> 1. Build an orthonormal Lanczos basis {v₁, v₂, …} with tridiagonal
!>    projection T_j, using full reorthogonalization for numerical
!>    robustness.
!> 2. After each new vector, diagonalize T_j → eigenvalues λ, eigenvectors U
!>    and evaluate the propagated coefficients
!>    φ = U · exp(−i Δt Λ) · U^† · e₁.
!> 3. Stop growing the basis as soon as the residual error estimate
!>    β_j · |φ_j| drops below `tolerance` (adaptive Krylov dimension).
!> 4. Map back: ψ_new = ||ψ|| · V_m · φ.
!>
!> ## Adaptive Sub-Stepping
!>
!> If the error estimate still exceeds `tolerance` at the maximum Krylov
!> dimension, the step is halved and retried up to `maxRestarts` times.
!>
!> ## Configuration (JSON)
!>
!> | Key           | Type    | Default | Description                        |
!> |---------------|---------|---------|------------------------------------|
!> | `krylovDim`   | integer | 32      | Maximum Krylov subspace dimension  |
!> | `tolerance`   | real    | 1e-10   | Error tolerance per sub-step       |
!> | `maxRestarts` | integer | 8       | Max step-halving retries           |
!> | `imagTimeQ`   | logical | false   | Imaginary time: exp(−Δτ Ĥ)         |
!>
!> With `imagTimeQ` the same Lanczos basis yields exp(−Δτ Ĥ)|ψ⟩ (the callback
!> keeps the real-time convention −i Ĥ v). The result is not normalized; the
!> caller renormalizes (e.g., orbital orthonormalization) as usual.
!>
!> ## Properties
!>
!> - **Accuracy**: Controlled by `tolerance`; the Krylov dimension and the
!>   sub-step size adapt automatically.
!> - **Stability**: Exact for time-independent Hermitian Hamiltonians; for
!>   time-dependent generators the derivative is evaluated at the sub-step
!>   midpoint (2nd-order accurate).
!> - **Cost**: One `TimeDerivative` evaluation per Lanczos vector.
!>
!> ## Metric
!>
!> The Lanczos recursion assumes the generator is self-adjoint with respect
!> to the inner product in use. On grids whose representation is not
!> weight-absorbed (e.g., FEDVR), operators are self-adjoint only in the
!> weighted metric; assign `metricWeights` on the integrator element to make
!> all inner products (norms, alpha/beta coefficients, reorthogonalization)
!> use sum(conjg(a) * w * b). Unallocated weights give the Euclidean product.
submodule(M_IntegratorList_Sil) S_IntegratorList_Sil

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine IntegratorList_Sil_Allocate(e, path)
    use M_Utils_UnusedVariables
    class(T_IntegratorList_E), allocatable, intent(out) :: e
    character(len=*), intent(in) :: path

    if (.false.) call UnusedVariables_Mark(path)

    allocate (T_IntegratorList_E_Sil :: e)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Fabricate(this)
    use M_Utils_Json
    use M_Utils_Say

    class(T_IntegratorList_E_Sil), intent(inout) :: this

    call Say_Fabricate(this % path//".sil")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    this % krylovDim = Json_Get("krylovDim", 32, path_=this % path)
    this % tolerance = Json_Get("tolerance", 1.0e-10_R64, path_=this % path)
    this % maxRestarts = Json_Get("maxRestarts", 8, path_=this % path)
    this % imagTimeQ = Json_Get("imagTimeQ", .false., path_=this % path)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Setup(this)
    use M_Utils_Say

    class(T_IntegratorList_E_Sil), intent(inout) :: this

    call Say_Setup(this % path//".sil")

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Advance the state using the Short Iterative Lanczos method.
  !>
  !> Builds a Krylov subspace via symmetric Lanczos with full
  !> reorthogonalization, growing the subspace only until the residual error
  !> estimate drops below tolerance. If the maximum Krylov dimension is
  !> insufficient, the sub-step is halved and retried.
  !>
  !> The Krylov basis and the work vectors persist in the integrator element
  !> (reallocated only when the state size changes). The vector operations
  !> are fused OpenMP loops, the projections and the final combination BLAS
  !> calls, so the cost per Lanczos vector beyond the TimeDerivative call is
  !> a few passes over the state. TimeDerivative must depend on its argument
  !> only: the state is updated after every sub-step.
  !>
  !> @param[inout] this   The SIL integrator instance.
  !> @param[inout] state  Complex state vector (modified in-place).
  !> @param[in]    t0     Starting time.
  !> @param[in]    t1     Target ending time.
  module subroutine Integrate(this, state, t0, t1)
    use M_Utils_BlasLib, only: ZGEMV

    class(T_IntegratorList_E_Sil), intent(inout) :: this
    complex(R64), intent(inout), contiguous :: state(:)
    real(R64), intent(in) :: t0
    real(R64), intent(in) :: t1

    real(R64), parameter :: zeroThreshold = 1.0e-14_R64
    ! states below this size are not worth the OpenMP fork/join
    integer, parameter :: minParallelSize = 8192
    complex(R64), parameter :: one = (1.0_R64, 0.0_R64)
    complex(R64), parameter :: zero = (0.0_R64, 0.0_R64)

    complex(R64), allocatable :: phi(:)
    real(R64), allocatable :: alpha(:), beta(:)
    integer :: nState, maxDim

    if (abs(t1 - t0) <= zeroThreshold) return

    nState = size(state)
    maxDim = max(1, min(this % krylovDim, nState))

    if (allocated(this % basis)) then
      if (size(this % basis, 1) .ne. nState .or. size(this % basis, 2) .ne. maxDim) &
        deallocate (this % basis, this % w, this % deriv, this % wMetric)
    end if
    if (.not. allocated(this % basis)) &
      allocate (this % basis(nState, maxDim), this % w(nState), this % deriv(nState), this % wMetric(nState))

    allocate (phi(maxDim), alpha(maxDim), beta(maxDim))

    if (allocated(this % metricWeights)) then
      if (size(this % metricWeights) .ne. nState) then
        error stop "IntegratorList_Sil: metricWeights size does not match state size"
      end if
      call Propagate(this % metricWeights)
    else
      if (allocated(this % unitWeights)) then
        if (size(this % unitWeights) .ne. nState) deallocate (this % unitWeights)
      end if
      if (.not. allocated(this % unitWeights)) allocate (this % unitWeights(nState), source=1.0_R64)
      call Propagate(this % unitWeights)
    end if

  contains

    !> Sub-step loop from t0 to t1 in the metric with weights wts.
    subroutine Propagate(wts)
      real(R64), intent(in), contiguous :: wts(:)

      real(R64) :: dtRemaining, dtStep, tCurrent, errorEstimate, stateNorm
      integer :: mDim, restartCount
      logical :: stepConverged

      tCurrent = t0
      dtRemaining = t1 - t0

      do while (abs(dtRemaining) > zeroThreshold)
        dtStep = dtRemaining
        restartCount = 0

        stateNorm = sqrt(MetricNormSq(state, wts))
        if (stateNorm <= zeroThreshold) return
        call ScaleInto(this % basis(:, 1), state, 1.0_R64 / stateNorm)

        do
          call PerformLanczosStep(mDim, errorEstimate, stepConverged, tCurrent, dtStep, wts)
          if (stepConverged) exit
          if (abs(dtStep) <= zeroThreshold .or. restartCount >= this % maxRestarts) exit
          dtStep = 0.5_R64 * dtStep
          restartCount = restartCount + 1
        end do

        if (.not. stepConverged) then
          write (*, '(A,ES13.6,A,ES13.6)') "Warning: Short iterative Lanczos step starting at t=", tCurrent, &
            " exceeded tolerance. Estimated error=", errorEstimate
        end if

        ! state = ||state|| V_m phi
        call ZGEMV('N', nState, mDim, cmplx(stateNorm, 0.0_R64, kind=R64), this % basis, nState, phi, 1, zero, state, 1)

        tCurrent = tCurrent + dtStep
        dtRemaining = dtRemaining - dtStep
      end do

    end subroutine

    !> Perform one Lanczos sub-step from v_1 = basis(:, 1): grow the Krylov
    !> basis until the error estimate is below tolerance and leave the
    !> propagated coefficients in phi(1:mDimOut).
    !>
    !> @param[out] mDimOut    Krylov dimension used.
    !> @param[out] errEst     Estimated truncation error.
    !> @param[out] converged  True if the error is below tolerance.
    !> @param[in]  tBase      Sub-step starting time.
    !> @param[in]  dt         Time step to propagate.
    !> @param[in]  wts        Metric weights.
    subroutine PerformLanczosStep(mDimOut, errEst, converged, tBase, dt, wts)
      integer, intent(out) :: mDimOut
      real(R64), intent(out) :: errEst
      logical, intent(out) :: converged
      real(R64), intent(in) :: tBase
      real(R64), intent(in) :: dt
      real(R64), intent(in), contiguous :: wts(:)

      real(R64) :: timeEval, betaPrev
      integer :: j
      logical :: breakdown

      timeEval = tBase + 0.5_R64 * dt
      mDimOut = maxDim

      do j = 1, maxDim
        ! w = Ĥ v_j - beta_(j-1) v_(j-1) (TimeDerivative supplies −i Ĥ v_j)
        call this % TimeDerivative(this % deriv, this % basis(:, j), timeEval)
        betaPrev = 0.0_R64
        if (j > 1) betaPrev = beta(j - 1)
        alpha(j) = RecurrenceStep(j, betaPrev, wts)

        ! full reorthogonalization against all previous Lanczos vectors; after
        ! the three-term recurrence w is orthogonal up to round-off, so one
        ! classical Gram-Schmidt pass (two BLAS sweeps over the basis instead
        ! of j sweeps) is as accurate as the modified variant
        call ZGEMV('C', nState, j, one, this % basis, nState, this % wMetric, 1, zero, phi, 1)
        call ZGEMV('N', nState, j, -one, this % basis, nState, phi, 1, one, this % w, 1)

        beta(j) = sqrt(MetricNormSq(this % w, wts))
        breakdown = beta(j) <= zeroThreshold

        call EvaluateExpInSubspace(phi(1:j), errEst, dt, j)
        if (breakdown) errEst = 0.0_R64
        converged = errEst <= this % tolerance

        if (converged .or. breakdown .or. j .eq. maxDim) then
          mDimOut = j
          exit
        end if

        call ScaleInto(this % basis(:, j + 1), this % w, 1.0_R64 / beta(j))
      end do

    end subroutine

    !> One fused pass for the three-term recurrence of Lanczos vector j:
    !> w = i deriv - betaPrev v_(j-1), alpha = Re<v_j|w>, w = w - alpha v_j,
    !> and wMetric = W w for the reorthogonalization (W the metric weights).
    function RecurrenceStep(j, betaPrev, wts) result(alphaJ)
      integer, intent(in) :: j
      real(R64), intent(in) :: betaPrev
      real(R64), intent(in), contiguous :: wts(:)
      real(R64) :: alphaJ

      integer :: jPrev, i

      jPrev = max(j - 1, 1)
      alphaJ = 0.0_R64
      !$omp parallel do default(shared) private(i) reduction(+:alphaJ) schedule(static) if (nState >= minParallelSize)
      do i = 1, nState
        this % w(i) = cmplx(-aimag(this % deriv(i)), real(this % deriv(i), kind=R64), kind=R64) &
                      - betaPrev * this % basis(i, jPrev)
        alphaJ = alphaJ + wts(i) * real(conjg(this % basis(i, j)) * this % w(i), kind=R64)
      end do
      !$omp end parallel do
      !$omp parallel do default(shared) private(i) schedule(static) if (nState >= minParallelSize)
      do i = 1, nState
        this % w(i) = this % w(i) - alphaJ * this % basis(i, j)
        this % wMetric(i) = wts(i) * this % w(i)
      end do
      !$omp end parallel do

    end function

    !> Squared norm of v in the metric with weights wts.
    function MetricNormSq(v, wts) result(res)
      complex(R64), intent(in), contiguous :: v(:)
      real(R64), intent(in), contiguous :: wts(:)
      real(R64) :: res

      integer :: i

      res = 0.0_R64
      !$omp parallel do default(shared) private(i) reduction(+:res) schedule(static) if (nState >= minParallelSize)
      do i = 1, nState
        res = res + wts(i) * (real(v(i), kind=R64)**2 + aimag(v(i))**2)
      end do
      !$omp end parallel do

    end function

    !> y = s x.
    subroutine ScaleInto(y, x, s)
      complex(R64), intent(out), contiguous :: y(:)
      complex(R64), intent(in), contiguous :: x(:)
      real(R64), intent(in) :: s

      integer :: i

      !$omp parallel do default(shared) private(i) schedule(static) if (nState >= minParallelSize)
      do i = 1, nState
        y(i) = s * x(i)
      end do
      !$omp end parallel do

    end subroutine

    !> Evaluate φ = U exp(−i dt Λ) U^† e₁ in the current Krylov subspace and
    !> the associated residual error estimate β_m · |φ_m|.
    !>
    !> @param[out] phiOut  Propagated coefficients in the Krylov basis.
    !> @param[out] errEst  Residual-based truncation error estimate.
    !> @param[in]  dt      Time step to propagate.
    !> @param[in]  mDim Current Krylov subspace dimension.
    subroutine EvaluateExpInSubspace(phiOut, errEst, dt, mDim)
      use M_Utils_LapackLib, only: LapackLib_DiagonalizeGeneric

      complex(R64), intent(out) :: phiOut(:)
      real(R64), intent(out) :: errEst
      real(R64), intent(in) :: dt
      integer, intent(in) :: mDim

      complex(R64), parameter :: eye = (0.0_R64, 1.0_R64)

      real(R64), allocatable :: triDiag(:, :)
      real(R64), allocatable :: eigVals(:), eigVecs(:, :)
      complex(R64) :: coeffs(mDim)
      integer :: j

      allocate (triDiag(mDim, mDim))
      triDiag(:, :) = 0.0_R64

      do j = 1, mDim
        triDiag(j, j) = alpha(j)
      end do
      do j = 1, mDim - 1
        triDiag(j, j + 1) = beta(j)
        triDiag(j + 1, j) = beta(j)
      end do

      call LapackLib_DiagonalizeGeneric(eigVals, eigVecs, triDiag, .true.)

      if (this % imagTimeQ) then
        do j = 1, mDim
          coeffs(j) = eigVecs(1, j) * exp(-dt * eigVals(j))
        end do
      else
        do j = 1, mDim
          coeffs(j) = eigVecs(1, j) * exp(-eye * dt * eigVals(j))
        end do
      end if
      phiOut(:) = matmul(eigVecs, coeffs)

      errEst = beta(mDim) * abs(phiOut(mDim))

    end subroutine

  end subroutine

end submodule
