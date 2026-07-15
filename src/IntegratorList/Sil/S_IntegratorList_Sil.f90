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
!>
!> ## Properties
!>
!> - **Accuracy**: Controlled by `tolerance`; the Krylov dimension and the
!>   sub-step size adapt automatically.
!> - **Stability**: Exact for time-independent Hermitian Hamiltonians; for
!>   time-dependent generators the derivative is evaluated at the sub-step
!>   midpoint (2nd-order accurate).
!> - **Cost**: One `TimeDerivative` evaluation per Lanczos vector.
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
  !> @param[inout] this   The SIL integrator instance.
  !> @param[inout] state  Complex state vector (modified in-place).
  !> @param[in]    t0     Starting time.
  !> @param[in]    t1     Target ending time.
  module subroutine Integrate(this, state, t0, t1)
    use M_Utils_LapackLib, only: LapackLib_DiagonalizeGeneric

    class(T_IntegratorList_E_Sil), intent(inout) :: this
    complex(R64), intent(inout), contiguous :: state(:)
    real(R64), intent(in) :: t0
    real(R64), intent(in) :: t1

    real(R64), parameter :: zeroThreshold = 1.0e-14_R64

    complex(R64), allocatable :: workState(:), trialState(:)
    complex(R64), allocatable, target :: basis(:, :)
    complex(R64), allocatable :: deriv(:), w(:), phi(:)
    real(R64), allocatable :: alpha(:), beta(:)
    real(R64) :: dtRemaining, dtStep, tCurrent, errorEstimate
    integer :: nState, maxDim, restartCount
    logical :: stepConverged

    if (abs(t1 - t0) <= zeroThreshold) return

    nState = size(state)
    maxDim = max(1, min(this % krylovDim, nState))

    allocate (workState, source=state)
    allocate (trialState, mold=state)
    allocate (basis(nState, maxDim))
    allocate (deriv(nState))
    allocate (w(nState))
    allocate (phi(maxDim))
    allocate (alpha(maxDim))
    allocate (beta(maxDim))

    tCurrent = t0
    dtRemaining = t1 - t0

    do while (abs(dtRemaining) > zeroThreshold)
      dtStep = dtRemaining
      restartCount = 0

      do
        call PerformLanczosStep(trialState, errorEstimate, stepConverged, workState, tCurrent, dtStep)
        if (stepConverged) exit
        if (abs(dtStep) <= zeroThreshold .or. restartCount >= this % maxRestarts) exit
        dtStep = 0.5_R64 * dtStep
        restartCount = restartCount + 1
      end do

      if (.not. stepConverged) then
        write (*, '(A,ES13.6,A,ES13.6)') "Warning: Short iterative Lanczos step starting at t=", tCurrent, &
          " exceeded tolerance. Estimated error=", errorEstimate
      end if

      workState = trialState
      tCurrent = tCurrent + dtStep
      dtRemaining = dtRemaining - dtStep
    end do

    state = workState

  contains

    !> Perform one Lanczos sub-step: grow the Krylov basis until the error
    !> estimate is below tolerance, then apply the exponential in the
    !> reduced basis.
    !>
    !> @param[out] stateOut   Output state after propagation.
    !> @param[out] errEst     Estimated truncation error.
    !> @param[out] converged  True if the error is below tolerance.
    !> @param[in]  stateIn    Input state vector.
    !> @param[in]  tBase      Sub-step starting time.
    !> @param[in]  dt         Time step to propagate.
    subroutine PerformLanczosStep(stateOut, errEst, converged, stateIn, tBase, dt)
      complex(R64), intent(out), contiguous :: stateOut(:)
      real(R64), intent(out) :: errEst
      logical, intent(out) :: converged
      complex(R64), intent(in), contiguous :: stateIn(:)
      real(R64), intent(in) :: tBase
      real(R64), intent(in) :: dt

      complex(R64), parameter :: eye = (0.0_R64, 1.0_R64)

      complex(R64), pointer :: vCurrent(:)
      real(R64) :: stateNorm, timeEval
      integer :: j, k, mDim
      logical :: breakdown

      stateNorm = sqrt(real(dot_product(stateIn, stateIn), kind=R64))
      if (stateNorm <= zeroThreshold) then
        stateOut = stateIn
        errEst = 0.0_R64
        converged = .true.
        return
      end if

      basis(:, 1) = stateIn / stateNorm
      timeEval = tBase + 0.5_R64 * dt

      mDim = maxDim

      do j = 1, maxDim
        vCurrent(1:nState) => basis(:, j)

        ! w = Ĥ v_j (TimeDerivative supplies −i Ĥ v_j)
        call this % TimeDerivative(deriv, vCurrent, timeEval)
        w(:) = eye * deriv(:)

        if (j > 1) w(:) = w(:) - beta(j - 1) * basis(:, j - 1)
        alpha(j) = real(dot_product(vCurrent, w), kind=R64)
        w(:) = w(:) - alpha(j) * vCurrent(:)

        ! full reorthogonalization against all previous Lanczos vectors
        do k = 1, j
          w(:) = w(:) - dot_product(basis(:, k), w) * basis(:, k)
        end do

        beta(j) = sqrt(real(dot_product(w, w), kind=R64))
        breakdown = beta(j) <= zeroThreshold

        call EvaluateExpInSubspace(phi(1:j), errEst, dt, j)
        if (breakdown) errEst = 0.0_R64
        converged = errEst <= this % tolerance

        if (converged .or. breakdown .or. j .eq. maxDim) then
          mDim = j
          exit
        end if

        basis(:, j + 1) = w(:) / beta(j)
      end do

      stateOut(:) = stateNorm * matmul(basis(:, 1:mDim), phi(1:mDim))

    end subroutine

    !> Evaluate φ = U exp(−i dt Λ) U^† e₁ in the current Krylov subspace and
    !> the associated residual error estimate β_m · |φ_m|.
    !>
    !> @param[out] phiOut  Propagated coefficients in the Krylov basis.
    !> @param[out] errEst  Residual-based truncation error estimate.
    !> @param[in]  dt      Time step to propagate.
    !> @param[in]  mDim    Current Krylov subspace dimension.
    subroutine EvaluateExpInSubspace(phiOut, errEst, dt, mDim)
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

      do j = 1, mDim
        coeffs(j) = eigVecs(1, j) * exp(-eye * dt * eigVals(j))
      end do
      phiOut(:) = matmul(eigVecs, coeffs)

      errEst = beta(mDim) * abs(phiOut(mDim))

    end subroutine

  end subroutine

end submodule
