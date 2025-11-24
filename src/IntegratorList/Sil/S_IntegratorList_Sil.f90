! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Setup(this)
    use M_Utils_Say
    use M_Utils_UnusedVariables

    class(T_IntegratorList_E_Sil), intent(inout) :: this

    call Say_Setup(this % path//".sil")

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Integrate(this, state, t0, t1)
    use M_Utils_LapackLib, only: LapackLib_DiagonalizeGeneric
    use M_Utils_UnusedVariables

    class(T_IntegratorList_E_Sil), intent(inout) :: this
    complex(R64), intent(inout), contiguous :: state(:)
    real(R64), intent(in) :: t0       ! Starting time
    real(R64), intent(in) :: t1       ! Target ending time

    real(R64), parameter :: tolerance = 1.0e-10_R64
    real(R64), parameter :: zeroThreshold = 1.0e-14_R64
    integer, parameter :: maxKrylovDimDefault = 32
    integer, parameter :: maxRestarts = 8
    logical, parameter :: debugMode = .true.

    complex(R64), allocatable :: workState(:), trialState(:)
    real(R64) :: dtTotal, dtRemaining, dtStep, tCurrent, errorEstimate
    integer :: restartCount
    integer :: attemptIndex
    logical :: stepConverged

    if (.false.) call UnusedVariables_Mark(this)

    dtTotal = t1 - t0

    if (debugMode) then
      write (*, '(A)') '=== Debug: Sil Integrate ==='
      write (*, '(A,I0)') '  state size           = ', size(state)
      write (*, '(A,ES13.6)') '  initial time         = ', t0
      write (*, '(A,ES13.6)') '  target time          = ', t1
      write (*, '(A,ES13.6)') '  total dt             = ', dtTotal
    end if

    if (abs(dtTotal) <= zeroThreshold) then
      if (debugMode) then
        write (*, '(A)') '  Early exit: |total dt| <= zeroThreshold'
        write (*, '(A)') '=== End Debug: Sil Integrate ==='
      end if
      return
    end if

    allocate (workState, source=state)
    allocate (trialState, mold=state)

    tCurrent = t0
    dtRemaining = dtTotal

    do while (abs(dtRemaining) > zeroThreshold)
      dtStep = dtRemaining
      restartCount = 0

      if (debugMode) then
        write (*, '("  Step start: t=",ES13.6,", dtRemaining=",ES13.6,", dtStep=",ES13.6)') &
          tCurrent, dtRemaining, dtStep
      end if

      do
        attemptIndex = restartCount + 1
        call PerformLanczosStep(workState, tCurrent, dtStep, trialState, errorEstimate, stepConverged)
        if (debugMode) then
          write (*, '("    Attempt ",I0,": dt=",ES13.6,", converged=",L1,", err_est=",ES13.6)') &
            attemptIndex, dtStep, stepConverged, errorEstimate
        end if
        if (stepConverged) exit
        if (abs(dtStep) <= zeroThreshold .or. restartCount >= maxRestarts) exit
        dtStep = 0.5_R64 * dtStep
        restartCount = restartCount + 1
        if (debugMode) then
          write (*, '("    Reducing step: new dt=",ES13.6,", restartCount=",I0)') dtStep, restartCount
        end if
      end do

      if (.not. stepConverged .and. errorEstimate > tolerance) then
        write (*, '(A,ES13.6,A,ES13.6)') "Warning: Short iterative Lanczos step starting at t=", tCurrent, &
          " exceeded tolerance. Estimated error=", errorEstimate
      end if

      workState = trialState
      if (debugMode) then
        write (*, '("  Step finish: dtUsed=",ES13.6,", tNext=",ES13.6,", restarts=",I0,", err_est=",ES13.6)') &
          dtStep, tCurrent + dtStep, restartCount, errorEstimate
      end if
      tCurrent = tCurrent + dtStep
      dtRemaining = dtRemaining - dtStep
    end do

    if (debugMode) then
      write (*, '(A,ES13.6)') '  final time           = ', tCurrent
      write (*, '(A)') '=== End Debug: Sil Integrate ==='
    end if

    state = workState

    deallocate (workState)
    deallocate (trialState)

  contains

    subroutine PerformLanczosStep(stateIn, tBase, dt, stateOut, err_est, converged)
      complex(R64), intent(in), contiguous :: stateIn(:)
      real(R64), intent(in) :: tBase
      real(R64), intent(in) :: dt
      complex(R64), intent(out), contiguous :: stateOut(:)
      real(R64), intent(out) :: err_est
      logical, intent(out) :: converged

      integer :: n, maxDim, j, mDim
      real(R64) :: stateNorm, betaTail, vectorNorm
      real(R64) :: derivNorm, residualNorm, stateOutNorm
      real(R64), allocatable :: alpha(:), beta(:), triDiag(:, :)
      real(R64), allocatable :: eigVals(:), eigVecs(:, :)
      complex(R64), allocatable, target :: basis(:, :)
      complex(R64), allocatable, target :: deriv(:)
      complex(R64), allocatable :: w(:)
      complex(R64), allocatable :: phi(:), coeffs(:)
      real(R64) :: timeEval
      complex(R64), pointer :: vCurrent(:)
      complex(R64), pointer :: vPrev(:)
      real(R64), parameter :: orthoTolerance = 1.0e-12_R64
      complex(R64), parameter :: eye = (0.0_R64, 1.0_R64)

      if (abs(dt) <= zeroThreshold) then
        stateOut = stateIn
        err_est = 0.0_R64
        converged = .true.
        if (debugMode) then
          write (*, '(A)') '    [Lanczos] skipped: |dt| <= zeroThreshold'
        end if
        return
      end if

      n = size(stateIn)
      maxDim = min(maxKrylovDimDefault, n)

      stateNorm = sqrt(max(real(sum(conjg(stateIn) * stateIn), kind=R64), 0.0_R64))
      if (stateNorm <= zeroThreshold) then
        stateOut = stateIn
        err_est = 0.0_R64
        converged = .true.
        if (debugMode) then
          write (*, '(A)') '    [Lanczos] skipped: state norm <= zeroThreshold'
        end if
        return
      end if

      allocate (basis(n, maxDim))
      allocate (deriv(n))
      allocate (w(n))
      allocate (alpha(maxDim))
      if (maxDim > 1) allocate (beta(maxDim - 1))

      basis(:, :) = (0.0_R64, 0.0_R64)
      deriv(:) = (0.0_R64, 0.0_R64)
      w(:) = (0.0_R64, 0.0_R64)
      alpha(:) = 0.0_R64
      if (allocated(beta)) beta(:) = 0.0_R64

      basis(:, 1) = stateIn / stateNorm
      timeEval = tBase + 0.5_R64 * dt
      betaTail = 0.0_R64
      mDim = maxDim

      if (debugMode) then
        write (*, '("    [Lanczos] n=",I0,", dt=",ES13.6,", maxDim=",I0,", timeEval=",ES13.6)') &
          n, dt, maxDim, timeEval
      end if

      do j = 1, maxDim
        vCurrent(1:n) => basis(:, j)

        call this % TimeDerivative(deriv, vCurrent, timeEval)
        if (debugMode) then
          derivNorm = sqrt(max(real(sum(conjg(deriv) * deriv), kind=R64), 0.0_R64))
          write (*, '("      deriv norm(",I0,")=",ES13.6)') j, derivNorm
        end if
        w(:) = eye * deriv(:)

        if (j > 1) then
          vPrev(1:n) => basis(:, j - 1)
          w(:) = w(:) - beta(j - 1) * vPrev(:)
        end if

        alpha(j) = real(sum(conjg(vCurrent) * w(:)), kind=R64)
        if (debugMode) then
          write (*, '("      alpha(",I0,")=",ES13.6)') j, alpha(j)
        end if
        w(:) = w(:) - alpha(j) * vCurrent(:)
        if (debugMode) then
          residualNorm = sqrt(max(real(sum(conjg(w) * w), kind=R64), 0.0_R64))
          write (*, '("      residual norm(",I0,")=",ES13.6)') j, residualNorm
        end if

        if (j .eq. maxDim) then
          vectorNorm = sqrt(max(real(sum(conjg(w) * w), kind=R64), 0.0_R64))
          betaTail = vectorNorm
          mDim = j
          if (debugMode) then
            write (*, '("      betaTail=",ES13.6," (maxDim reached)")') betaTail
          end if
          exit
        end if

        vectorNorm = sqrt(max(real(sum(conjg(w) * w), kind=R64), 0.0_R64))
        beta(j) = vectorNorm
        if (debugMode) then
          write (*, '("      beta(",I0,")=",ES13.6)') j, beta(j)
        end if

        if (vectorNorm <= orthoTolerance) then
          mDim = j
          if (debugMode) then
            write (*, '("      early termination: beta(",I0,")=",ES13.6)') j, vectorNorm
          end if
          exit
        end if

        basis(:, j + 1) = w(:) / beta(j)
      end do

      allocate (triDiag(mDim, mDim))
      triDiag(:, :) = 0.0_R64
      triDiag(1:mDim, 1:mDim) = 0.0_R64

      do j = 1, mDim
        triDiag(j, j) = alpha(j)
      end do

      if (mDim > 1) then
        do j = 1, mDim - 1
          triDiag(j, j + 1) = beta(j)
          triDiag(j + 1, j) = beta(j)
        end do
      end if

      allocate (eigVals(mDim))
      allocate (eigVecs(mDim, mDim))
      call LapackLib_DiagonalizeGeneric(eigVals, eigVecs, triDiag, .true.)

      allocate (coeffs(mDim))
      allocate (phi(mDim))
      do j = 1, mDim
        coeffs(j) = cmplx(eigVecs(1, j), 0.0_R64, kind=R64) * exp(-eye * dt * eigVals(j))
      end do

      phi(:) = matmul(eigVecs(:, 1:mDim), coeffs(:))

      stateOut(:) = stateNorm * matmul(basis(:, 1:mDim), phi(:))

      if (mDim .eq. 1) then
        err_est = 0.0_R64
      else if (mDim .eq. maxDim) then
        err_est = betaTail * abs(phi(mDim))
      else
        err_est = beta(mDim) * abs(phi(mDim))
      end if

      converged = (err_est <= tolerance)

      if (debugMode) then
        stateOutNorm = sqrt(max(real(sum(conjg(stateOut) * stateOut), kind=R64), 0.0_R64))
        write (*, '("    [Lanczos] mDim=",I0,", err_est=",ES13.6,", converged=",L1)') mDim, err_est, converged
        write (*, '("    [Lanczos] stateOut norm=",ES13.6)') stateOutNorm
      end if

      if (allocated(beta)) deallocate (beta)
      deallocate (phi)
      deallocate (coeffs)
      deallocate (eigVals)
      deallocate (eigVecs)
      deallocate (triDiag)
      deallocate (w)
      deallocate (deriv)
      deallocate (basis)
      deallocate (alpha)

    end subroutine

  end subroutine

end submodule
