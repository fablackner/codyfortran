! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> EXPOKIT wrappers for matrix exponential actions.
!>
!> Provides typed wrappers around zgexpv/zhexpv/dgexpv/dsexpv to apply exp(tA)
!> to vectors for general/hermitian complex and general/symmetric real cases.
module M_Utils_ExpokitLib
  use M_Utils_Types

  implicit none

  ! Status codes (can be kept if needed, or removed if not used by wrappers)
  integer, parameter, public :: EXPOKIT_SUCCESS = 0
  integer, parameter, public :: EXPOKIT_FAILURE = -1

  ! Abstract interfaces for matrix-vector product
  abstract interface
    subroutine I_ApplyMatOnVecComplex(dState, state)
      import :: R64
      complex(R64), intent(out), contiguous, target :: dState(:)
      complex(R64), intent(in), contiguous, target  :: state(:)
    end subroutine

    subroutine I_ApplyMatOnVecReal(dState, state)
      import :: R64
      real(R64), intent(out), contiguous, target :: dState(:)
      real(R64), intent(in), contiguous, target  :: state(:)
    end subroutine
  end interface

  ! Generic interface for high-level wrappers
  interface ExpokitLib_IntegrateSym
    procedure WrapZhexpv
    procedure WrapDsexpv
  end interface

  ! Generic interface for high-level wrappers
  interface ExpokitLib_IntegrateGeneric
    procedure WrapZgexpv
    procedure WrapDgexpv
  end interface

  ! Private implementation details
  private :: WrapZgexpv, WrapZhexpv, WrapDgexpv, WrapDsexpv
  private :: I_ApplyMatOnVecComplex, I_ApplyMatOnVecReal

  ! Interfaces to Expokit FORTRAN routines using specific kinds
  interface
    subroutine zgexpv(n, m, t, v, w, tol, anorm, wsp, lwsp, iwsp, liwsp, matvec, itrace, iflag)
      import :: I32, R64
      implicit none
      integer(I32), intent(in) :: n                      ! Problem dimension
      integer(I32), intent(in) :: m                      ! Krylov subspace dimension
      real(R64), intent(in) :: t                    ! Time step
      complex(R64), intent(in) :: v(n)              ! Input vector
      complex(R64), intent(out) :: w(n)             ! Output vector
      real(R64), intent(in) :: tol                  ! Error tolerance
      real(R64), intent(in) :: anorm                ! Estimate of 1-norm of A
      complex(R64), intent(inout) :: wsp(*)         ! Complex workspace
      integer(I32), intent(in) :: lwsp                   ! Length of wsp
      integer(I32), intent(inout) :: iwsp(*)             ! Integer workspace
      integer(I32), intent(in) :: liwsp                  ! Length of iwsp
      external :: matvec                            ! Matrix-vector product subroutine
      integer(I32), intent(in) :: itrace                 ! Trace flag (0=silent, 1=trace)
      integer(I32), intent(out) :: iflag                 ! Error flag
    end subroutine
  end interface

  interface
    subroutine zhexpv(n, m, t, v, w, tol, anorm, wsp, lwsp, iwsp, liwsp, matvec, itrace, iflag)
      import :: I32, R64
      implicit none
      integer(I32), intent(in) :: n                      ! Problem dimension
      integer(I32), intent(in) :: m                      ! Krylov subspace dimension
      real(R64), intent(in) :: t                    ! Time step
      complex(R64), intent(in) :: v(n)              ! Input vector
      complex(R64), intent(out) :: w(n)             ! Output vector
      real(R64), intent(in) :: tol                  ! Error tolerance
      real(R64), intent(in) :: anorm                ! Estimate of 1-norm of A
      complex(R64), intent(inout) :: wsp(*)         ! Complex workspace
      integer(I32), intent(in) :: lwsp                   ! Length of wsp
      integer(I32), intent(inout) :: iwsp(*)             ! Integer workspace
      integer(I32), intent(in) :: liwsp                  ! Length of iwsp
      external :: matvec                            ! Matrix-vector product subroutine
      integer(I32), intent(in) :: itrace                 ! Trace flag (0=silent, 1=trace)
      integer(I32), intent(out) :: iflag                 ! Error flag
    end subroutine
  end interface

  interface
    subroutine dgexpv(n, m, t, v, w, tol, anorm, wsp, lwsp, iwsp, liwsp, matvec, itrace, iflag)
      import :: I32, R64
      implicit none
      integer(I32), intent(in) :: n                  ! Problem dimension
      integer(I32), intent(in) :: m                  ! Krylov subspace dimension
      real(R64), intent(in) :: t                ! Time step
      real(R64), intent(in) :: v(n)             ! Input vector
      real(R64), intent(out) :: w(n)            ! Output vector
      real(R64), intent(in) :: tol              ! Error tolerance
      real(R64), intent(in) :: anorm            ! Estimate of 1-norm of A
      real(R64), intent(inout) :: wsp(*)        ! Real workspace
      integer(I32), intent(in) :: lwsp               ! Length of wsp
      integer(I32), intent(inout) :: iwsp(*)         ! Integer workspace
      integer(I32), intent(in) :: liwsp              ! Length of iwsp
      external :: matvec                        ! Matrix-vector product subroutine
      integer(I32), intent(in) :: itrace             ! Trace flag (0=silent, 1=trace)
      integer(I32), intent(out) :: iflag             ! Error flag
    end subroutine
  end interface

  interface
    subroutine dsexpv(n, m, t, v, w, tol, anorm, wsp, lwsp, iwsp, liwsp, matvec, itrace, iflag)
      import :: I32, R64
      implicit none
      integer(I32), intent(in) :: n                  ! Problem dimension
      integer(I32), intent(in) :: m                  ! Krylov subspace dimension
      real(R64), intent(in) :: t                ! Time step
      real(R64), intent(in) :: v(n)             ! Input vector
      real(R64), intent(out) :: w(n)            ! Output vector
      real(R64), intent(in) :: tol              ! Error tolerance
      real(R64), intent(in) :: anorm            ! Estimate of 1-norm of A
      real(R64), intent(inout) :: wsp(*)        ! Real workspace
      integer(I32), intent(in) :: lwsp               ! Length of wsp
      integer(I32), intent(inout) :: iwsp(*)         ! Integer workspace
      integer(I32), intent(in) :: liwsp              ! Length of iwsp
      external :: matvec                        ! Matrix-vector product subroutine
      integer(I32), intent(in) :: itrace             ! Trace flag (0=silent, 1=trace)
      integer(I32), intent(out) :: iflag             ! Error flag
    end subroutine
  end interface

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Wrapper for zgexpv (General Complex)
  subroutine WrapZgexpv(state, t, ApplyMatOnVec, krylovDim_, tol_, nSteps_, iflagOut_)
    complex(R64), intent(inout), contiguous, target :: state(:) ! State vector to integrate
    real(R64), intent(in) :: t                     ! Time step
    procedure(I_ApplyMatOnVecComplex) :: ApplyMatOnVec ! Matrix-vector product
    integer(I32), intent(in), optional :: krylovDim_    ! Krylov subspace dimension
    real(R64), intent(in), optional :: tol_         ! Error tolerance
    integer(I32), intent(out), optional :: nSteps_       ! Number of steps taken
    integer(I32), intent(out), optional :: iflagOut_    ! Error flag output

    ! Local variables
    integer(I32) :: n, m, lWsp, lIwsp, maxSteps, iTrace, iFlag
    real(R64) :: aNorm, errorTol
    complex(R64), allocatable :: wSp(:), resultVec(:)
    integer(I32), allocatable :: iWsp(:)

    ! Set dimensions and parameters
    n = size(state)
    m = 30                       ! Default Krylov dimension
    if (present(krylovDim_)) m = krylovDim_

    errorTol = 1.0e-7_R64       ! Default tolerance
    if (present(tol_)) errorTol = tol_

    maxSteps = 1000             ! Default max steps (adjust as needed)
    iTrace = 0                   ! Silent mode by default

    ! Calculate workspace sizes
    ! Note: Expokit documentation formulas might vary slightly. Verify these.
    lWsp = n * (m + 2) + 5 * (m + 2)**2 + 7 ! Complex workspace size for zgexpv/zhexpv
    lIwsp = m + 2                  ! Integer workspace size (check Expokit docs, might need max_steps)

    allocate (wSp(lWsp), iWsp(lIwsp), resultVec(n))

    ! Estimate matrix norm (simplified - replace with better estimator if needed)
    aNorm = 1.0_R64 ! Placeholder - Expokit_EstimateNorm could be adapted or another used

    ! Call zgexpv to perform the propagation
    call zgexpv(n, m, t, state, resultVec, errorTol, aNorm, wSp, lWsp, &
                iWsp, lIwsp, MatVecWrapper, iTrace, iFlag)

    ! Copy result back to state
    state = resultVec

    ! Return optional outputs
    ! Note: Expokit's iwsp might store different info than ARPACK's iparam.
    ! Check Expokit documentation for where step count is stored if needed.
    if (present(nSteps_)) nSteps_ = -1 ! Placeholder - Update if step count location is known
    if (present(iflagOut_)) iflagOut_ = iFlag

    ! Check for errors reported by Expokit
    if (iFlag .ne. 0) then
      write (*, *) 'Expokit zgexpv failed with iflag = ', iFlag
      ! Consider adding more specific error messages based on iflag values
    end if

    ! Clean up
    deallocate (wSp, iWsp, resultVec)

  contains
    ! Wrapper to adapt our procedure pointer to Expokit's external interface
    subroutine MatVecWrapper(nDimIn, inVec, outVec)
      integer(I32), intent(in) :: nDimIn
      complex(R64), intent(in) :: inVec(nDimIn)
      complex(R64), intent(out) :: outVec(nDimIn)
      ! Ensure contiguous target attribute matches abstract interface
      complex(R64), target :: inVecTarget(nDimIn)
      complex(R64), target :: outVecTarget(nDimIn)
      inVecTarget = inVec
      call ApplyMatOnVec(outVecTarget, inVecTarget)
      outVec = outVecTarget
    end subroutine
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Wrapper for zhexpv (Hermitian Complex)
  subroutine WrapZhexpv(state, t, ApplyMatOnVec, krylovDim_, tol_, nSteps_, iflagOut_)
    complex(R64), intent(inout), contiguous, target :: state(:) ! State vector to integrate
    real(R64), intent(in) :: t                     ! Time step
    procedure(I_ApplyMatOnVecComplex) :: ApplyMatOnVec ! Matrix-vector product
    integer(I32), intent(in), optional :: krylovDim_    ! Krylov subspace dimension
    real(R64), intent(in), optional :: tol_         ! Error tolerance
    integer(I32), intent(out), optional :: nSteps_       ! Number of steps taken
    integer(I32), intent(out), optional :: iflagOut_    ! Error flag output

    ! Local variables (similar to WrapZgexpv)
    integer(I32) :: n, m, lWsp, lIwsp, maxSteps, iTrace, iFlag
    real(R64) :: aNorm, errorTol
    complex(R64), allocatable :: wSp(:), resultVec(:)
    integer(I32), allocatable :: iWsp(:)

    ! Set dimensions and parameters
    n = size(state)
    m = 30                       ! Default Krylov dimension
    if (present(krylovDim_)) m = krylovDim_

    errorTol = 1.0e-7_R64       ! Default tolerance
    if (present(tol_)) errorTol = tol_

    maxSteps = 1000             ! Default max steps
    iTrace = 0                   ! Silent mode

    ! Calculate workspace sizes (same as zgexpv)
    lWsp = n * (m + 2) + 5 * (m + 2)**2 + 7
    lIwsp = m + 2                  ! Check Expokit docs

    allocate (wSp(lWsp), iWsp(lIwsp), resultVec(n))

    ! Estimate matrix norm (simplified)
    aNorm = 1.0_R64 ! Placeholder

    ! Call zhexpv
    call zhexpv(n, m, t, state, resultVec, errorTol, aNorm, wSp, lWsp, &
                iWsp, lIwsp, MatVecWrapper, iTrace, iFlag)

    ! Copy result back
    state = resultVec

    ! Return optional outputs
    if (present(nSteps_)) nSteps_ = -1 ! Placeholder
    if (present(iflagOut_)) iflagOut_ = iFlag

    ! Error check
    if (iFlag .ne. 0) then
      write (*, *) 'Expokit zhexpv failed with iflag = ', iFlag
    end if

    ! Clean up
    deallocate (wSp, iWsp, resultVec)

  contains
    ! Wrapper (identical to the one in WrapZgexpv)
    subroutine MatVecWrapper(nDimIn, inVec, outVec)
      integer(I32), intent(in) :: nDimIn
      complex(R64), intent(in) :: inVec(nDimIn)
      complex(R64), intent(out) :: outVec(nDimIn)
      complex(R64), target :: inVecTarget(nDimIn)
      complex(R64), target :: outVecTarget(nDimIn)
      inVecTarget = inVec
      call ApplyMatOnVec(outVecTarget, inVecTarget)
      outVec = outVecTarget
    end subroutine
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Wrapper for dgexpv (General Real)
  subroutine WrapDgexpv(state, t, ApplyMatOnVec, krylovDim_, tol_, nSteps_, iflagOut_)
    real(R64), intent(inout), contiguous, target :: state(:) ! State vector to integrate
    real(R64), intent(in) :: t                    ! Time step
    procedure(I_ApplyMatOnVecReal) :: ApplyMatOnVec ! Matrix-vector product
    integer(I32), intent(in), optional :: krylovDim_   ! Krylov subspace dimension
    real(R64), intent(in), optional :: tol_        ! Error tolerance
    integer(I32), intent(out), optional :: nSteps_      ! Number of steps taken
    integer(I32), intent(out), optional :: iflagOut_   ! Error flag output

    ! Local variables
    integer(I32) :: n, m, lWsp, lIwsp, maxSteps, iTrace, iFlag
    real(R64) :: aNorm, errorTol
    real(R64), allocatable :: wSp(:), resultVec(:)
    integer(I32), allocatable :: iWsp(:)

    ! Set dimensions and parameters
    n = size(state)
    m = 30                       ! Default Krylov dimension
    if (present(krylovDim_)) m = krylovDim_

    errorTol = 1.0e-7_R64       ! Default tolerance
    if (present(tol_)) errorTol = tol_

    maxSteps = 1000             ! Default max steps
    iTrace = 0                   ! Silent mode

    ! Calculate workspace sizes for dgexpv/dsexpv
    ! Check Expokit documentation for correct formulas
    lWsp = n * (m + 1) + m * m + 5 * m + 7 ! Example formula, verify!
    lIwsp = m + 2                  ! Example formula, verify!

    allocate (wSp(lWsp), iWsp(lIwsp), resultVec(n))

    ! Estimate matrix norm (simplified)
    aNorm = 1.0_R64 ! Placeholder

    ! Call dgexpv
    call dgexpv(n, m, t, state, resultVec, errorTol, aNorm, wSp, lWsp, &
                iWsp, lIwsp, MatVecWrapper, iTrace, iFlag)

    ! Copy result back
    state = resultVec

    ! Return optional outputs
    if (present(nSteps_)) nSteps_ = -1 ! Placeholder
    if (present(iflagOut_)) iflagOut_ = iFlag

    ! Error check
    if (iFlag .ne. 0) then
      write (*, *) 'Expokit dgexpv failed with iflag = ', iFlag
    end if

    ! Clean up
    deallocate (wSp, iWsp, resultVec)

  contains
    ! Wrapper for real matrix-vector product
    subroutine MatVecWrapper(nDimIn, inVec, outVec)
      integer(I32), intent(in) :: nDimIn
      real(R64), intent(in) :: inVec(nDimIn)
      real(R64), intent(out) :: outVec(nDimIn)
      real(R64), target :: inVecTarget(nDimIn)
      real(R64), target :: outVecTarget(nDimIn)
      inVecTarget = inVec
      call ApplyMatOnVec(outVecTarget, inVecTarget)
      outVec = outVecTarget
    end subroutine
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Wrapper for dsexpv (Symmetric Real)
  subroutine WrapDsexpv(state, t, ApplyMatOnVec, krylovDim_, tol_, nSteps_, iflagOut_)
    real(R64), intent(inout), contiguous, target :: state(:) ! State vector to integrate
    real(R64), intent(in) :: t                    ! Time step
    procedure(I_ApplyMatOnVecReal) :: ApplyMatOnVec ! Matrix-vector product
    integer(I32), intent(in), optional :: krylovDim_   ! Krylov subspace dimension
    real(R64), intent(in), optional :: tol_        ! Error tolerance
    integer(I32), intent(out), optional :: nSteps_      ! Number of steps taken
    integer(I32), intent(out), optional :: iflagOut_   ! Error flag output

    ! Local variables (similar to WrapDgexpv)
    integer(I32) :: n, m, lWsp, lIwsp, maxSteps, iTrace, iFlag
    real(R64) :: aNorm, errorTol
    real(R64), allocatable :: wSp(:), resultVec(:)
    integer(I32), allocatable :: iWsp(:)

    ! Set dimensions and parameters
    n = size(state)
    m = 30                       ! Default Krylov dimension
    if (present(krylovDim_)) m = krylovDim_

    errorTol = 1.0e-7_R64       ! Default tolerance
    if (present(tol_)) errorTol = tol_

    maxSteps = 1000             ! Default max steps
    iTrace = 0                   ! Silent mode

    ! Calculate workspace sizes (same as dgexpv - check docs)
    lWsp = n * (m + 1) + m * m + 5 * m + 7 ! Example formula, verify!
    lIwsp = m + 2                  ! Example formula, verify!

    allocate (wSp(lWsp), iWsp(lIwsp), resultVec(n))

    ! Estimate matrix norm (simplified)
    aNorm = 1.0_R64 ! Placeholder

    ! Call dsexpv
    call dsexpv(n, m, t, state, resultVec, errorTol, aNorm, wSp, lWsp, &
                iWsp, lIwsp, MatVecWrapper, iTrace, iFlag)

    ! Copy result back
    state = resultVec

    ! Return optional outputs
    if (present(nSteps_)) nSteps_ = -1 ! Placeholder
    if (present(iflagOut_)) iflagOut_ = iFlag

    ! Error check
    if (iFlag .ne. 0) then
      write (*, *) 'Expokit dsexpv failed with iflag = ', iFlag
    end if

    ! Clean up
    deallocate (wSp, iWsp, resultVec)

  contains
    ! Wrapper (identical to the one in WrapDgexpv)
    subroutine MatVecWrapper(nDimIn, inVec, outVec)
      integer(I32), intent(in) :: nDimIn
      real(R64), intent(in) :: inVec(nDimIn)
      real(R64), intent(out) :: outVec(nDimIn)
      real(R64), target :: inVecTarget(nDimIn)
      real(R64), target :: outVecTarget(nDimIn)
      inVecTarget = inVec
      call ApplyMatOnVec(outVecTarget, inVecTarget)
      outVec = outVecTarget
    end subroutine
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end module
