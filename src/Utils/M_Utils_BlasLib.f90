! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> BLAS helper interfaces for typed linear algebra operations.
!>
!> Provides Fortran interfaces and convenience wrappers for common BLAS level
!> 1/2/3 routines in double precision real/complex, simplifying calls and
!> enabling allocatable output variants used across the codebase.
module M_Utils_BlasLib
  use M_Utils_Types

  implicit none

  interface BlasLib_CalcNorm
    procedure WrapDZNRM2
    procedure WrapDNRM2
  end interface

  interface BlasLib_CalcMatOnMat
    procedure WrapZGEMM
    procedure WrapDGEMM
    procedure WrapAllocatableZGEMM
    procedure WrapAllocatableDGEMM
  end interface

  interface BlasLib_CalcMatOnVec
    procedure WrapZGEMV
    procedure WrapDGEMV
    procedure WrapAllocatableZGEMV
    procedure WrapAllocatableDGEMV
  end interface

  interface

    !> Euclidean 2-norm of a complex vector.
    real(R64) function DZNRM2(N, X, INCX)
      import :: I32, R64
      !> Number of elements in the vector X.
      integer(I32), intent(in) :: N
      !> Input vector X. Length at least 1+(N-1)*abs(INCX).
      complex(R64), intent(in) :: X(*)
      !> Storage increment for elements of X.
      integer(I32), intent(in) :: INCX
    end function

    !> Euclidean 2-norm of a real vector.
    real(R64) function DNRM2(N, X, INCX)
      import :: I32, R64
      !> Number of elements in the vector X.
      integer(I32), intent(in) :: N
      !> Input vector X. Length at least 1+(N-1)*abs(INCX).
      real(R64), intent(in) :: X(*)
      !> Storage increment for elements of X.
      integer(I32), intent(in) :: INCX
    end function

  end interface

  interface

    !> BLAS ZGEMM: matrix-matrix multiplication for complex numbers.
    subroutine ZGEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
      import :: I32, R64
      !> Operation flag for A: 'N' (no transpose), 'T' (transpose), 'C' (conjugate transpose).
      character(len=1), intent(in) :: TRANSA
      !> Operation flag for B: 'N' (no transpose), 'T' (transpose), 'C' (conjugate transpose).
      character(len=1), intent(in) :: TRANSB
      !> Number of rows of op(A) and of C.
      integer(I32), intent(in) :: M
      !> Number of columns of op(B) and of C.
      integer(I32), intent(in) :: N
      !> Number of columns of op(A) and rows of op(B).
      integer(I32), intent(in) :: K
      !> Scalar multiplier for op(A)*op(B).
      complex(R64), intent(in) :: ALPHA
      !> Leading dimension (row stride) of A.
      integer(I32), intent(in) :: LDA
      !> Input matrix A, size at least (LDA, ka) with ka = K if TRANSA='N', else ka = M.
      complex(R64), intent(in) :: A(LDA, *)
      !> Leading dimension (row stride) of B.
      integer(I32), intent(in) :: LDB
      !> Input matrix B, size at least (LDB, kb) with kb = N if TRANSB='N', else kb = K.
      complex(R64), intent(in) :: B(LDB, *)
      !> Scalar multiplier for C.
      complex(R64), intent(in) :: BETA
      !> Leading dimension (row stride) of C.
      integer(I32), intent(in) :: LDC
      !> In-place output matrix C on entry, alpha*op(A)*op(B) + beta*C on exit.
      complex(R64), intent(inout) :: C(LDC, *)
    end subroutine

    !> BLAS DGEMM: matrix-matrix multiplication for real numbers.
    subroutine DGEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
      import :: I32, R64
      !> Operation flag for A: 'N' (no transpose), 'T' (transpose).
      character(len=1), intent(in) :: TRANSA
      !> Operation flag for B: 'N' (no transpose), 'T' (transpose).
      character(len=1), intent(in) :: TRANSB
      !> Number of rows of op(A) and of C.
      integer(I32), intent(in) :: M
      !> Number of columns of op(B) and of C.
      integer(I32), intent(in) :: N
      !> Number of columns of op(A) and rows of op(B).
      integer(I32), intent(in) :: K
      !> Scalar multiplier for op(A)*op(B).
      real(R64), intent(in) :: ALPHA
      !> Leading dimension (row stride) of A.
      integer(I32), intent(in) :: LDA
      !> Input matrix A.
      real(R64), intent(in) :: A(LDA, *)
      !> Leading dimension (row stride) of B.
      integer(I32), intent(in) :: LDB
      !> Input matrix B.
      real(R64), intent(in) :: B(LDB, *)
      !> Scalar multiplier for C.
      real(R64), intent(in) :: BETA
      !> Leading dimension (row stride) of C.
      integer(I32), intent(in) :: LDC
      !> In-place output matrix C on entry, alpha*op(A)*op(B) + beta*C on exit.
      real(R64), intent(inout) :: C(LDC, *)
    end subroutine

    !> BLAS ZGEMV: matrix-vector multiplication for complex numbers.
    subroutine ZGEMV(TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)
      import :: I32, R64
      !> Operation flag for A: 'N' (y=alpha*A*x+beta*y), 'T' or 'C'.
      character(len=1), intent(in) :: TRANS
      !> Matrix dimensions and leading dimension of A, and increments for X/Y.
      integer(I32), intent(in)     :: M, N, LDA, INCX, INCY
      !> Scalar multiplier for A*x.
      complex(R64), intent(in)     :: ALPHA
      !> Input matrix A.
      complex(R64), intent(in)     :: A(LDA, *)
      !> Input vector x.
      complex(R64), intent(in)     :: X(*)
      !> Scalar multiplier for y.
      complex(R64), intent(in)     :: BETA
      !> In-place output vector y on entry, alpha*A*x+beta*y on exit.
      complex(R64), intent(inout)  :: Y(*)
    end subroutine

    !> BLAS DGEMV: matrix-vector multiplication for real numbers.
    subroutine DGEMV(TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)
      import :: I32, R64
      !> Operation flag for A: 'N' (y=alpha*A*x+beta*y), 'T'.
      character(len=1), intent(in) :: TRANS
      !> Matrix dimensions and leading dimension of A, and increments for X/Y.
      integer(I32), intent(in)     :: M, N, LDA, INCX, INCY
      !> Scalar multiplier for A*x.
      real(R64), intent(in)        :: ALPHA
      !> Input matrix A.
      real(R64), intent(in)        :: A(LDA, *)
      !> Input vector x.
      real(R64), intent(in)        :: X(*)
      !> Scalar multiplier for y.
      real(R64), intent(in)        :: BETA
      !> In-place output vector y on entry, alpha*A*x+beta*y on exit.
      real(R64), intent(inout)     :: Y(*)
    end subroutine

  end interface

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function WrapDZNRM2(in) result(res)
    real(R64)                :: res
    !> Input complex vector to compute the Euclidean 2-norm.
    complex(R64), intent(in), contiguous  :: in(:)

    integer(I32) :: n
    integer(I32), parameter :: INCX = 1

    n = size(in, 1)
    res = DZNRM2(n, in, INCX)

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function WrapDNRM2(in) result(res)
    real(R64)                :: res
    !> Input real vector to compute the Euclidean 2-norm.
    real(R64), intent(in), contiguous  :: in(:)

    integer(I32) :: n
    integer(I32), parameter :: INCX = 1

    n = size(in, 1)
    res = DNRM2(n, in, INCX)

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine WrapZGEMM(outputMatrix, inputMatrixA, inputMatrixB, transA, transB)

    !> Output matrix C = op(A)*op(B). Must be sized (m,n) consistent with transA/transB.
    complex(R64), intent(out), contiguous             :: outputMatrix(:, :)
    !> Input matrix A.
    complex(R64), intent(in), contiguous  :: inputMatrixA(:, :)
    !> Input matrix B.
    complex(R64), intent(in), contiguous  :: inputMatrixB(:, :)
    !> Operation flag for A: 'N', 'T' or 'C'.
    character, intent(in)                 :: transA
    !> Operation flag for B: 'N', 'T' or 'C'.
    character, intent(in)                 :: transB

    integer(I32) :: mA, nA, mB, nB
    integer(I32) :: m, n, k
    complex(R64), parameter :: ALPHA = (1.0_R64, 0.0_R64)
    complex(R64), parameter :: BETA = (0.0_R64, 0.0_R64)

    ! Initialize variables
    mA = size(inputMatrixA, 1)
    nA = size(inputMatrixA, 2)
    mB = size(inputMatrixB, 1)
    nB = size(inputMatrixB, 2)

    if (transA .eq. 'N') then
      m = mA
      k = nA
    else
      m = nA
      k = mA
    end if

    if (transB .eq. 'N') then
      k = mB
      n = nB
    else
      k = nB
      n = mB
    end if

    ! check matrix size
    if (size(outputMatrix, 1) .ne. m) error stop "WrapZGEMM wrong matrix result size"
    if (size(outputMatrix, 2) .ne. n) error stop "WrapZGEMM wrong matrix result size"

    ! Step 1: Compute A^H * v (using zgemm)
    call ZGEMM(transA, transB, m, n, k, ALPHA, inputMatrixA, mA, inputMatrixB, mB, BETA, outputMatrix, m)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine WrapDGEMM(outputMatrix, inputMatrixA, inputMatrixB, transA, transB)

    !> Output matrix C = op(A)*op(B). Must be sized (m,n) consistent with transA/transB.
    real(R64), intent(out), contiguous             :: outputMatrix(:, :)
    !> Input matrix A.
    real(R64), intent(in), contiguous  :: inputMatrixA(:, :)
    !> Input matrix B.
    real(R64), intent(in), contiguous  :: inputMatrixB(:, :)
    !> Operation flag for A: 'N' or 'T'.
    character, intent(in)                 :: transA
    !> Operation flag for B: 'N' or 'T'.
    character, intent(in)                 :: transB

    integer(I32) :: mA, nA, mB, nB
    integer(I32) :: m, n, k
    real(R64), parameter :: ALPHA = 1.0_R64
    real(R64), parameter :: BETA = 0.0_R64

    ! Initialize variables
    mA = size(inputMatrixA, 1)
    nA = size(inputMatrixA, 2)
    mB = size(inputMatrixB, 1)
    nB = size(inputMatrixB, 2)

    if (transA .eq. 'N') then
      m = mA
      k = nA
    else
      m = nA
      k = mA
    end if

    if (transB .eq. 'N') then
      k = mB
      n = nB
    else
      k = nB
      n = mB
    end if

    ! check matrix size
    if (size(outputMatrix, 1) .ne. m) error stop "WrapDGEMM wrong matrix result size"
    if (size(outputMatrix, 2) .ne. n) error stop "WrapDGEMM wrong matrix result size"

    ! Step 1: Compute A^H * v (using zgemm)
    call DGEMM(transA, transB, m, n, k, ALPHA, inputMatrixA, mA, inputMatrixB, mB, BETA, outputMatrix, m)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine WrapAllocatableZGEMM(outputMatrix, inputMatrixA, inputMatrixB, transA, transB, allocatableQ)
    use M_Utils_UnusedVariables

    !> Output matrix C = op(A)*op(B). Allocated to size (m,n).
    complex(R64), intent(out), allocatable :: outputMatrix(:, :)
    !> Input matrix A.
    complex(R64), intent(in), contiguous  :: inputMatrixA(:, :)
    !> Input matrix B.
    complex(R64), intent(in), contiguous  :: inputMatrixB(:, :)
    !> Operation flag for A: 'N', 'T' or 'C'.
    character, intent(in)                 :: transA
    !> Operation flag for B: 'N', 'T' or 'C'.
    character, intent(in)                 :: transB
    !> Ignored hint to indicate allocatable output (kept for API symmetry).
    logical, intent(in)                   :: allocatableQ

    integer(I32) :: mA, nA, mB, nB
    integer(I32) :: m, n, k
    complex(R64), parameter :: ALPHA = (1.0_R64, 0.0_R64)
    complex(R64), parameter :: BETA = (0.0_R64, 0.0_R64)

    if (.false.) call UnusedVariables_Mark(allocatableQ)

    ! Initialize variables
    mA = size(inputMatrixA, 1)
    nA = size(inputMatrixA, 2)
    mB = size(inputMatrixB, 1)
    nB = size(inputMatrixB, 2)

    if (transA .eq. 'N') then
      m = mA
      k = nA
    else
      m = nA
      k = mA
    end if

    if (transB .eq. 'N') then
      k = mB
      n = nB
    else
      k = nB
      n = mB
    end if

    ! Allocate temporary matrix
    if (.not. allocated(outputMatrix)) allocate (outputMatrix(m, n))

    ! Step 1: Compute A^H * v (using zgemm)
    call ZGEMM(transA, transB, m, n, k, ALPHA, inputMatrixA, mA, inputMatrixB, mB, BETA, outputMatrix, m)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine WrapAllocatableDGEMM(outputMatrix, inputMatrixA, inputMatrixB, transA, transB, allocatableQ)
    use M_Utils_UnusedVariables

    !> Output matrix C = op(A)*op(B). Allocated to size (m,n).
    real(R64), intent(out), allocatable :: outputMatrix(:, :)
    !> Input matrix A.
    real(R64), intent(in), contiguous  :: inputMatrixA(:, :)
    !> Input matrix B.
    real(R64), intent(in), contiguous  :: inputMatrixB(:, :)
    !> Operation flag for A: 'N' or 'T'.
    character, intent(in)                 :: transA
    !> Operation flag for B: 'N' or 'T'.
    character, intent(in)                 :: transB
    !> Ignored hint to indicate allocatable output (kept for API symmetry).
    logical, intent(in)                   :: allocatableQ

    integer(I32) :: mA, nA, mB, nB
    integer(I32) :: m, n, k
    real(R64), parameter :: ALPHA = (1.0_R64, 0.0_R64)
    real(R64), parameter :: BETA = (0.0_R64, 0.0_R64)

    if (.false.) call UnusedVariables_Mark(allocatableQ)

    ! Initialize variables
    mA = size(inputMatrixA, 1)
    nA = size(inputMatrixA, 2)
    mB = size(inputMatrixB, 1)
    nB = size(inputMatrixB, 2)

    if (transA .eq. 'N') then
      m = mA
      k = nA
    else
      m = nA
      k = mA
    end if

    if (transB .eq. 'N') then
      k = mB
      n = nB
    else
      k = nB
      n = mB
    end if

    ! Allocate temporary matrix
    if (.not. allocated(outputMatrix)) allocate (outputMatrix(m, n))

    ! Step 1: Compute A^H * v (using zgemm)
    call DGEMM(transA, transB, m, n, k, ALPHA, inputMatrixA, mA, inputMatrixB, mB, BETA, outputMatrix, m)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine WrapZGEMV(outputVector, inputMatrixA, inputVectorX, transA)

    !> Output vector y = op(A)*x. Must be sized m consistent with transA.
    complex(R64), intent(out), contiguous :: outputVector(:)
    !> Input matrix A.
    complex(R64), intent(in), contiguous  :: inputMatrixA(:, :)
    !> Input vector x.
    complex(R64), intent(in), contiguous  :: inputVectorX(:)
    !> Operation flag for A: 'N', 'T' or 'C'.
    character, intent(in)                 :: transA

    integer(I32) :: mA, nA
    integer(I32) :: m, n
    integer(I32), parameter :: INCX = 1, INCY = 1
    complex(R64), parameter :: ALPHA = (1.0_R64, 0.0_R64)
    complex(R64), parameter :: BETA = (0.0_R64, 0.0_R64)

    mA = size(inputMatrixA, 1)
    nA = size(inputMatrixA, 2)

    if (transA .eq. 'N') then
      m = mA
      n = nA
    else
      m = nA
      n = mA
    end if

    if (size(inputVectorX, 1) .ne. n) error stop "WrapZGEMV wrong input vector size"
    if (size(outputVector, 1) .ne. m) error stop "WrapZGEMV wrong result vector size"

    call ZGEMV(transA, m, n, ALPHA, inputMatrixA, mA, inputVectorX, INCX, BETA, outputVector, INCY)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine WrapDGEMV(outputVector, inputMatrixA, inputVectorX, transA)

    !> Output vector y = op(A)*x. Must be sized m consistent with transA.
    real(R64), intent(out), contiguous :: outputVector(:)
    !> Input matrix A.
    real(R64), intent(in), contiguous  :: inputMatrixA(:, :)
    !> Input vector x.
    real(R64), intent(in), contiguous  :: inputVectorX(:)
    !> Operation flag for A: 'N' or 'T'.
    character, intent(in)              :: transA

    integer(I32) :: mA, nA
    integer(I32) :: m, n
    integer(I32), parameter :: INCX = 1, INCY = 1
    real(R64), parameter :: ALPHA = 1.0_R64
    real(R64), parameter :: BETA = 0.0_R64

    mA = size(inputMatrixA, 1)
    nA = size(inputMatrixA, 2)

    if (transA .eq. 'N') then
      m = mA
      n = nA
    else
      m = nA
      n = mA
    end if

    if (size(inputVectorX, 1) .ne. n) error stop "WrapDGEMV wrong input vector size"
    if (size(outputVector, 1) .ne. m) error stop "WrapDGEMV wrong result vector size"

    call DGEMV(transA, m, n, ALPHA, inputMatrixA, mA, inputVectorX, INCX, BETA, outputVector, INCY)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine WrapAllocatableZGEMV(outputVector, inputMatrixA, inputVectorX, transA, allocatableQ)
    use M_Utils_UnusedVariables

    !> Output vector y = op(A)*x. Allocated to size m.
    complex(R64), intent(out), allocatable :: outputVector(:)
    !> Input matrix A.
    complex(R64), intent(in), contiguous   :: inputMatrixA(:, :)
    !> Input vector x.
    complex(R64), intent(in), contiguous   :: inputVectorX(:)
    !> Operation flag for A: 'N', 'T' or 'C'.
    character, intent(in)                  :: transA
    !> Ignored hint to indicate allocatable output (kept for API symmetry).
    logical, intent(in)                    :: allocatableQ

    integer(I32) :: mA, nA
    integer(I32) :: m, n
    integer(I32), parameter :: INCX = 1, INCY = 1
    complex(R64), parameter :: ALPHA = (1.0_R64, 0.0_R64)
    complex(R64), parameter :: BETA = (0.0_R64, 0.0_R64)

    if (.false.) call UnusedVariables_Mark(allocatableQ)

    mA = size(inputMatrixA, 1)
    nA = size(inputMatrixA, 2)

    if (transA .eq. 'N') then
      m = mA
      n = nA
    else
      m = nA
      n = mA
    end if

    if (size(inputVectorX, 1) .ne. n) error stop "WrapAllocatableZGEMV wrong input vector size"

    if (.not. allocated(outputVector)) allocate (outputVector(m))

    call ZGEMV(transA, m, n, ALPHA, inputMatrixA, mA, inputVectorX, INCX, BETA, outputVector, INCY)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine WrapAllocatableDGEMV(outputVector, inputMatrixA, inputVectorX, transA, allocatableQ)
    use M_Utils_UnusedVariables

    !> Output vector y = op(A)*x. Allocated to size m.
    real(R64), intent(out), allocatable :: outputVector(:)
    !> Input matrix A.
    real(R64), intent(in), contiguous   :: inputMatrixA(:, :)
    !> Input vector x.
    real(R64), intent(in), contiguous   :: inputVectorX(:)
    !> Operation flag for A: 'N' or 'T'.
    character, intent(in)               :: transA
    !> Ignored hint to indicate allocatable output (kept for API symmetry).
    logical, intent(in)                 :: allocatableQ

    integer(I32) :: mA, nA
    integer(I32) :: m, n
    integer(I32), parameter :: INCX = 1, INCY = 1
    real(R64), parameter :: ALPHA = 1.0_R64
    real(R64), parameter :: BETA = 0.0_R64

    if (.false.) call UnusedVariables_Mark(allocatableQ)

    mA = size(inputMatrixA, 1)
    nA = size(inputMatrixA, 2)

    if (transA .eq. 'N') then
      m = mA
      n = nA
    else
      m = nA
      n = mA
    end if

    if (size(inputVectorX, 1) .ne. n) error stop "WrapAllocatableDGEMV wrong input vector size"

    if (.not. allocated(outputVector)) allocate (outputVector(m))

    call DGEMV(transA, m, n, ALPHA, inputMatrixA, mA, inputVectorX, INCX, BETA, outputVector, INCY)

  end subroutine

end module
