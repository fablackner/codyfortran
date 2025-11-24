! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> LAPACK helper interfaces for factorizations, solves, and eigenproblems.
!>
!> Collects typed Fortran interfaces and thin wrappers (real/complex, double
!> precision) around selected LAPACK routines used by the project (e.g.,
!> linear solves, symmetric/hermitian eigensolvers, and utility operations).
module M_Utils_LapackLib
  use M_Utils_Types

  implicit none

  ! Interface for DGETRF (LU factorization)
  interface
    subroutine DGETRF(m, n, A, lda, ipiv, info)
      import :: I32, R64
      !> Number of rows of A; number of columns of A; leading dimension of A.
      integer(I32), intent(in) :: m, n, lda
      !> Matrix to factorize; overwritten with L and U factors.
      real(R64), intent(inout) :: A(lda, *)
      !> Pivot indices.
      integer(I32), intent(out) :: ipiv(*)
      !> Status info (0 on success; >0 singular).
      integer(I32), intent(out) :: info
    end subroutine
  end interface

  ! Interface for DGETRS (Solve using LU factorization)
  interface
    subroutine DGETRS(trans, n, nrhs, A, lda, ipiv, B, ldb, info)
      import :: I32, R64
      !> Operation, 'N' no transpose, 'T' transpose, 'C' conjugate transpose.
      character(1), intent(in) :: trans
      !> Order of A; number of RHS; leading dims of A and B.
      integer(I32), intent(in) :: n, nrhs, lda, ldb
      !> LU factors from DGETRF.
      real(R64), intent(in) :: A(lda, *)
      !> Pivot indices from DGETRF.
      integer(I32), intent(in) :: ipiv(*)
      !> On entry RHS B; on exit solution X.
      real(R64), intent(inout) :: B(ldb, *)
      !> Status info (0 on success).
      integer(I32), intent(out) :: info
    end subroutine
  end interface

  ! Interface for ZGETRF (Complex LU factorization)
  interface
    subroutine ZGETRF(m, n, A, lda, ipiv, info)
      import :: I32, R64
      !> Number of rows of A; number of columns of A; leading dimension of A.
      integer(I32), intent(in) :: m, n, lda
      !> Matrix to factorize; overwritten with L and U factors.
      complex(R64), intent(inout) :: A(lda, *)
      !> Pivot indices.
      integer(I32), intent(out) :: ipiv(*)
      !> Status info (0 on success; >0 singular).
      integer(I32), intent(out) :: info
    end subroutine
  end interface

  ! Interface for ZGETRS (Solve complex system using LU factorization)
  interface
    subroutine ZGETRS(trans, n, nrhs, A, lda, ipiv, B, ldb, info)
      import :: I32, R64
      !> Operation, 'N' no transpose, 'T' transpose, 'C' conjugate transpose.
      character(1), intent(in) :: trans
      !> Order of A; number of RHS; leading dims of A and B.
      integer(I32), intent(in) :: n, nrhs, lda, ldb
      !> LU factors from ZGETRF.
      complex(R64), intent(in) :: A(lda, *)
      !> Pivot indices from ZGETRF.
      integer(I32), intent(in) :: ipiv(*)
      !> On entry RHS B; on exit solution X.
      complex(R64), intent(inout) :: B(ldb, *)
      !> Status info (0 on success).
      integer(I32), intent(out) :: info
    end subroutine
  end interface

  ! Interface for ZPOTRF (Complex Cholesky factorization)
  interface
    subroutine ZPOTRF(uplo, n, A, lda, info)
      import :: I32, R64
      !> 'U' upper triangle, 'L' lower triangle.
      character(1), intent(in) :: uplo
      !> Order of A; leading dimension of A.
      integer(I32), intent(in) :: n, lda
      !> Matrix to factorize; overwritten with Cholesky factor.
      complex(R64), intent(inout) :: A(lda, *)
      !> Status info (0 on success; >0 not positive definite).
      integer(I32), intent(out) :: info
    end subroutine
  end interface

  ! New: Interface for ZPSTRF (Hermitian PSD pivoted Cholesky)
  interface
    subroutine ZPSTRF(uplo, n, A, lda, piv, rank, tol, work, info)
      import :: I32, R64
      character(len=1), intent(in) :: uplo
      integer(I32), intent(in)     :: n, lda
      complex(R64), intent(inout)  :: A(lda, *)
      integer(I32), intent(out)    :: piv(*)
      integer(I32), intent(out)    :: rank
      real(R64), intent(in)        :: tol
      real(R64), intent(out)       :: work(*)
      integer(I32), intent(out)    :: info
    end subroutine
  end interface

  interface LapackLib_DiagonalizeSym
    procedure WrapZHEEVR
    procedure WrapDSYEVR
  end interface

  interface LapackLib_DiagonalizeGeneric
    procedure WrapZHEEV
    procedure WrapDSYEV
  end interface

  interface LapackLib_FactorizeQR
    procedure WrapZGEQRF
    procedure WrapDGEQRF
  end interface

  interface LapackLib_SolveLinearSystem
    procedure WrapZGESV
    procedure WrapDGESV
  end interface

  interface LapackLib_FactorizeLU
    procedure WrapDGETRF
    procedure WrapZGETRF
  end interface

  interface LapackLib_FactorizeCholesky
    procedure WrapZPOTRF
    procedure WrapZPSTRF
  end interface

  interface LapackLib_SolveFactorized
    procedure WrapDGETRS
    procedure WrapDGETRSsingle
    procedure WrapZGETRS
    procedure WrapZGETRSsingle
  end interface

  private WrapZHEEVR, WrapDSYEVR, WrapZHEEV, WrapDSYEV, WrapZGEQRF, WrapDGEQRF
  private WrapZGESV, WrapDGESV, WrapDGETRF, WrapDGETRS, WrapZGETRF, WrapZGETRS
  private WrapDGETRSsingle, WrapZGETRSsingle, WrapZPOTRF, WrapZPSTRF

  interface
    subroutine ZHEEVR(JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, &
                      LWORK, &
                      RWORK, LRWORK, &
                      IWORK, LIWORK, INFO)
      import :: I32, R64
      !> 'N' eigenvalues only, 'V' eigenvalues and eigenvectors.
      character(len=1), intent(in) :: JOBZ
      !> 'A' all, 'V' by value range [VL,VU], 'I' by index [IL,IU].
      character(len=1), intent(in) :: RANGE
      !> Triangle of A to reference: 'U' upper, 'L' lower.
      character(len=1), intent(in) :: UPLO
      !> Order of A; leading dimension of A.
      integer(I32), intent(in)     :: N
      integer(I32), intent(in)     :: LDA
      !> Hermitian matrix A on entry; contents destroyed on exit.
      complex(R64), intent(inout)  :: A(LDA, *)
      !> Lower/upper bounds if RANGE='V'.
      real(R64), intent(in)        :: VL
      real(R64), intent(in)        :: VU
      !> Lower/upper indices if RANGE='I'.
      integer(I32), intent(in)     :: IL
      integer(I32), intent(in)     :: IU
      !> Absolute error tolerance (<=0 uses default).
      real(R64), intent(in)        :: ABSTOL
      !> Number of eigenvalues found; eigenvalues; eigenvectors leading dim.
      integer(I32), intent(out)    :: M
      real(R64), intent(out)       :: W(*)
      integer(I32), intent(in)     :: LDZ
      !> Eigenvectors Z if JOBZ='V'.
      complex(R64), intent(out)    :: Z(LDZ, *)
      !> Support of eigenvectors.
      integer(I32), intent(out)    :: ISUPPZ(*)
      !> Workspace; length; real workspace; its length.
      complex(R64), intent(out)    :: WORK(*)
      integer(I32), intent(in)     :: LWORK
      real(R64), intent(out)       :: RWORK(*)
      integer(I32), intent(in)     :: LRWORK
      !> Integer workspace; its length; status info.
      integer(I32), intent(out)    :: IWORK(*)
      integer(I32), intent(in)     :: LIWORK
      integer(I32), intent(out)    :: INFO
    end subroutine

    subroutine DSYEVR(JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, &
                      LWORK, &
                      IWORK, LIWORK, INFO)
      import :: I32, R64
      !> 'N' eigenvalues only, 'V' eigenvalues and eigenvectors.
      character(len=1), intent(in) :: JOBZ
      !> 'A' all, 'V' by value range [VL,VU], 'I' by index [IL,IU].
      character(len=1), intent(in) :: RANGE
      !> Triangle of A to reference: 'U' upper, 'L' lower.
      character(len=1), intent(in) :: UPLO
      !> Order of A; leading dimension of A.
      integer(I32), intent(in)     :: N
      integer(I32), intent(in)     :: LDA
      !> Symmetric matrix A on entry; contents destroyed on exit.
      real(R64), intent(inout)  :: A(LDA, *)
      !> Lower/upper bounds if RANGE='V'.
      real(R64), intent(in)        :: VL
      real(R64), intent(in)        :: VU
      !> Lower/upper indices if RANGE='I'.
      integer(I32), intent(in)     :: IL
      integer(I32), intent(in)     :: IU
      !> Absolute error tolerance (<=0 uses default).
      real(R64), intent(in)        :: ABSTOL
      !> Number of eigenvalues found; eigenvalues; eigenvectors leading dim.
      integer(I32), intent(out)    :: M
      real(R64), intent(out)       :: W(*)
      integer(I32), intent(in)     :: LDZ
      !> Eigenvectors Z if JOBZ='V'.
      real(R64), intent(out)    :: Z(LDZ, *)
      !> Support of eigenvectors.
      integer(I32), intent(out)    :: ISUPPZ(*)
      !> Workspace; length; integer workspace; its length; status info.
      real(R64), intent(out)    :: WORK(*)
      integer(I32), intent(in)     :: LWORK
      integer(I32), intent(out)    :: IWORK(*)
      integer(I32), intent(in)     :: LIWORK
      integer(I32), intent(out)    :: INFO
    end subroutine

  end interface

  interface
    subroutine ZHEEV(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, &
                     RWORK, &
                     INFO)
      import :: I32, R64
      !> 'N' eigenvalues only, 'V' eigenvalues and eigenvectors.
      character(len=1), intent(in) :: JOBZ
      !> Triangle of A to reference: 'U' upper, 'L' lower.
      character(len=1), intent(in) :: UPLO
      !> Order of A; leading dimension of A.
      integer(I32), intent(in)     :: N
      integer(I32), intent(in)     :: LDA
      !> Hermitian matrix A on entry; eigenvectors if JOBZ='V'.
      complex(R64), intent(inout)  :: A(LDA, *)
      !> Eigenvalues output.
      real(R64), intent(out)       :: W(*)
      !> Workspace and its length; real workspace; info status.
      complex(R64), intent(out)    :: WORK(*)
      integer(I32), intent(in)     :: LWORK
      real(R64), intent(out)       :: RWORK(*)
      integer(I32), intent(inout)  :: INFO
    end subroutine

    subroutine DSYEV(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, &
                     INFO)
      import :: I32, R64
      !> 'N' eigenvalues only, 'V' eigenvalues and eigenvectors.
      character(len=1), intent(in) :: JOBZ
      !> Triangle of A to reference: 'U' upper, 'L' lower.
      character(len=1), intent(in) :: UPLO
      !> Order of A; leading dimension of A.
      integer(I32), intent(in)     :: N
      integer(I32), intent(in)     :: LDA
      !> Symmetric matrix A on entry; eigenvectors if JOBZ='V'.
      real(R64), intent(inout)  :: A(LDA, *)
      !> Eigenvalues output.
      real(R64), intent(out)       :: W(*)
      !> Workspace and its length; info status.
      real(R64), intent(out)    :: WORK(*)
      integer(I32), intent(in)     :: LWORK
      integer(I32), intent(inout)  :: INFO
    end subroutine

  end interface

  interface
    subroutine ZGEQRF(M, N, A, LDA, TAU, WORK, LWORK, INFO)
      import :: I32, R64
      !> Rows of A; columns of A; leading dimension of A.
      integer(I32), intent(in)     :: M
      integer(I32), intent(in)     :: N
      integer(I32), intent(in)     :: LDA
      !> Matrix A; overwritten with R in the upper triangle and reflectors below.
      complex(R64), intent(inout)  :: A(LDA, *)
      !> Scalar factors of the elementary reflectors.
      complex(R64), intent(out)    :: TAU(*)
      !> Workspace and its length; status info.
      complex(R64), intent(out)    :: WORK(*)
      integer(I32), intent(in)     :: LWORK
      integer(I32), intent(out)    :: INFO
    end subroutine

    subroutine DGEQRF(M, N, A, LDA, TAU, WORK, LWORK, INFO)
      import :: I32, R64
      !> Rows of A; columns of A; leading dimension of A.
      integer(I32), intent(in)     :: M
      integer(I32), intent(in)     :: N
      integer(I32), intent(in)     :: LDA
      !> Matrix A; overwritten with R in the upper triangle and reflectors below.
      real(R64), intent(inout)  :: A(LDA, *)
      !> Scalar factors of the elementary reflectors.
      real(R64), intent(out)    :: TAU(*)
      !> Workspace and its length; status info.
      real(R64), intent(out)    :: WORK(*)
      integer(I32), intent(in)     :: LWORK
      integer(I32), intent(out)    :: INFO
    end subroutine

  end interface

  interface
    subroutine ZUNGQR(M, N, K, A, LDA, TAU, WORK, LWORK, INFO)
      import :: I32, R64
      !> Rows and columns of Q; number of reflectors K; leading dimension of A.
      integer(I32), intent(in)     :: M
      integer(I32), intent(in)     :: N
      integer(I32), intent(in)     :: K
      integer(I32), intent(in)     :: LDA
      !> On entry: R and reflectors from GEQRF; on exit: Q.
      complex(R64), intent(inout)  :: A(LDA, *)
      !> Scalar factors of reflectors.
      complex(R64), intent(in)     :: TAU(*)
      !> Workspace and its length; status info.
      complex(R64), intent(out)    :: WORK(*)
      integer(I32), intent(in)     :: LWORK
      integer(I32), intent(out)    :: INFO
    end subroutine

    subroutine DORGQR(M, N, K, A, LDA, TAU, WORK, LWORK, INFO)
      import :: I32, R64
      !> Rows and columns of Q; number of reflectors K; leading dimension of A.
      integer(I32), intent(in)     :: M
      integer(I32), intent(in)     :: N
      integer(I32), intent(in)     :: K
      integer(I32), intent(in)     :: LDA
      !> On entry: R and reflectors from GEQRF; on exit: Q.
      real(R64), intent(inout)  :: A(LDA, *)
      !> Scalar factors of reflectors.
      real(R64), intent(in)     :: TAU(*)
      !> Workspace and its length; status info.
      real(R64), intent(out)    :: WORK(*)
      integer(I32), intent(in)     :: LWORK
      integer(I32), intent(out)    :: INFO
    end subroutine

  end interface

  interface
    subroutine ZGESV(N, NRHS, A, LDA, IPIV, B, LDB, INFO)
      import :: I32, R64
      !> Order of A; number of RHS; leading dimensions of A and B.
      integer(I32), intent(in)     :: N
      integer(I32), intent(in)     :: NRHS
      integer(I32), intent(in)     :: LDA
      integer(I32), intent(in)     :: LDB
      !> Matrix A (overwritten by LU factorization); pivot indices; RHS/solution.
      complex(R64), intent(inout)  :: A(LDA, *)
      integer(I32), intent(out)    :: IPIV(*)
      complex(R64), intent(inout)  :: B(LDB, *)
      !> Status info (0 on success).
      integer(I32), intent(out)    :: INFO
    end subroutine

    subroutine DGESV(N, NRHS, A, LDA, IPIV, B, LDB, INFO)
      import :: I32, R64
      !> Order of A; number of RHS; leading dimensions of A and B.
      integer(I32), intent(in)     :: N
      integer(I32), intent(in)     :: NRHS
      integer(I32), intent(in)     :: LDA
      integer(I32), intent(in)     :: LDB
      !> Matrix A (overwritten by LU factorization); pivot indices; RHS/solution.
      real(R64), intent(inout)     :: A(LDA, *)
      integer(I32), intent(out)    :: IPIV(*)
      real(R64), intent(inout)     :: B(LDB, *)
      !> Status info (0 on success).
      integer(I32), intent(out)    :: INFO
    end subroutine
  end interface

  interface
    subroutine DLASRT(ID, N, D, INFO)
      import :: I32, R64
      !> 'I' increasing order, 'D' decreasing order.
      character(len=1), intent(in) :: ID
      !> Length of array D.
      integer(I32), intent(in)     :: N
      !> Array to be sorted in-place.
      real(R64), intent(inout)  :: D(*)
      !> Status info (0 on success).
      integer(I32), intent(out)    :: INFO
    end subroutine
  end interface

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine WrapZHEEVR(evals, evecs, nFound, A, evecsQ, imin_, imax_, min_, max_)
    !> Output eigenvalues; allocated to nFound on exit.
    real(R64), intent(out), allocatable              :: evals(:)
    !> Output eigenvectors, size (N, nFound) if evecsQ.
    complex(R64), intent(out), allocatable           :: evecs(:, :)
    !> Number of converged eigenpairs on exit.
    integer(I32), intent(out)                        :: nFound
    !> Input Hermitian matrix (destroyed by LAPACK).
    complex(R64), intent(inout), contiguous                       :: A(:, :)!destroyed
    !> Whether to return eigenvectors.
    logical, intent(in)                              :: evecsQ
    !> Optional index range [imin_, imax_] if selected by index.
    integer(I32), intent(in), optional               :: imin_
    integer(I32), intent(in), optional               :: imax_
    !> Optional value range [min_, max_] if selected by value.
    real(R64), intent(in), optional                  :: min_
    real(R64), intent(in), optional                  :: max_

    !local
    character                 :: range
    complex(R64), allocatable :: WORK(:)
    integer(I32), allocatable :: iWORK(:)
    complex(R64)              :: OPT(1), evecsDummy(1, 1)
    integer(I32)              :: iOPT(1)
    real(R64), parameter      :: ABSTOL = 1e-15_R64
    integer(I32)              :: N, lwork, ilwork, info
    integer(I32), allocatable :: ISUPPZ(:)
    real(R64), allocatable    :: rWORK(:)
    integer(I32)              :: rlwork
    real(R64)                 :: rOPT(1)

    N = size(A, 1)

    if (.not. allocated(evals)) allocate (evals(N))

    allocate (ISUPPZ(2 * N))

    if (present(imin_)) then
      if (.not. present(imax_)) error stop "invalid arguments for WrapZHEEVR"
      range = 'I'
    else if (present(min_)) then
      if (.not. present(max_)) error stop "invalid arguments for WrapZHEEVR"
      range = 'V'
    else
      range = 'A'
    end if

    if (evecsQ) then
      if (.not. allocated(evecs)) allocate (evecs(N, N))
      call ZHEEVR('V', range, 'U', N, A, N, min_, max_, imin_, imax_, ABSTOL, nFound, evals, evecs, N, &
                  ISUPPZ, OPT, -1, rOPT, -1, iOPT, -1, info)
    else
      call ZHEEVR('N', range, 'U', N, A, N, min_, max_, imin_, imax_, ABSTOL, nFound, evals, evecsDummy, 1, &
                  ISUPPZ, OPT, -1, rOPT, -1, iOPT, -1, info)
    end if

    lwork = int(OPT(1))
    rlwork = int(rOPT(1))
    ilwork = int(iOPT(1))
    allocate (WORK(lwork))
    allocate (rWORK(rlwork))
    allocate (iWORK(ilwork))

    if (evecsQ) then
      call ZHEEVR('V', range, 'U', N, A, N, min_, max_, imin_, imax_, ABSTOL, nFound, evals, evecs, N, &
                  ISUPPZ, WORK, lwork, rWORK, rlwork, iWORK, ilwork, info)
    else
      call ZHEEVR('N', range, 'U', N, A, N, min_, max_, imin_, imax_, ABSTOL, nFound, evals, evecsDummy, 1, &
                  ISUPPZ, WORK, lwork, rWORK, rlwork, iWORK, ilwork, info)
    end if

    if (info .ne. 0) then
      write (6, "('diagonalizer_lapack_ZHEEVR: info = ', i20)") info
      error stop 'error in Wraputil_ZHEEVR'
    end if

    evals = evals(1:nFound)
    evecs = evecs(1:N, 1:nFound)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine WrapDSYEVR(evals, evecs, nFound, A, evecsQ, imin_, imax_, min_, max_)
    !> Output eigenvalues; allocated to nFound on exit.
    real(R64), intent(out), allocatable              :: evals(:)
    !> Output eigenvectors, size (N, nFound) if evecsQ.
    real(R64), intent(out), allocatable           :: evecs(:, :)
    !> Number of converged eigenpairs on exit.
    integer(I32), intent(out)                        :: nFound
    !> Input symmetric matrix (destroyed by LAPACK).
    real(R64), intent(inout), contiguous                       :: A(:, :)!destroyed
    !> Whether to return eigenvectors.
    logical, intent(in)                              :: evecsQ
    !> Optional index range [imin_, imax_] if selected by index.
    integer(I32), intent(in), optional               :: imin_, imax_
    !> Optional value range [min_, max_] if selected by value.
    real(R64), intent(in), optional                  :: min_, max_

    !local
    character                 :: range
    real(R64), allocatable :: WORK(:)
    integer(I32), allocatable :: iWORK(:)
    real(R64)              :: OPT(1), evecsDummy(1, 1)
    integer(I32)              :: iOPT(1)
    real(R64), parameter      :: ABSTOL = 1e-15_R64
    integer(I32)              :: N, lwork, ilwork, info
    integer(I32), allocatable :: ISUPPZ(:)

    N = size(A, 1)

    if (.not. allocated(evals)) allocate (evals(N))

    allocate (ISUPPZ(2 * N))

    if (present(imin_)) then
      if (.not. present(imax_)) error stop "invalid arguments for WrapDSYEVR"
      range = 'I'
    else if (present(min_)) then
      if (.not. present(max_)) error stop "invalid arguments for WrapDSYEVR"
      range = 'V'
    else
      range = 'A'
    end if

    if (evecsQ) then
      if (.not. allocated(evecs)) allocate (evecs(N, N))
      call DSYEVR('V', range, 'U', N, A, N, min_, max_, imin_, imax_, ABSTOL, nFound, evals, evecs, N, &
                  ISUPPZ, OPT, -1, iOPT, -1, info)
    else
      call DSYEVR('N', range, 'U', N, A, N, min_, max_, imin_, imax_, ABSTOL, nFound, evals, evecsDummy, 1, &
                  ISUPPZ, OPT, -1, iOPT, -1, info)
    end if

    lwork = int(OPT(1))
    ilwork = int(iOPT(1))
    allocate (WORK(lwork))
    allocate (iWORK(ilwork))

    if (evecsQ) then
      call DSYEVR('V', range, 'U', N, A, N, min_, max_, imin_, imax_, ABSTOL, nFound, evals, evecs, N, &
                  ISUPPZ, WORK, lwork, iWORK, ilwork, info)
    else
      call DSYEVR('N', range, 'U', N, A, N, min_, max_, imin_, imax_, ABSTOL, nFound, evals, evecsDummy, 1, &
                  ISUPPZ, WORK, lwork, iWORK, ilwork, info)
    end if

    if (info .ne. 0) then
      write (6, "('diagonalizer_lapack_DSYEVR: info = ', i20)") info
      error stop 'error in Wraputil_DSYEVR'
    end if

    evals = evals(1:nFound)
    evecs = evecs(1:N, 1:nFound)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine WrapZHEEV(evals, evecs, A, evecsQ)
    !> Output eigenvalues; allocated to size N.
    real(R64), intent(out), allocatable              :: evals(:)
    !> Output eigenvectors if evecsQ.
    complex(R64), intent(out), allocatable           :: evecs(:, :)
    !> Input Hermitian matrix; not destroyed (copied internally).
    complex(R64), intent(in), contiguous                          :: A(:, :)
    !> Whether to return eigenvectors.
    logical, intent(in)                              :: evecsQ

    !local
    character                 :: jobz
    complex(R64), allocatable :: ATmp(:, :)
    complex(R64), allocatable :: WORK(:)
    complex(R64)              :: OPT(1)
    integer(I32)              :: N, lwork, info
    real(R64), allocatable    :: rWORK(:)
    integer(I32)              :: rlwork

    N = size(A, 1)

    if (.not. allocated(evals)) allocate (evals(N))

    if (evecsQ) then
      jobz = 'V'
    else
      jobz = 'N'
    end if

    allocate (ATmp, source=A)

    rlwork = 3 * N - 2
    allocate (rWORK(rlwork))

    call ZHEEV(jobz, 'U', N, ATmp, N, evals, OPT, -1, rWORK, info)

    lwork = int(OPT(1))
    allocate (WORK(lwork))

    call ZHEEV(jobz, 'U', N, ATmp, N, evals, WORK, lwork, rWORK, info)

    if (info .ne. 0) then
      write (6, "('LapackLib_DiagonalizeGeneric: info = ', i20)") info
      error stop 'error in Wraputil_ZHEEV'
    end if

    if (evecsQ) evecs = ATmp

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine WrapDSYEV(evals, evecs, A, evecsQ)
    !> Output eigenvalues; allocated to size N.
    real(R64), intent(out), allocatable              :: evals(:)
    !> Output eigenvectors if evecsQ.
    real(R64), intent(out), allocatable           :: evecs(:, :)
    !> Input symmetric matrix; not destroyed (copied internally).
    real(R64), intent(in), contiguous                          :: A(:, :)
    !> Whether to return eigenvectors.
    logical, intent(in)                              :: evecsQ

    !local
    character                 :: jobz
    real(R64), allocatable :: ATmp(:, :)
    real(R64), allocatable :: WORK(:)
    real(R64)              :: OPT(1)
    integer(I32)              :: N, lwork, info

    N = size(A, 1)

    if (.not. allocated(evals)) allocate (evals(N))

    if (evecsQ) then
      jobz = 'V'
    else
      jobz = 'N'
    end if

    allocate (ATmp, source=A)

    call DSYEV(jobz, 'U', N, ATmp, N, evals, OPT, -1, info)

    lwork = int(OPT(1))
    allocate (WORK(lwork))

    call DSYEV(jobz, 'U', N, ATmp, N, evals, WORK, lwork, info)

    if (info .ne. 0) then
      write (6, "('LapackLib_DiagonalizeGeneric: info = ', i20)") info
      error stop 'error in Wraputil_DSYEV'
    end if

    if (evecsQ) evecs = ATmp

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine WrapZGEQRF(tri, A)
    !> Output upper-triangular R factor.
    complex(R64), intent(out), allocatable :: tri(:, :)
    !> Input matrix A; overwritten with Q on exit.
    complex(R64), intent(inout), contiguous             :: A(:, :)

    !local
    complex(R64), allocatable :: TAU(:)
    complex(R64), allocatable :: WORK(:)
    complex(R64) :: OPT(1)
    integer(I32) :: N, M, LWORK, INFO, LDA

    M = size(A, 1)
    N = size(A, 2)

    allocate (TAU(N))

    if (M < N) error stop 'ZGEQRF does not work for M<N'

    LDA = M

    call ZGEQRF(M, N, A, LDA, TAU, OPT, -1, INFO)

    LWORK = int(OPT(1))
    allocate (WORK(lwork))

    call ZGEQRF(M, N, A, LDA, TAU, WORK, LWORK, INFO)

    if (INFO .ne. 0) error stop 'fail in ZGEQRF'

    tri = A

    call ZUNGQR(M, N, N, A, M, TAU, OPT, -1, INFO)

    LWORK = int(OPT(1))
    deallocate (WORK)
    allocate (WORK(lwork))

    call ZUNGQR(M, N, N, A, M, TAU, WORK, LWORK, INFO)

    if (INFO .ne. 0) error stop 'fail in ZUNGQR'

    ! on output the A contains the orthogonal matrix Q and the
    ! tri matrix contains the triangular matrix R of the QR factorization
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine WrapDGEQRF(tri, A)
    !> Output upper-triangular R factor.
    real(R64), intent(out), allocatable :: tri(:, :)
    !> Input matrix A; overwritten with Q on exit.
    real(R64), intent(inout), contiguous             :: A(:, :)

    !local
    real(R64), allocatable :: TAU(:)
    real(R64), allocatable :: WORK(:)
    real(R64) :: OPT(1)
    integer(I32) :: N, M, LWORK, INFO, LDA

    M = size(A, 1)
    N = size(A, 2)

    allocate (TAU(N))

    if (M < N) error stop 'DGEQRF does not work for M<N'

    LDA = M

    call DGEQRF(M, N, A, LDA, TAU, OPT, -1, INFO)

    LWORK = int(OPT(1))
    allocate (WORK(lwork))

    call DGEQRF(M, N, A, LDA, TAU, WORK, LWORK, INFO)

    if (INFO .ne. 0) error stop 'fail in DGEQRF'

    tri = A

    call DORGQR(M, N, N, A, M, TAU, OPT, -1, INFO)

    LWORK = int(OPT(1))
    deallocate (WORK)
    allocate (WORK(lwork))

    call DORGQR(M, N, N, A, M, TAU, WORK, LWORK, INFO)

    if (INFO .ne. 0) error stop 'fail in DORGQR'

    ! on output the A contains the orthogonal matrix Q and the
    ! tri matrix contains the triangular matrix R of the QR factorization
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine WrapZGESV(A, B)
    !> Matrix A (overwritten with LU factors).
    complex(R64), intent(inout), contiguous :: A(:, :)
    !> Right-hand sides, overwritten with solution.
    complex(R64), intent(inout), contiguous :: B(:, :)

    integer(I32) :: N, NRHS, LDA, LDB, INFO
    integer(I32), allocatable :: IPIV(:)

    N = size(A, 1)
    NRHS = size(B, 2)
    LDA = N
    LDB = size(B, 1)

    allocate (IPIV(N))

    call ZGESV(N, NRHS, A, LDA, IPIV, B, LDB, INFO)

    if (INFO .ne. 0) then
      if (INFO < 0) then
        error stop 'ZGESV: Illegal parameter value'
      else
        error stop 'ZGESV: U(i,i) is exactly zero'
      end if
    end if
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine WrapDGESV(A, B)
    !> Matrix A (overwritten with LU factors).
    real(R64), intent(inout), contiguous :: A(:, :)
    !> Right-hand sides, overwritten with solution.
    real(R64), intent(inout), contiguous :: B(:, :)

    integer(I32) :: N, NRHS, LDA, LDB, INFO
    integer(I32), allocatable :: IPIV(:)

    N = size(A, 1)
    NRHS = size(B, 2)
    LDA = N
    LDB = size(B, 1)

    allocate (IPIV(N))

    call DGESV(N, NRHS, A, LDA, IPIV, B, LDB, INFO)

    if (INFO .ne. 0) then
      if (INFO < 0) then
        error stop 'DGESV: Illegal parameter value'
      else
        error stop 'DGESV: U(i,i) is exactly zero'
      end if
    end if
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine LapackLib_SortIntegerArray(array)

    !> Integer array to be sorted in-place in decreasing order.
    integer(I32), intent(inout), contiguous  :: array(:)

    real(R64) :: sorted(size(array))
    integer(I32) :: info

    sorted = real(array)
    call DLASRT('D', size(array), sorted, info)
    array = int(sorted)

    if (info .ne. 0) error stop 'fail in sortintegerarray'

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine WrapDGETRF(A, ipiv)
    !> Matrix to be factorized; replaced with L and U factors.
    real(R64), intent(inout), contiguous  :: A(:, :)
    !> Pivot indices.
    integer(I32), intent(out), contiguous :: ipiv(:)

    integer(I32) :: m, n, lda, info

    m = size(A, 1)
    n = size(A, 2)
    lda = max(1, m)

    call DGETRF(m, n, A, lda, ipiv, info)

    if (info .ne. 0) then
      if (info < 0) then
        write (6, "('DGETRF: INFO = ', i0, ' → the ', i0, '-th argument had an illegal value')") info, -info
      else
        write (6, "('DGETRF: INFO = ', i0, ' → U(', i0, ',', i0, ') is exactly zero; matrix is singular')") info, info, info
      end if
      error stop 'Error in WrapDGETRF'
    end if
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine WrapDGETRS(A, ipiv, B, trans)
    !> LU factorization from DGETRF.
    real(R64), intent(in), contiguous    :: A(:, :)
    !> Pivot indices from DGETRF.
    integer(I32), intent(in), contiguous :: ipiv(:)
    !> Right-hand sides, replaced with solution.
    real(R64), intent(inout), contiguous :: B(:, :)
    !> Operation, 'N', 'T', or 'C'.
    character(len=1), intent(in)         :: trans

    integer(I32) :: n, nrhs, lda, ldb, info

    n = size(A, 1)
    nrhs = size(B, 2)
    lda = max(1, n)
    ldb = max(1, size(B, 1))

    call DGETRS(trans, n, nrhs, A, lda, ipiv, B, ldb, info)

    if (info .ne. 0) then
      if (info < 0) then
        write (6, "('DGETRS: INFO = ', i0, ' → the ', i0, '-th argument had an illegal value')") info, -info
      else
        write (6, "('DGETRS: INFO = ', i0, ' (unexpected positive value)')") info
      end if
      error stop 'Error in WrapDGETRS'
    end if
  end subroutine

  ! New: single RHS wrapper delegating to WrapDGETRS without copying
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine WrapDGETRSsingle(A, ipiv, b, trans)
    !> LU factorization from DGETRF.
    real(R64), intent(in), contiguous                 :: A(:, :)
    !> Pivot indices from DGETRF.
    integer(I32), intent(in), contiguous              :: ipiv(:)
    !> RHS vector, replaced with solution.
    real(R64), intent(inout), contiguous, target      :: b(:)
    !> Operation, 'N', 'T', or 'C'.
    character(len=1), intent(in)                      :: trans

    real(R64), contiguous, pointer                    :: B2D(:, :)

    B2D(1:size(b), 1:1) => b
    call WrapDGETRS(A, ipiv, B2D, trans)
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine WrapZGETRF(A, ipiv)
    !> Matrix to be factorized; replaced with L and U factors.
    complex(R64), intent(inout), contiguous :: A(:, :)
    !> Pivot indices.
    integer(I32), intent(out), contiguous   :: ipiv(:)

    integer(I32) :: m, n, lda, info

    m = size(A, 1)
    n = size(A, 2)
    lda = max(1, m)

    call ZGETRF(m, n, A, lda, ipiv, info)

    if (info .ne. 0) then
      if (info < 0) then
        write (6, "('ZGETRF: INFO = ', i0, ' → the ', i0, '-th argument had an illegal value')") info, -info
      else
        write (6, "('ZGETRF: INFO = ', i0, ' → U(', i0, ',', i0, ') is exactly zero; matrix is singular')") info, info, info
      end if
      error stop 'Error in WrapZGETRF'
    end if
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine WrapZGETRS(A, ipiv, B, trans)
    !> LU factorization from ZGETRF.
    complex(R64), intent(in), contiguous     :: A(:, :)
    !> Pivot indices from ZGETRF.
    integer(I32), intent(in), contiguous     :: ipiv(:)
    !> Right-hand sides, replaced with solution.
    complex(R64), intent(inout), contiguous  :: B(:, :)
    !> Operation, 'N', 'T', or 'C'.
    character(len=1), intent(in)             :: trans

    integer(I32) :: n, nrhs, lda, ldb, info

    n = size(A, 1)
    nrhs = size(B, 2)
    lda = max(1, n)
    ldb = max(1, size(B, 1))

    call ZGETRS(trans, n, nrhs, A, lda, ipiv, B, ldb, info)

    if (info .ne. 0) then
      if (info < 0) then
        write (6, "('ZGETRS: INFO = ', i0, ' → the ', i0, '-th argument had an illegal value')") info, -info
      else
        write (6, "('ZGETRS: INFO = ', i0, ' (unexpected positive value)')") info
      end if
      error stop 'Error in WrapZGETRS'
    end if
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine WrapZGETRSsingle(A, ipiv, b, trans)
    !> LU factorization from ZGETRF.
    complex(R64), intent(in), contiguous              :: A(:, :)
    !> Pivot indices from ZGETRF.
    integer(I32), intent(in), contiguous              :: ipiv(:)
    !> RHS vector, replaced with solution.
    complex(R64), intent(inout), contiguous, target   :: b(:)
    !> Operation, 'N', 'T', or 'C'.
    character(len=1), intent(in)                      :: trans

    complex(R64), contiguous, pointer                 :: B2D(:, :)

    B2D(1:size(b), 1:1) => b
    call WrapZGETRS(A, ipiv, B2D, trans)
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine WrapZPOTRF(A, uplo)
    !> Matrix to be factorized; replaced with Cholesky factor.
    complex(R64), intent(inout), contiguous :: A(:, :)
    !> 'U' for upper triangle, 'L' for lower.
    character(len=1), intent(in)             :: uplo

    integer(I32) :: n, lda, info

    n = size(A, 1)
    lda = max(1, n)

    call ZPOTRF(uplo, n, A, lda, info)

    if (info .ne. 0) then
      if (info < 0) then
        write (6, "('ZPOTRF: INFO = ', i0, ' → the ', i0, '-th argument had an illegal value')") info, -info
      else
        write (6, "('ZPOTRF: INFO = ', i0, ' → the leading minor of order ', i0, ' is not positive definite')") info, info
      end if
      error stop 'Error in WrapZPOTRF'
    end if
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine WrapZPSTRF(A, piv, rank, uplo, tol)
    !> Pivoted Cholesky (complete pivoting) for Hermitian PSD matrices.
    complex(R64), intent(inout), contiguous :: A(:, :)
    integer(I32), intent(out), contiguous   :: piv(:)
    integer(I32), intent(out)               :: rank
    character(len=1), intent(in)            :: uplo
    real(R64), intent(in), optional         :: tol

    integer(I32) :: n, lda, info
    real(R64), allocatable :: work(:)
    real(R64) :: tol_eff

    n = size(A, 1)
    if (size(A, 2) .ne. n) error stop 'WrapZPSTRF: A must be square'
    if (size(piv) < n) error stop 'WrapZPSTRF: piv size < N'

    lda = max(1, n)
    allocate (work(n))

    if (present(tol)) then
      tol_eff = tol
    else
      tol_eff = -1.0_R64
    end if

    call ZPSTRF(uplo, n, A, lda, piv, rank, tol_eff, work, info)

    if (info .ne. 0) then
      if (info < 0) then
        write (6, "('ZPSTRF: INFO = ', i0, ' → the ', i0, '-th argument had an illegal value')") info, -info
      else
        write (6, "('ZPSTRF: INFO = ', i0, ' (unexpected positive value)')") info
      end if
      error stop 'Error in WrapZPSTRF'
    end if
  end subroutine

end module
