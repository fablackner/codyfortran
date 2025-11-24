! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> ARPACK wrappers for iterative sparse eigenproblems (real/complex).
!>
!> Provides typed interfaces and helper drivers to run ARPACK’s reverse
!> communication interface for the most common use-cases in this codebase.
module M_Utils_ArpackLib
  use M_Utils_Types

  implicit none

  interface ArpackLib_Diagonalize
    procedure WrapZNAUPD
    procedure WrapDSAUPD
  end interface

  private WrapZNAUPD, WrapDSAUPD
  private DSAUPD, DSEUPD, ZNAUPD, ZNEUPD, I_ApplyMatOnVecComplex, I_ApplyMatOnVecReal

  abstract interface
    subroutine I_ApplyMatOnVecComplex(dState, state)
      import :: R64
      !> Output y = A*x corresponding to the applied operator.
      complex(R64), intent(out), contiguous, target :: dState(:)
      !> Input vector x to be multiplied by the operator.
      complex(R64), intent(in), contiguous, target  :: state(:)
    end subroutine
  end interface

  abstract interface
    subroutine I_ApplyMatOnVecReal(dState, state)
      import :: R64
      !> Output y = A*x corresponding to the applied operator.
      real(R64), intent(out), contiguous, target :: dState(:)
      !> Input vector x to be multiplied by the operator.
      real(R64), intent(in), contiguous, target  :: state(:)
    end subroutine
  end interface

  interface
    subroutine ZNAUPD(IDO, BMAT, N, WHICH, NEV, TOL, &
                      RESID, NCV, V, LDV, IPARAM, IPNTR, WORKD, WORKL, LWORKL, &
                      RWORK, &
                      INFO)
      import :: I32, R64
      !> Reverse communication flag (set by ARPACK).
      integer(I32), intent(inout)  :: IDO
      !> Type of B in Ax = λ Bx ('I' identity, 'G' general).
      character(len=1), intent(in) :: BMAT
      !> Dimension of the eigenproblem.
      integer(I32), intent(in)     :: N
      !> Which eigenvalues to compute (e.g., 'LM','SM','LR','SR','LI','SI').
      character(len=2), intent(in) :: WHICH
      !> Number of eigenvalues to compute.
      integer(I32), intent(in)     :: NEV
      !> Convergence tolerance (0 => machine precision).
      real(R64), intent(in)        :: TOL
      !> Residual / initial vector (modified by ARPACK).
      complex(R64), intent(inout)  :: RESID(*)
      !> Number of Arnoldi vectors (subspace size).
      integer(I32), intent(in)     :: NCV
      !> Leading dimension of V.
      integer(I32), intent(in)     :: LDV
      !> Arnoldi basis vectors (input/output).
      complex(R64), intent(inout)  :: V(LDV, *)
      !> Integer parameters per ARPACK documentation.
      integer(I32), intent(inout)  :: IPARAM(11)
      !> Reverse-communication pointers.
      integer(I32), intent(inout)  :: IPNTR(14)
      !> Work array used in reverse communication.
      complex(R64), intent(inout)  :: WORKD(3 * N)
      !> Length of WORKL.
      integer(I32), intent(in)     :: LWORKL
      !> Internal work array used by ARPACK.
      complex(R64), intent(inout)  :: WORKL(LWORKL)
      !> Real work array.
      real(R64), intent(out)       :: RWORK(*)
      !> Status/diagnostic info (0 on success).
      integer(I32), intent(inout)  :: INFO
    end subroutine

    subroutine DSAUPD(IDO, BMAT, N, WHICH, NEV, TOL, &
                      RESID, NCV, V, LDV, IPARAM, IPNTR, WORKD, WORKL, LWORKL, &
                      INFO)
      import :: I32, R64
      !> Reverse communication flag (set by ARPACK).
      integer(I32), intent(inout)  :: IDO
      !> Type of B in Ax = λ Bx ('I' identity, 'G' general).
      character(len=1), intent(in) :: BMAT
      !> Dimension of the eigenproblem.
      integer(I32), intent(in)     :: N
      !> Which eigenvalues to compute (e.g., 'LM','SM','LR','SR').
      character(len=2), intent(in) :: WHICH
      !> Number of eigenvalues to compute.
      integer(I32), intent(in)     :: NEV
      !> Convergence tolerance (0 => machine precision).
      real(R64), intent(in)        :: TOL
      !> Residual / initial vector (modified by ARPACK).
      real(R64), intent(inout)  :: RESID(*)
      !> Number of Arnoldi vectors (subspace size).
      integer(I32), intent(in)     :: NCV
      !> Leading dimension of V.
      integer(I32), intent(in)     :: LDV
      !> Arnoldi basis vectors (input/output).
      real(R64), intent(inout)  :: V(LDV, *)
      !> Integer parameters per ARPACK documentation.
      integer(I32), intent(inout)  :: IPARAM(11)
      !> Reverse-communication pointers.
      integer(I32), intent(inout)  :: IPNTR(14)
      !> Work array used in reverse communication.
      real(R64), intent(inout)  :: WORKD(3 * N)
      !> Length of WORKL.
      integer(I32), intent(in)     :: LWORKL
      !> Internal work array used by ARPACK.
      real(R64), intent(inout)  :: WORKL(LWORKL)
      !> Status/diagnostic info (0 on success).
      integer(I32), intent(inout)  :: INFO
    end subroutine

  end interface

  interface
    subroutine ZNEUPD(RVEC, HOWMNY, SELECT, D, Z, LDZ, SIGMA, &
                      WORKEV, &
                      BMAT, N, WHICH, NEV, TOL, &
                      RESID, NCV, V, LDV, IPARAM, IPNTR, WORKD, WORKL, LWORKL, &
                      RWORK, &
                      INFO)
      import :: I32, R64
      !> Compute eigenvectors? If .true., Z will contain them.
      logical, intent(in)          :: RVEC
      !> How many eigenvectors ('A' all, 'P' selected by SELECT).
      character(len=1), intent(in) :: HOWMNY
      !> Selection mask used when HOWMNY='P'.
      logical, intent(in)          :: SELECT(*)
      !> Output eigenvalues (complex).
      complex(R64), intent(out)    :: D(*)
      !> Leading dimension of Z.
      integer(I32), intent(in)     :: LDZ
      !> Output eigenvectors (columns) if RVEC.
      complex(R64), intent(out)    :: Z(LDZ, *)
      !> Shift for the implicit restart.
      complex(R64), intent(in)     :: SIGMA
      !> Work array for eigenvalue computations.
      complex(R64), intent(out)    :: WORKEV(*)
      !> Type of B in Ax = λ Bx ('I' identity, 'G' general).
      character(len=1), intent(in) :: BMAT
      !> Dimension of the eigenproblem.
      integer(I32), intent(in)     :: N
      !> Which eigenvalues to compute.
      character(len=2), intent(in) :: WHICH
      !> Number of eigenvalues to compute.
      integer(I32), intent(in)     :: NEV
      !> Convergence tolerance.
      real(R64), intent(in)        :: TOL
      !> Residual / initial vector (modified by ARPACK).
      complex(R64), intent(inout)  :: RESID(*)
      !> Number of Arnoldi vectors (subspace size).
      integer(I32), intent(in)     :: NCV
      !> Leading dimension of V.
      integer(I32), intent(in)     :: LDV
      !> Arnoldi basis vectors.
      complex(R64), intent(inout)  :: V(LDV, *)
      !> Integer parameters per ARPACK documentation.
      integer(I32), intent(inout)  :: IPARAM(11)
      !> Reverse-communication pointers.
      integer(I32), intent(inout)  :: IPNTR(14)
      !> Work array used in reverse communication.
      complex(R64), intent(inout)  :: WORKD(3 * N)
      !> Length of WORKL.
      integer(I32), intent(in)     :: LWORKL
      !> Internal work array used by ARPACK.
      complex(R64), intent(inout)  :: WORKL(LWORKL)
      !> Real work array.
      real(R64), intent(out)       :: RWORK(*)
      !> Status/diagnostic info (0 on success).
      integer(I32), intent(inout)  :: INFO
    end subroutine

    subroutine DSEUPD(RVEC, HOWMNY, SELECT, D, Z, LDZ, SIGMA, &
                      BMAT, N, WHICH, NEV, TOL, &
                      RESID, NCV, V, LDV, IPARAM, IPNTR, WORKD, WORKL, LWORKL, &
                      INFO)
      import :: I32, R64
      !> Compute eigenvectors? If .true., Z will contain them.
      logical, intent(in)          :: RVEC
      !> How many eigenvectors ('A' all, 'P' selected by SELECT).
      character(len=1), intent(in) :: HOWMNY
      !> Selection mask used when HOWMNY='P'.
      logical, intent(in)          :: SELECT(*)
      !> Output eigenvalues (real).
      real(R64), intent(out)    :: D(*)
      !> Leading dimension of Z.
      integer(I32), intent(in)     :: LDZ
      !> Output eigenvectors (columns) if RVEC.
      real(R64), intent(out)    :: Z(LDZ, *)
      !> Shift for the implicit restart.
      real(R64), intent(in)     :: SIGMA
      !> Type of B in Ax = λ Bx ('I' identity, 'G' general).
      character(len=1), intent(in) :: BMAT
      !> Dimension of the eigenproblem.
      integer(I32), intent(in)     :: N
      !> Which eigenvalues to compute.
      character(len=2), intent(in) :: WHICH
      !> Number of eigenvalues to compute.
      integer(I32), intent(in)     :: NEV
      !> Convergence tolerance.
      real(R64), intent(in)        :: TOL
      !> Residual / initial vector (modified by ARPACK).
      real(R64), intent(inout)  :: RESID(*)
      !> Number of Arnoldi vectors (subspace size).
      integer(I32), intent(in)     :: NCV
      !> Leading dimension of V.
      integer(I32), intent(in)     :: LDV
      !> Arnoldi basis vectors.
      real(R64), intent(inout)  :: V(LDV, *)
      !> Integer parameters per ARPACK documentation.
      integer(I32), intent(inout)  :: IPARAM(11)
      !> Reverse-communication pointers.
      integer(I32), intent(inout)  :: IPNTR(14)
      !> Work array used in reverse communication.
      real(R64), intent(inout)  :: WORKD(3 * N)
      !> Length of WORKL.
      integer(I32), intent(in)     :: LWORKL
      !> Internal work array used by ARPACK.
      real(R64), intent(inout)  :: WORKL(LWORKL)
      !> Status/diagnostic info (0 on success).
      integer(I32), intent(inout)  :: INFO
    end subroutine

  end interface

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine WrapZNAUPD(evals, evecs, nFound, ApplyMatOnVec, dim, evecsQ, nEvals, which, bmat, nKry)
    use M_Utils_Combinatorics

    !> Output eigenvalues (real). Allocated to size nFound on exit.
    real(R64), intent(out), allocatable              :: evals(:)
    !> Output eigenvectors (complex), size (dim, nFound) if evecsQ.
    complex(R64), intent(out), allocatable           :: evecs(:, :)
    !> Number of converged eigenpairs on exit.
    integer(I32), intent(out)                        :: nFound
    !> Callback to apply the matrix/operator: y = A*x.
    procedure(I_ApplyMatOnVecComplex)                 :: ApplyMatOnVec
    !> Problem dimension.
    integer(I32), intent(in)                         :: dim
    !> Whether to return eigenvectors.
    logical, intent(in)                              :: evecsQ
    !> Requested number of eigenvalues.
    integer(I32), intent(in)                         :: nEvals
    !> Which eigenvalues to compute (e.g., 'LM').
    character(2), intent(in)                         :: which
    !> Type of B in Ax = λ Bx ('I' identity, 'G' general).
    character(1), intent(in)                         :: bmat
    !> Arnoldi subspace size (NCV).
    integer(I32), intent(in)                         :: nKry

    complex(R64), allocatable :: workl(:), d(:), resid(:)
    complex(R64), allocatable :: v(:, :), workd(:), workev(:)
    real(R64), allocatable :: rwork(:)
    logical, allocatable :: select(:)
    integer(I32), allocatable :: permut(:)
    integer(I32) :: iparam(11), ipntr(14)
    integer(I32) :: ldv, ido, lworkl, info, ierr, j, maxitr, mode, ishfts
    integer(I32) :: count
    real(R64) :: tol
    complex(R64) :: sigma
    real(R64), parameter :: zero = 0.0

    ldv = dim
    lworkl = 3 * nKry**2 + 5 * nKry

    allocate (d(nKry))
    allocate (resid(dim))
    allocate (select(nKry))
    allocate (workl(lworkl))
    allocate (workd(3 * dim))
    allocate (v(dim, nKry))
    allocate (rwork(nKry))
    allocate (workev(2 * nKry))

    tol = zero
    info = 0
    ido = 0
    ishfts = 1
    maxitr = 300
    mode = 1
    iparam(1) = ishfts
    iparam(3) = maxitr
    iparam(7) = mode

    count = 0
    do
      count = count + 1

      call ZNAUPD(ido, bmat, dim, which, nEvals, tol, resid, nKry, v, ldv, iparam, ipntr, workd, workl, &
                  lworkl, rwork, info)

      if (ido .ne. -1 .and. ido .ne. 1) then
        exit
      end if

      call ApplyMatOnVec(workd(ipntr(2):ipntr(2) + dim - 1), &
                         workd(ipntr(1):ipntr(1) + dim - 1))

    end do

    !  Either we have convergence or there is an error.
    if (info < 0) then
      write (*, *) "Fatal error with ZNAUPD, INFO = ", info
      error stop
    end if

    call ZNEUPD(evecsQ, "A", select, d, v, ldv, sigma, workev, bmat, dim, which, &
                nEvals, tol, resid, nKry, v, ldv, iparam, ipntr, workd, workl, lworkl, rwork, ierr)

    if (ierr .ne. 0) then
      write (*, *) "Fatal error with ZNEUPD, IERR = ", ierr
      error stop

    else

      nFound = iparam(5) ! number of converged eigenvalues

      if (.not. allocated(evals)) allocate (evals(nFound))

      do j = 1, nFound
        evals(j) = dble(d(j))
      end do
      call Combinatorics_SortRealArray(evals, nFound, permut_=permut)

      if (evecsQ) then
        if (.not. allocated(evecs)) allocate (evecs(dim, nFound))

        do j = 1, nFound
          evecs(:, j) = v(:, permut(j))
        end do

      end if
    end if

    if (info .eq. 1) then
      write (*, *) "Maximum number of iterations reached, try increasing nKry."
      error stop
    else if (info .eq. 3) then
      write (*, *) "No shifts could be applied during implicit"//" Arnoldi update, try increasing nKry."
      error stop
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine WrapDSAUPD(evals, evecs, nFound, ApplyMatOnVec, dim, evecsQ, nEvals, which, bmat, nKry)
    use M_Utils_Combinatorics

    !> Output eigenvalues (real). Allocated to size nFound on exit.
    real(R64), intent(out), allocatable              :: evals(:)
    !> Output eigenvectors (real), size (dim, nFound) if evecsQ.
    real(R64), intent(out), allocatable           :: evecs(:, :)
    !> Number of converged eigenpairs on exit.
    integer(I32), intent(out)                        :: nFound
    !> Callback to apply the matrix/operator: y = A*x.
    procedure(I_ApplyMatOnVecReal)                 :: ApplyMatOnVec
    !> Problem dimension.
    integer(I32), intent(in)                         :: dim
    !> Whether to return eigenvectors.
    logical, intent(in)                              :: evecsQ
    !> Requested number of eigenvalues.
    integer(I32), intent(in)                         :: nEvals
    !> Which eigenvalues to compute (e.g., 'LM').
    character(2), intent(in)                         :: which
    !> Type of B in Ax = λ Bx ('I' identity, 'G' general).
    character(1), intent(in)                         :: bmat
    !> Arnoldi subspace size (NCV).
    integer(I32), intent(in)                         :: nKry

    real(R64), allocatable :: workl(:), d(:), resid(:)
    real(R64), allocatable :: v(:, :), workd(:), workev(:)
    real(R64), allocatable :: rwork(:)
    logical, allocatable :: select(:)
    integer(I32), allocatable :: permut(:)
    integer(I32) :: iparam(11), ipntr(14)
    integer(I32) :: ldv, ido, lworkl, info, ierr, j, maxitr, mode, ishfts
    integer(I32) :: count
    real(R64) :: tol
    real(R64) :: sigma
    real(R64), parameter :: zero = 0.0

    ldv = dim
    lworkl = 3 * nKry**2 + 5 * nKry

    allocate (d(nKry))
    allocate (resid(dim))
    allocate (select(nKry))
    allocate (workl(lworkl))
    allocate (workd(3 * dim))
    allocate (v(dim, nKry))
    allocate (rwork(nKry))
    allocate (workev(2 * nKry))

    tol = zero
    info = 0
    ido = 0
    ishfts = 1
    maxitr = 300
    mode = 1
    iparam(1) = ishfts
    iparam(3) = maxitr
    iparam(7) = mode

    count = 0
    do
      count = count + 1

      call DSAUPD(ido, bmat, dim, which, nEvals, tol, resid, nKry, v, ldv, iparam, ipntr, workd, workl, &
                  lworkl, info)

      if (ido .ne. -1 .and. ido .ne. 1) then
        exit
      end if

      call ApplyMatOnVec(workd(ipntr(2):ipntr(2) + dim - 1), &
                         workd(ipntr(1):ipntr(1) + dim - 1))

    end do

    !  Either we have convergence or there is an error.
    if (info < 0) then
      write (*, *) "Fatal error with DSAUPD, INFO = ", info
      error stop
    end if

    call DSEUPD(evecsQ, "All", select, d, v, ldv, sigma, bmat, dim, which, &
                nEvals, tol, resid, nKry, v, ldv, iparam, ipntr, workd, workl, lworkl, ierr)

    if (ierr .ne. 0) then
      write (*, *) "Fatal error with DSEUPD, IERR = ", ierr
      error stop

    else

      nFound = iparam(5) ! number of converged eigenvalues

      if (.not. allocated(evals)) allocate (evals(nFound))

      do j = 1, nFound
        evals(j) = dble(d(j))
      end do
      call Combinatorics_SortRealArray(evals, nFound, permut_=permut)

      if (evecsQ) then
        if (.not. allocated(evecs)) allocate (evecs(dim, nFound))

        do j = 1, nFound
          evecs(:, j) = v(:, permut(j))
        end do

      end if
    end if

    if (info .eq. 1) then
      write (*, *) "Maximum number of iterations reached, try increasing nKry."
      error stop
    else if (info .eq. 3) then
      write (*, *) "No shifts could be applied during implicit"//" Arnoldi update, try increasing nKry."
      error stop
    end if

  end subroutine

end module
