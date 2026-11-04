! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Restarted, flexible, right-preconditioned GMRES(m) for a complex linear
!> system A y = b in an inner product of the caller (e.g. Grid_InnerProduct,
!> the grid metric), matrix-free: the operator, the preconditioner and the
!> inner product are callbacks.
!>
!> Algorithm (flexible GMRES, Saad 1993):
!> - Initial guess y0 = M^-1 b.
!> - Every convergence decision, the first one included, is made on the true
!>   residual b - A y with A applied explicitly (ApplyMatOnVec), so the
!>   round-off of M M^-1 = 1 or an inexact preconditioner cannot be accepted
!>   unchecked: converged if ||b - A y|| <= tolerance ||b|| (norm of the
!>   inner product). It is computed for y0 and at the end of every cycle.
!> - A cycle builds the Arnoldi basis v_1 .. v_m of A M^-1 from the
!>   residual with modified Gram-Schmidt. The preconditioned vectors
!>   z_j = M^-1 v_j are kept (flexible form), so the update y += Z c needs no
!>   extra preconditioner solve, and the operator is applied to z_j with v_j
!>   given (ApplyMatOnPrecVec), so that A = M + N can be applied as
!>   v_j + N z_j without a second solve. Complex Givens rotations reduce the
!>   Hessenberg matrix; the cycle ends after m iterations, once the estimated
!>   residual |g(j+1)| <= tolerance ||b||, or when the new Arnoldi vector
!>   vanishes (happy breakdown: the solution lies in the Krylov space).
!> - At most 1 + maxRestarts cycles. The solver neither prints nor stops: it
!>   returns convergedQ, the relative true residual and the number of
!>   Arnoldi iterations. b = 0 gives y = 0 (converged, no iteration).
!>
!> Context (T_Gmres_Ctx, owned by the caller): the controls and the
!> workspace of 2m + 2 vectors (the Arnoldi vectors, their preconditioned
!> images and the residual). The first solve allocates the workspace for
!> size(b) and later solves keep it (its pages stay mapped from call to
!> call); a solve of another size reallocates it. A solve uses only its
!> context, its arguments and the callbacks: solves with different contexts
!> may run concurrently, one context serves one solve at a time.
!>
!> The callbacks are module procedures of the caller (internal procedures as
!> actual arguments need trampolines, i.e. an executable stack), which may
!> read data the caller prepared for the solve. InnerProduct has the
!> characteristics of I_Grid_InnerProduct, so Grid_InnerProduct can be
!> passed directly:
!>
!>   call Gmres_CreateCtx(ctx, krylovDim, maxRestarts, tolerance)
!>   call Gmres_Solve(y, nIterations, relativeResidual, convergedQ, ctx, b, &
!>                    ApplyPrecOnVec, ApplyMatOnPrecVec, ApplyMatOnVec, Grid_InnerProduct)
!>
!> The vector passes are OpenMP parallel when the solve is called outside a
!> parallel region (serial inside one). User: the Cayley step of the exact
!> field-free part of Method on grid.fem
!> (S_Method_Mb_OrbBased_FieldFreeCayley).
module M_Utils_Gmres
  use M_Utils_Types

  implicit none

  private
  public :: T_Gmres_Ctx
  public :: Gmres_CreateCtx
  public :: Gmres_DestroyCtx
  public :: Gmres_Solve
  public :: I_Gmres_ApplyPrecOnVec
  public :: I_Gmres_ApplyMatOnPrecVec
  public :: I_Gmres_ApplyMatOnVec
  public :: I_Gmres_InnerProduct

  !> Controls and workspace of the solver (Gmres_CreateCtx).
  type :: T_Gmres_Ctx
    !> Restart length m (Arnoldi iterations per cycle).
    integer(I32) :: krylovDim = 0
    !> Restarts: at most 1 + maxRestarts cycles.
    integer(I32) :: maxRestarts = 0
    !> Relative true residual ||b - A y|| / ||b|| to reach.
    real(R64) :: tolerance = 0.0_R64
    !> Workspace (allocated by the first Gmres_Solve, kept): Arnoldi vectors
    !> (n x m+1), their preconditioned images M^-1 v_j (n x m) and the
    !> residual (n). Pages are touched only for the iterations done.
    complex(R64), allocatable :: basis(:, :)
    complex(R64), allocatable :: precBasis(:, :)
    complex(R64), allocatable :: residual(:)
  end type

  abstract interface
    !> z = M^-1 v: the right preconditioner.
    subroutine I_Gmres_ApplyPrecOnVec(z, v)
      import :: R64
      complex(R64), intent(out), contiguous :: z(:)
      complex(R64), intent(in), contiguous :: v(:)
    end subroutine

    !> w = A z for the preconditioned vector z = M^-1 v (the Arnoldi
    !> operator A M^-1 v); v is given so that an operator A = M + N can
    !> form v + N z.
    subroutine I_Gmres_ApplyMatOnPrecVec(w, v, z)
      import :: R64
      complex(R64), intent(out), contiguous :: w(:)
      complex(R64), intent(in), contiguous :: v(:)
      complex(R64), intent(in), contiguous :: z(:)
    end subroutine

    !> w = A y, applied explicitly (the true residuals).
    subroutine I_Gmres_ApplyMatOnVec(w, y)
      import :: R64
      complex(R64), intent(out), contiguous :: w(:)
      complex(R64), intent(in), contiguous :: y(:)
    end subroutine

    !> Inner product <fConjg|f> (the characteristics of I_Grid_InnerProduct).
    pure function I_Gmres_InnerProduct(fConjg, f) result(res)
      import :: R64
      complex(R64) :: res
      complex(R64), intent(in), contiguous :: fConjg(:)
      complex(R64), intent(in), contiguous :: f(:)
    end function
  end interface

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Stores the controls in a fresh context (the workspace of an earlier
  !> context is released; the first solve allocates it).
  subroutine Gmres_CreateCtx(ctx, krylovDim, maxRestarts, tolerance)
    type(T_Gmres_Ctx), intent(out) :: ctx
    !> Restart length m >= 1.
    integer(I32), intent(in) :: krylovDim
    !> Restarts >= 0.
    integer(I32), intent(in) :: maxRestarts
    !> Relative true residual > 0.
    real(R64), intent(in) :: tolerance

    if (krylovDim < 1) error stop "Gmres_CreateCtx: krylovDim must be >= 1"
    if (maxRestarts < 0) error stop "Gmres_CreateCtx: maxRestarts must be >= 0"
    if (.not. (tolerance > 0.0_R64)) error stop "Gmres_CreateCtx: tolerance must be > 0"

    ctx % krylovDim = krylovDim
    ctx % maxRestarts = maxRestarts
    ctx % tolerance = tolerance

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Releases the workspace and resets the controls.
  subroutine Gmres_DestroyCtx(ctx)
    type(T_Gmres_Ctx), intent(inout) :: ctx

    call ReleaseWorkspace(ctx)
    ctx % krylovDim = 0
    ctx % maxRestarts = 0
    ctx % tolerance = 0.0_R64

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> y = A^-1 b by right-preconditioned flexible GMRES(m) (module header):
  !> initial guess M^-1 b, preconditioned Arnoldi vectors kept (precBasis),
  !> modified Gram-Schmidt, complex Givens rotations; every convergence
  !> decision on the true residual (initial guess and end of every cycle).
  subroutine Gmres_Solve(y, nIterations, relativeResidual, convergedQ, ctx, b, &
                         ApplyPrecOnVec, ApplyMatOnPrecVec, ApplyMatOnVec, InnerProduct)
    !> Solution (same size as b).
    complex(R64), intent(out), contiguous :: y(:)
    !> Arnoldi iterations over all cycles.
    integer(I32), intent(out) :: nIterations
    !> Last true residual ||b - A y|| / ||b|| (0 for b = 0).
    real(R64), intent(out) :: relativeResidual
    !> Whether relativeResidual <= ctx % tolerance.
    logical, intent(out) :: convergedQ
    !> Controls and workspace (Gmres_CreateCtx).
    type(T_Gmres_Ctx), intent(inout) :: ctx
    !> Right-hand side.
    complex(R64), intent(in), contiguous :: b(:)
    procedure(I_Gmres_ApplyPrecOnVec) :: ApplyPrecOnVec
    procedure(I_Gmres_ApplyMatOnPrecVec) :: ApplyMatOnPrecVec
    procedure(I_Gmres_ApplyMatOnVec) :: ApplyMatOnVec
    procedure(I_Gmres_InnerProduct) :: InnerProduct

    complex(R64), allocatable :: hessenberg(:, :), rhsSmall(:), cosines(:), sines(:), coefficients(:)
    complex(R64) :: temp
    real(R64) :: bNorm, residualNorm, target, hNext
    integer(I32) :: m, iCycle, j, i, nUsed

    if (ctx % krylovDim < 1) error stop "Gmres_Solve: the context has no controls (Gmres_CreateCtx)"
    if (size(y) .ne. size(b)) error stop "Gmres_Solve: y and b must have the same size"
    call AllocateWorkspace(ctx, size(b))

    associate (basis => ctx % basis, precBasis => ctx % precBasis, residual => ctx % residual)

      m = ctx % krylovDim
      allocate (hessenberg(m + 1, m), rhsSmall(m + 1), cosines(m), sines(m), coefficients(m))

      bNorm = MetricNorm(b, InnerProduct)
      nIterations = 0

      ! A zero right-hand side has the solution zero
      if (.not. (bNorm > 0.0_R64)) then
        y(:) = (0.0_R64, 0.0_R64)
        relativeResidual = 0.0_R64
        convergedQ = .true.
        return
      end if
      target = ctx % tolerance * bNorm

      ! y0 = M^-1 b and its true residual r0 = b - A y0
      call ApplyPrecOnVec(y, b)
      call ApplyResidual(residual, b, y, ApplyMatOnVec)
      residualNorm = MetricNorm(residual, InnerProduct)

      convergedQ = residualNorm <= target
      do iCycle = 0, ctx % maxRestarts
        if (convergedQ) exit

        basis(:, 1) = residual(:)
        call Scale(basis(:, 1), 1.0_R64 / residualNorm)
        rhsSmall(:) = (0.0_R64, 0.0_R64)
        rhsSmall(1) = cmplx(residualNorm, 0.0_R64, kind=R64)
        nUsed = 0

        do j = 1, m
          nIterations = nIterations + 1
          nUsed = j

          ! w = A M^-1 v_j = A z_j with z_j = M^-1 v_j, orthogonalized
          ! against v_1 .. v_j (modified Gram-Schmidt)
          call ApplyPrecOnVec(precBasis(:, j), basis(:, j))
          call ApplyMatOnPrecVec(basis(:, j + 1), basis(:, j), precBasis(:, j))
          do i = 1, j
            hessenberg(i, j) = InnerProduct(basis(:, i), basis(:, j + 1))
            call Axpy(basis(:, j + 1), -hessenberg(i, j), basis(:, i))
          end do
          hNext = MetricNorm(basis(:, j + 1), InnerProduct)
          hessenberg(j + 1, j) = cmplx(hNext, 0.0_R64, kind=R64)

          ! Previous rotations on the new column, then the rotation that
          ! annihilates its subdiagonal entry
          do i = 1, j - 1
            temp = cosines(i) * hessenberg(i, j) + sines(i) * hessenberg(i + 1, j)
            hessenberg(i + 1, j) = -conjg(sines(i)) * hessenberg(i, j) + cosines(i) * hessenberg(i + 1, j)
            hessenberg(i, j) = temp
          end do
          call FillGivensRotation(cosines(j), sines(j), hessenberg(j, j), hessenberg(j + 1, j))
          hessenberg(j, j) = cosines(j) * hessenberg(j, j) + sines(j) * hessenberg(j + 1, j)
          hessenberg(j + 1, j) = (0.0_R64, 0.0_R64)
          rhsSmall(j + 1) = -conjg(sines(j)) * rhsSmall(j)
          rhsSmall(j) = cosines(j) * rhsSmall(j)

          ! Estimated residual |g(j+1)|; a vanishing new vector ends the
          ! Krylov space (the solution lies in it)
          if (abs(rhsSmall(j + 1)) <= target .or. .not. (hNext > 0.0_R64)) exit
          call Scale(basis(:, j + 1), 1.0_R64 / hNext)
        end do

        ! Least-squares solution of the triangular system, y += Z c
        do i = nUsed, 1, -1
          coefficients(i) = (rhsSmall(i) - sum(hessenberg(i, i + 1:nUsed) * coefficients(i + 1:nUsed))) / hessenberg(i, i)
        end do
        do i = 1, nUsed
          call Axpy(y, coefficients(i), precBasis(:, i))
        end do

        ! True residual of the cycle
        call ApplyResidual(residual, b, y, ApplyMatOnVec)
        residualNorm = MetricNorm(residual, InnerProduct)
        convergedQ = residualNorm <= target
      end do

      relativeResidual = residualNorm / bNorm

    end associate

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Allocates the workspace for vectors of size n and the Krylov dimension
  !> of the context, unless it exists with these sizes.
  subroutine AllocateWorkspace(ctx, n)
    type(T_Gmres_Ctx), intent(inout) :: ctx
    integer(I32), intent(in) :: n

    if (allocated(ctx % basis)) then
      if (size(ctx % basis, 1) .eq. n .and. size(ctx % basis, 2) .eq. ctx % krylovDim + 1) return
      call ReleaseWorkspace(ctx)
    end if

    allocate (ctx % basis(n, ctx % krylovDim + 1), ctx % precBasis(n, ctx % krylovDim))
    allocate (ctx % residual(n))

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Deallocates the workspace.
  subroutine ReleaseWorkspace(ctx)
    type(T_Gmres_Ctx), intent(inout) :: ctx

    if (allocated(ctx % basis)) deallocate (ctx % basis)
    if (allocated(ctx % precBasis)) deallocate (ctx % precBasis)
    if (allocated(ctx % residual)) deallocate (ctx % residual)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Complex Givens rotation [c s; -conjg(s) c] (c real) that maps (a, b) to
  !> (r, 0).
  pure subroutine FillGivensRotation(c, s, a, b)
    complex(R64), intent(out) :: c
    complex(R64), intent(out) :: s
    complex(R64), intent(in) :: a
    complex(R64), intent(in) :: b

    real(R64) :: absA, norm

    absA = abs(a)
    norm = hypot(absA, abs(b))
    if (.not. (norm > 0.0_R64)) then
      c = (1.0_R64, 0.0_R64)
      s = (0.0_R64, 0.0_R64)
    else if (.not. (absA > 0.0_R64)) then
      c = (0.0_R64, 0.0_R64)
      s = conjg(b) / abs(b)
    else
      c = cmplx(absA / norm, 0.0_R64, kind=R64)
      s = (a / absA) * conjg(b) / norm
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> r = b - A y with A applied explicitly (the true residual).
  subroutine ApplyResidual(r, b, y, ApplyMatOnVec)
    !$ use omp_lib, only: omp_in_parallel

    complex(R64), intent(out), contiguous :: r(:)
    complex(R64), intent(in), contiguous :: b(:)
    complex(R64), intent(in), contiguous :: y(:)
    procedure(I_Gmres_ApplyMatOnVec) :: ApplyMatOnVec

    integer(I32) :: i

    call ApplyMatOnVec(r, y)

    !$omp parallel do default(shared) private(i) schedule(static) if (.not. omp_in_parallel())
    do i = 1, size(r)
      r(i) = b(i) - r(i)
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> y = y + a x.
  subroutine Axpy(y, a, x)
    !$ use omp_lib, only: omp_in_parallel

    complex(R64), intent(inout), contiguous :: y(:)
    complex(R64), intent(in) :: a
    complex(R64), intent(in), contiguous :: x(:)

    integer(I32) :: i

    !$omp parallel do default(shared) private(i) schedule(static) if (.not. omp_in_parallel())
    do i = 1, size(y)
      y(i) = y(i) + a * x(i)
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> y = a y.
  subroutine Scale(y, a)
    !$ use omp_lib, only: omp_in_parallel

    complex(R64), intent(inout), contiguous :: y(:)
    real(R64), intent(in) :: a

    integer(I32) :: i

    !$omp parallel do default(shared) private(i) schedule(static) if (.not. omp_in_parallel())
    do i = 1, size(y)
      y(i) = a * y(i)
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Norm of the inner product, sqrt(<f|f>).
  function MetricNorm(f, InnerProduct) result(res)
    complex(R64), intent(in), contiguous :: f(:)
    procedure(I_Gmres_InnerProduct) :: InnerProduct
    real(R64) :: res

    res = sqrt(max(real(InnerProduct(f, f), kind=R64), 0.0_R64))

  end function

end module
