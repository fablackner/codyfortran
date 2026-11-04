!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Systems and callbacks of T_Utils_Gmres (module procedures: callbacks need
!> no trampolines), dimension N_DIM, in the diagonal metric
!> <f|g> = sum conj(f_i) w_i g_i (w in (0.5, 1.5)):
!> - A = M + N, non-Hermitian: M = diag(dM) complex (the preconditioner),
!>   N a dense random complex perturbation; the Arnoldi operator in the form
!>   A M^-1 v = v + N z (z = M^-1 v), A applied explicitly as dM y + N y
!> - A = M (the preconditioner is exact)
!> - A = 2 with the identity preconditioner, Euclidean inner product
module M_T_Utils_Gmres
  use M_Utils_Types
  use M_Utils_UnusedVariables

  implicit none

  !> Dimension of the test systems.
  integer(I32), parameter :: N_DIM = 40
  !> Amplitude of the entries of N (real and imaginary parts in (-1, 1)).
  real(R64), parameter :: OFF_AMP = 0.05_R64

  !> Diagonal dM of M, the perturbation N and the metric weights.
  complex(R64), allocatable :: dM(:), nMat(:, :)
  real(R64), allocatable :: wts(:)

  !> State of the deterministic generator (Park-Miller).
  integer(I64) :: seed = 12345_I64

contains

  !> Uniform pseudo-random number in (-1, 1).
  function RandomSigned() result(res)
    real(R64) :: res

    seed = mod(48271_I64 * seed, 2147483647_I64)
    res = 2.0_R64 * real(seed, kind=R64) / 2147483647.0_R64 - 1.0_R64

  end function

  !> dM_i = 2 + i/N + 0.5 i r_i (r in (-1, 1)), N with random entries of
  !> amplitude OFF_AMP, w in (0.5, 1.5).
  subroutine BuildSystem
    integer(I32) :: i, j
    real(R64) :: re, im

    allocate (dM(N_DIM), nMat(N_DIM, N_DIM), wts(N_DIM))
    do i = 1, N_DIM
      dM(i) = cmplx(2.0_R64 + real(i, kind=R64) / real(N_DIM, kind=R64), 0.5_R64 * RandomSigned(), kind=R64)
    end do
    do j = 1, N_DIM
      do i = 1, N_DIM
        re = RandomSigned()
        im = RandomSigned()
        nMat(i, j) = OFF_AMP * cmplx(re, im, kind=R64)
      end do
    end do
    do i = 1, N_DIM
      wts(i) = 1.0_R64 + 0.5_R64 * RandomSigned()
    end do

  end subroutine

  !> z = M^-1 v.
  subroutine ApplyPrecM(z, v)
    complex(R64), intent(out), contiguous :: z(:)
    complex(R64), intent(in), contiguous :: v(:)

    z(:) = v(:) / dM(1:size(v))

  end subroutine

  !> w = (M + N) z = v + N z for z = M^-1 v.
  subroutine ApplyMatOnPrecVecMN(w, v, z)
    complex(R64), intent(out), contiguous :: w(:)
    complex(R64), intent(in), contiguous :: v(:)
    complex(R64), intent(in), contiguous :: z(:)

    w(:) = v(:) + matmul(nMat, z)

  end subroutine

  !> w = (M + N) y.
  subroutine ApplyMatMN(w, y)
    complex(R64), intent(out), contiguous :: w(:)
    complex(R64), intent(in), contiguous :: y(:)

    w(:) = dM(1:size(y)) * y(:) + matmul(nMat, y)

  end subroutine

  !> w = M z = v for z = M^-1 v.
  subroutine ApplyMatOnPrecVecM(w, v, z)
    complex(R64), intent(out), contiguous :: w(:)
    complex(R64), intent(in), contiguous :: v(:)
    complex(R64), intent(in), contiguous :: z(:)

    if (.false.) call UnusedVariables_Mark(z)

    w(:) = v(:)

  end subroutine

  !> w = M y.
  subroutine ApplyMatM(w, y)
    complex(R64), intent(out), contiguous :: w(:)
    complex(R64), intent(in), contiguous :: y(:)

    w(:) = dM(1:size(y)) * y(:)

  end subroutine

  !> z = v (identity preconditioner).
  subroutine ApplyPrecIdentity(z, v)
    complex(R64), intent(out), contiguous :: z(:)
    complex(R64), intent(in), contiguous :: v(:)

    z(:) = v(:)

  end subroutine

  !> w = 2 z.
  subroutine ApplyMatOnPrecVecTwo(w, v, z)
    complex(R64), intent(out), contiguous :: w(:)
    complex(R64), intent(in), contiguous :: v(:)
    complex(R64), intent(in), contiguous :: z(:)

    if (.false.) call UnusedVariables_Mark(v)

    w(:) = 2.0_R64 * z(:)

  end subroutine

  !> w = 2 y.
  subroutine ApplyMatTwo(w, y)
    complex(R64), intent(out), contiguous :: w(:)
    complex(R64), intent(in), contiguous :: y(:)

    w(:) = 2.0_R64 * y(:)

  end subroutine

  !> Metric <fConjg|f> = sum conj(fConjg_i) w_i f_i.
  pure function InnerProductW(fConjg, f) result(res)
    complex(R64) :: res
    complex(R64), intent(in), contiguous :: fConjg(:)
    complex(R64), intent(in), contiguous :: f(:)

    res = sum(conjg(fConjg(:)) * wts(1:size(f)) * f(:))

  end function

  !> Euclidean <fConjg|f> = sum conj(fConjg_i) f_i.
  pure function InnerProductEuclidean(fConjg, f) result(res)
    complex(R64) :: res
    complex(R64), intent(in), contiguous :: fConjg(:)
    complex(R64), intent(in), contiguous :: f(:)

    res = sum(conjg(fConjg(:)) * f(:))

  end function

  !> ||b - (M + N) y|| / ||b|| in the metric, recomputed here.
  function RelativeResidualMN(y, b) result(res)
    complex(R64), intent(in), contiguous :: y(:)
    complex(R64), intent(in), contiguous :: b(:)
    real(R64) :: res

    complex(R64), allocatable :: r(:)

    allocate (r(size(b)))
    call ApplyMatMN(r, y)
    r(:) = b(:) - r(:)
    res = sqrt(real(InnerProductW(r, r), kind=R64) / real(InnerProductW(b, b), kind=R64))

  end function

end module

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Restarted flexible GMRES of M_Utils_Gmres on small dense complex systems
!> (M_T_Utils_Gmres; tolerance 1e-12 unless stated):
!> - A = M + N (n = 40, non-Hermitian, diagonal metric), M^-1 as
!>   preconditioner, one cycle (m = 40): converged, the solution matches a
!>   direct solve (ZGESV) to 1e-10, and the returned relative residual is
!>   the true one (recomputed here)
!> - the restart path: m = 2 and up to 200 restarts give the same solution
!>   to 1e-10 with more iterations than one cycle holds
!> - maxRestarts = 0 with m = 2: not converged after one cycle of two
!>   iterations; the solver returns (no stop) with the true residual of its
!>   iterate, which lies between the tolerance and 1
!> - exact preconditioner (A = M): the true residual of y0 = M^-1 b already
!>   meets the tolerance, so y = M^-1 b bitwise without an iteration
!> - happy breakdown: A = 2, identity preconditioner, Euclidean inner
!>   product, b a unit vector: the first Arnoldi vector spans an invariant
!>   space, the new vector vanishes exactly, and after one iteration y = b/2
!>   bitwise with a zero true residual
!> - b = 0: y = 0, converged, no iteration
!> - the context: Gmres_CreateCtx stores the controls without a workspace,
!>   the first solve allocates it for size(b) (2m + 2 vectors), a solve of
!>   another size reallocates it, Gmres_DestroyCtx releases it
program T_Utils_Gmres
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_LapackLib
  use M_Utils_Gmres
  use M_T_Utils_Gmres
  use testdrive, only: check, error_type

  implicit none

  real(R64), parameter :: TOL = 1.0e-12_R64

  type(error_type), allocatable :: error
  type(T_Gmres_Ctx) :: ctx
  complex(R64), allocatable :: b(:), y(:), yRef(:), yPrec(:), aMat(:, :), rhs(:, :)
  real(R64) :: relRes, relResOneCycle, relResRestart, relResFail, errOneCycle, errRestart
  integer(I32) :: nIter, nIterOneCycle, nIterRestart, nIterFail, i
  logical :: convergedQ

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call BuildSystem

  allocate (b(N_DIM), y(N_DIM))
  do i = 1, N_DIM
    b(i) = cmplx(RandomSigned(), RandomSigned(), kind=R64)
  end do

  ! Direct solve of (M + N) y = b
  allocate (aMat(N_DIM, N_DIM), rhs(N_DIM, 1))
  aMat(:, :) = nMat(:, :)
  do i = 1, N_DIM
    aMat(i, i) = aMat(i, i) + dM(i)
  end do
  rhs(:, 1) = b(:)
  call LapackLib_SolveLinearSystem(aMat, rhs)
  yRef = rhs(:, 1)

  !------------------------------------
  ! context without workspace
  !------------------------------------

  call Gmres_CreateCtx(ctx, N_DIM, 0, TOL)

  call check(error, ctx % krylovDim .eq. N_DIM .and. ctx % maxRestarts .eq. 0 .and. ctx % tolerance <= TOL &
             .and. ctx % tolerance >= TOL .and. .not. allocated(ctx % basis), "Gmres_CreateCtx stores the controls only")
  if (allocated(error)) error stop "T_Utils_Gmres context failure"

  !------------------------------------
  ! A = M + N, one cycle
  !------------------------------------

  call Gmres_Solve(y, nIterOneCycle, relResOneCycle, convergedQ, ctx, b, ApplyPrecM, ApplyMatOnPrecVecMN, ApplyMatMN, &
                   InnerProductW)
  errOneCycle = maxval(abs(y - yRef)) / maxval(abs(yRef))

  call check(error, convergedQ .and. nIterOneCycle >= 1 .and. nIterOneCycle <= N_DIM, "one cycle converges")
  if (allocated(error)) error stop "T_Utils_Gmres one-cycle convergence failure"
  call check(error, relResOneCycle <= TOL, "relative residual within the tolerance")
  if (allocated(error)) error stop "T_Utils_Gmres one-cycle residual failure"
  call check(error, abs(relResOneCycle - RelativeResidualMN(y, b)), 0.0_R64, thr=1.0e-14_R64)
  if (allocated(error)) error stop "T_Utils_Gmres one-cycle returned residual failure"
  call check(error, errOneCycle, 0.0_R64, thr=1.0e-10_R64)
  if (allocated(error)) error stop "T_Utils_Gmres one-cycle solution vs ZGESV failure"
  call check(error, allocated(ctx % basis) .and. allocated(ctx % precBasis) .and. allocated(ctx % residual), &
             "the first solve allocates the workspace")
  if (allocated(error)) error stop "T_Utils_Gmres workspace failure"
  call check(error, size(ctx % basis, 1) .eq. N_DIM .and. size(ctx % basis, 2) .eq. N_DIM + 1 .and. &
             size(ctx % precBasis, 2) .eq. N_DIM .and. size(ctx % residual) .eq. N_DIM, "workspace of 2m + 2 vectors")
  if (allocated(error)) error stop "T_Utils_Gmres workspace size failure"

  !------------------------------------
  ! A = M + N, restart path (m = 2)
  !------------------------------------

  call Gmres_CreateCtx(ctx, 2, 200, TOL)

  call check(error, .not. allocated(ctx % basis), "Gmres_CreateCtx releases the workspace of an earlier context")
  if (allocated(error)) error stop "T_Utils_Gmres context reset failure"

  call Gmres_Solve(y, nIterRestart, relResRestart, convergedQ, ctx, b, ApplyPrecM, ApplyMatOnPrecVecMN, ApplyMatMN, &
                   InnerProductW)
  errRestart = maxval(abs(y - yRef)) / maxval(abs(yRef))

  call check(error, convergedQ .and. nIterRestart > 2, "the restarted solve converges after several cycles")
  if (allocated(error)) error stop "T_Utils_Gmres restart convergence failure"
  call check(error, abs(relResRestart - RelativeResidualMN(y, b)), 0.0_R64, thr=1.0e-14_R64)
  if (allocated(error)) error stop "T_Utils_Gmres restart returned residual failure"
  call check(error, errRestart, 0.0_R64, thr=1.0e-10_R64)
  if (allocated(error)) error stop "T_Utils_Gmres restart solution vs ZGESV failure"
  call check(error, size(ctx % basis, 2) .eq. 3 .and. size(ctx % precBasis, 2) .eq. 2, "workspace of m = 2")
  if (allocated(error)) error stop "T_Utils_Gmres restart workspace failure"

  !------------------------------------
  ! no restart: not converged, no stop
  !------------------------------------

  call Gmres_CreateCtx(ctx, 2, 0, TOL)
  call Gmres_Solve(y, nIterFail, relResFail, convergedQ, ctx, b, ApplyPrecM, ApplyMatOnPrecVecMN, ApplyMatMN, &
                   InnerProductW)

  call check(error, .not. convergedQ .and. nIterFail .eq. 2, "one cycle of two iterations does not converge")
  if (allocated(error)) error stop "T_Utils_Gmres maxRestarts = 0 failure"
  call check(error, relResFail > TOL .and. relResFail < 1.0_R64, "the residual is reduced, above the tolerance")
  if (allocated(error)) error stop "T_Utils_Gmres maxRestarts = 0 residual failure"
  call check(error, abs(relResFail - RelativeResidualMN(y, b)) / relResFail, 0.0_R64, thr=1.0e-12_R64)
  if (allocated(error)) error stop "T_Utils_Gmres maxRestarts = 0 returned residual failure"

  !------------------------------------
  ! exact preconditioner (A = M)
  !------------------------------------

  call Gmres_Solve(y, nIter, relRes, convergedQ, ctx, b, ApplyPrecM, ApplyMatOnPrecVecM, ApplyMatM, InnerProductW)

  allocate (yPrec(N_DIM))
  call ApplyPrecM(yPrec, b)

  call check(error, convergedQ .and. nIter .eq. 0 .and. relRes <= TOL, "converged at the initial guess")
  if (allocated(error)) error stop "T_Utils_Gmres exact preconditioner failure"
  call check(error, maxval(abs(y - yPrec)) <= 0.0_R64, "y = M^-1 b bitwise")
  if (allocated(error)) error stop "T_Utils_Gmres exact preconditioner solution failure"

  !------------------------------------
  ! happy breakdown (A = 2, b a unit vector)
  !------------------------------------

  b(:) = (0.0_R64, 0.0_R64)
  b(7) = (1.0_R64, 0.0_R64)
  call Gmres_CreateCtx(ctx, 5, 0, TOL)
  call Gmres_Solve(y, nIter, relRes, convergedQ, ctx, b, ApplyPrecIdentity, ApplyMatOnPrecVecTwo, ApplyMatTwo, &
                   InnerProductEuclidean)

  call check(error, convergedQ .and. nIter .eq. 1 .and. .not. (relRes > 0.0_R64), "exact after one iteration")
  if (allocated(error)) error stop "T_Utils_Gmres happy breakdown failure"
  call check(error, maxval(abs(y - 0.5_R64 * b)) <= 0.0_R64, "y = b/2 bitwise")
  if (allocated(error)) error stop "T_Utils_Gmres happy breakdown solution failure"

  !------------------------------------
  ! b = 0
  !------------------------------------

  b(:) = (0.0_R64, 0.0_R64)
  y(:) = (1.0_R64, 1.0_R64)
  call Gmres_Solve(y, nIter, relRes, convergedQ, ctx, b, ApplyPrecM, ApplyMatOnPrecVecMN, ApplyMatMN, InnerProductW)

  call check(error, convergedQ .and. nIter .eq. 0 .and. .not. (relRes > 0.0_R64) .and. maxval(abs(y)) <= 0.0_R64, &
             "b = 0 gives y = 0")
  if (allocated(error)) error stop "T_Utils_Gmres zero right-hand side failure"

  !------------------------------------
  ! another size reallocates the workspace; destroy
  !------------------------------------

  deallocate (y)
  allocate (y(N_DIM / 2))
  call Gmres_Solve(y, nIter, relRes, convergedQ, ctx, yRef(1:N_DIM / 2), ApplyPrecM, ApplyMatOnPrecVecM, ApplyMatM, &
                   InnerProductW)

  call check(error, convergedQ .and. size(ctx % basis, 1) .eq. N_DIM / 2 .and. size(ctx % residual) .eq. N_DIM / 2, &
             "a solve of another size reallocates the workspace")
  if (allocated(error)) error stop "T_Utils_Gmres reallocation failure"

  call Gmres_DestroyCtx(ctx)

  call check(error, .not. (allocated(ctx % basis) .or. allocated(ctx % precBasis) .or. allocated(ctx % residual)) &
             .and. ctx % krylovDim .eq. 0, "Gmres_DestroyCtx releases the workspace")
  if (allocated(error)) error stop "T_Utils_Gmres destroy failure"

  !==========================================
  call Say_Section("printout results")
  !==========================================

  print *
  write (*, '(A, I0, A, ES10.3, A, ES10.3)') " A = M + N, one cycle (m = 40):    ", nIterOneCycle, &
    " iterations, relative residual ", relResOneCycle, ", error vs ZGESV ", errOneCycle
  write (*, '(A, I0, A, ES10.3, A, ES10.3)') " A = M + N, restarted (m = 2):     ", nIterRestart, &
    " iterations, relative residual ", relResRestart, ", error vs ZGESV ", errRestart
  write (*, '(A, I0, A, ES10.3)') " A = M + N, m = 2, no restart:     ", nIterFail, &
    " iterations, not converged, relative residual ", relResFail
  print *
  print *, "all GMRES checks passed"

  !==========================================
  call Say_Goodbye
  !==========================================

end program
