!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Operators, metric and preconditioners of T_DiagonalizerList_Lobpcg, all of
!> dimension N_DIM; except H itself, all are self-adjoint in the metric
!> ⟨f|g⟩ = Σ conj(f_i) w_i g_i (W = diag(w) positive), but not in the
!> Euclidean one:
!> - A = W⁻¹H with H complex Hermitian (dominant increasing diagonal 1..N plus
!>   small random off-diagonal entries), with a safeguarded diagonal
!>   preconditioner
!> - H itself, self-adjoint in the Euclidean metric (a solver without
!>   metricWeights), with the safeguarded diagonal preconditioner
!>   (diag(H) − θ)⁻¹
!> - AK = W^(-1/2) Q diag(dK) Q^H W^(1/2) with Q unitary (eigenvectors of a
!>   random Hermitian matrix) and the prescribed spectrum dK = 1, 2, 2, 2
!>   (a threefold cluster), 3, 3.37, 3.74, …; spectral functions of AK are
!>   cheap, which gives the fixed preconditioner (AK + 1)⁻¹ and the nearly
!>   exact shift-and-invert preconditioner (AK − θ + 1e-12)⁻¹
!> - AK + iε W⁻¹S with S real symmetric: slightly non-self-adjoint
module M_T_DiagonalizerList_Lobpcg
  use M_Utils_Types
  use M_Utils_UnusedVariables

  implicit none

  !> Dimension of the test operators.
  integer(I32), parameter :: N_DIM = 300
  !> Size ε of the anti-Hermitian perturbation iεW⁻¹S (entries of S in (−1, 1)).
  real(R64), parameter :: NON_HERMITIAN_AMP = 1e-6_R64
  !> Offset δ of the shift-and-invert preconditioner (AK − θ + δ)⁻¹.
  real(R64), parameter :: SHIFT_INVERT_OFFSET = 1e-12_R64

  !> Hermitian matrix H, metric weights W and the diagonal of A = W⁻¹H.
  complex(R64), allocatable :: hMat(:, :)
  real(R64), allocatable :: wts(:), diagA(:)

  !> AK: unitary Q (and Qᴴ), spectrum dK; S of the anti-Hermitian perturbation.
  complex(R64), allocatable :: qK(:, :), qKH(:, :)
  real(R64), allocatable :: dK(:), sK(:, :)

  !> State of the deterministic generator (Park–Miller).
  integer(I64) :: seed = 12345_I64

contains

  !> Uniform pseudo-random number in (−1, 1).
  function RandomSigned() result(res)
    real(R64) :: res

    seed = mod(48271_I64 * seed, 2147483647_I64)
    res = 2.0_R64 * real(seed, kind=R64) / 2147483647.0_R64 - 1.0_R64

  end function

  !> H = diag(1..N) + offAmp·(random Hermitian), W in (0.5, 1.5).
  subroutine BuildOperator(offAmp)
    real(R64), intent(in) :: offAmp

    integer(I32) :: i, j
    real(R64) :: re, im

    allocate (hMat(N_DIM, N_DIM), wts(N_DIM), diagA(N_DIM))
    do j = 1, N_DIM
      hMat(j, j) = cmplx(real(j, kind=R64), 0.0_R64, kind=R64)
      do i = 1, j - 1
        re = RandomSigned()
        im = RandomSigned()
        hMat(i, j) = offAmp * cmplx(re, im, kind=R64)
        hMat(j, i) = conjg(hMat(i, j))
      end do
    end do
    do i = 1, N_DIM
      wts(i) = 1.0_R64 + 0.5_R64 * RandomSigned()
    end do
    call UpdateDiagonal

  end subroutine

  !> AK (after BuildOperator, which sets the weights): Q from a random
  !> Hermitian matrix, dK = 1, 2, 2, 2, 3, 3.37, …; S real symmetric.
  subroutine BuildKnownSpectrumOperator
    use M_Utils_LapackLib

    complex(R64), allocatable :: rnd(:, :)
    real(R64), allocatable :: evalsRnd(:)
    integer(I32) :: i, j
    real(R64) :: re, im

    allocate (rnd(N_DIM, N_DIM), dK(N_DIM), sK(N_DIM, N_DIM))
    do j = 1, N_DIM
      rnd(j, j) = cmplx(RandomSigned(), 0.0_R64, kind=R64)
      sK(j, j) = RandomSigned()
      do i = 1, j - 1
        re = RandomSigned()
        im = RandomSigned()
        rnd(i, j) = cmplx(re, im, kind=R64)
        rnd(j, i) = conjg(rnd(i, j))
        sK(i, j) = RandomSigned()
        sK(j, i) = sK(i, j)
      end do
    end do
    call LapackLib_DiagonalizeGeneric(evalsRnd, qK, rnd, .true.)
    qKH = conjg(transpose(qK))

    dK(1) = 1.0_R64
    dK(2:4) = 2.0_R64
    do i = 5, N_DIM
      dK(i) = 3.0_R64 + 0.37_R64 * real(i - 5, kind=R64)
    end do

  end subroutine

  !> H ← H + amp·(random Hermitian).
  subroutine PerturbOperator(amp)
    real(R64), intent(in) :: amp

    integer(I32) :: i, j
    real(R64) :: re, im

    do j = 1, N_DIM
      hMat(j, j) = hMat(j, j) + amp * RandomSigned()
      do i = 1, j - 1
        re = RandomSigned()
        im = RandomSigned()
        hMat(i, j) = hMat(i, j) + amp * cmplx(re, im, kind=R64)
        hMat(j, i) = conjg(hMat(i, j))
      end do
    end do
    call UpdateDiagonal

  end subroutine

  !> diag(A) = diag(H)/w.
  subroutine UpdateDiagonal
    integer(I32) :: i

    do i = 1, N_DIM
      diagA(i) = real(hMat(i, i), kind=R64) / wts(i)
    end do

  end subroutine

  !> dState = W⁻¹ H state.
  subroutine ApplyA(dState, state, time)
    complex(R64), intent(out), contiguous, target :: dState(:)
    complex(R64), intent(in), contiguous, target :: state(:)
    real(R64), intent(in) :: time

    if (.false.) call UnusedVariables_Mark(time)

    dState(:) = matmul(hMat, state) / wts(:)

  end subroutine

  !> Diagonal preconditioner (diag(A) − shift)⁻¹, safeguarded.
  subroutine ApplyPrec(dState, state, shift, time)
    complex(R64), intent(out), contiguous, target :: dState(:)
    complex(R64), intent(in), contiguous, target :: state(:)
    real(R64), intent(in) :: shift
    real(R64), intent(in) :: time

    if (.false.) call UnusedVariables_Mark(time)

    dState(:) = state(:) / max(abs(diagA(:) - shift), 0.1_R64)

  end subroutine

  !> dState = H state (Euclidean self-adjoint).
  subroutine ApplyH(dState, state, time)
    complex(R64), intent(out), contiguous, target :: dState(:)
    complex(R64), intent(in), contiguous, target :: state(:)
    real(R64), intent(in) :: time

    if (.false.) call UnusedVariables_Mark(time)

    dState(:) = matmul(hMat, state)

  end subroutine

  !> Diagonal preconditioner (diag(H) − shift)⁻¹ of H, safeguarded.
  subroutine ApplyPrecH(dState, state, shift, time)
    complex(R64), intent(out), contiguous, target :: dState(:)
    complex(R64), intent(in), contiguous, target :: state(:)
    real(R64), intent(in) :: shift
    real(R64), intent(in) :: time

    integer(I32) :: i

    if (.false.) call UnusedVariables_Mark(time)

    do i = 1, N_DIM
      dState(i) = state(i) / max(abs(real(hMat(i, i), kind=R64) - shift), 0.1_R64)
    end do

  end subroutine

  !> dState = f(AK) state = W^(-1/2) Q diag(f) Qᴴ W^(1/2) state.
  subroutine ApplySpectralFunction(dState, state, f)
    complex(R64), intent(out), contiguous :: dState(:)
    complex(R64), intent(in), contiguous :: state(:)
    real(R64), intent(in), contiguous :: f(:)

    complex(R64), allocatable :: c(:)

    c = f(:) * matmul(qKH, sqrt(wts(:)) * state(:))
    dState(:) = matmul(qK, c) / sqrt(wts(:))

  end subroutine

  !> dState = AK state.
  subroutine ApplyAK(dState, state, time)
    complex(R64), intent(out), contiguous, target :: dState(:)
    complex(R64), intent(in), contiguous, target :: state(:)
    real(R64), intent(in) :: time

    if (.false.) call UnusedVariables_Mark(time)

    call ApplySpectralFunction(dState, state, dK)

  end subroutine

  !> dState = (AK + iε W⁻¹S) state: self-adjoint in the metric only up to ε.
  subroutine ApplyAKNonHermitian(dState, state, time)
    complex(R64), intent(out), contiguous, target :: dState(:)
    complex(R64), intent(in), contiguous, target :: state(:)
    real(R64), intent(in) :: time

    call ApplyAK(dState, state, time)
    dState(:) = dState(:) + cmplx(0.0_R64, NON_HERMITIAN_AMP, kind=R64) * matmul(sK, state) / wts(:)

  end subroutine

  !> Fixed preconditioner (AK + 1)⁻¹ (positive definite, like (T + σ)⁻¹).
  subroutine ApplyAKFixedPrec(dState, state, shift, time)
    complex(R64), intent(out), contiguous, target :: dState(:)
    complex(R64), intent(in), contiguous, target :: state(:)
    real(R64), intent(in) :: shift
    real(R64), intent(in) :: time

    if (.false.) call UnusedVariables_Mark(shift, time)

    call ApplySpectralFunction(dState, state, 1.0_R64 / (dK + 1.0_R64))

  end subroutine

  !> Nearly exact shift-and-invert preconditioner (AK − shift + δ)⁻¹: the
  !> preconditioned residual of a Ritz pair is x_j up to O(δ), i.e. it has
  !> numerically no component outside the current basis.
  subroutine ApplyAKShiftInvert(dState, state, shift, time)
    complex(R64), intent(out), contiguous, target :: dState(:)
    complex(R64), intent(in), contiguous, target :: state(:)
    real(R64), intent(in) :: shift
    real(R64), intent(in) :: time

    if (.false.) call UnusedVariables_Mark(time)

    call ApplySpectralFunction(dState, state, 1.0_R64 / (dK - shift + SHIFT_INVERT_OFFSET))

  end subroutine

  !> Metric ⟨fConjg|f⟩ = Σ conj(fConjg_i) w_i f_i.
  pure function InnerProductW(fConjg, f) result(res)
    complex(R64) :: res
    complex(R64), intent(in), contiguous :: fConjg(:)
    complex(R64), intent(in), contiguous :: f(:)

    res = sum(conjg(fConjg(:)) * wts(:) * f(:))

  end function

  !> ⟨fConjg|f⟩ in the metric of a solve: the weighted one (InnerProductW),
  !> or the Euclidean one Σ conj(fConjg_i) f_i if euclideanQ (a solver
  !> without metricWeights).
  pure function InnerProductOfSolve(fConjg, f, euclideanQ) result(res)
    complex(R64) :: res
    complex(R64), intent(in), contiguous :: fConjg(:)
    complex(R64), intent(in), contiguous :: f(:)
    logical, intent(in) :: euclideanQ

    if (euclideanQ) then
      res = sum(conjg(fConjg(:)) * f(:))
    else
      res = InnerProductW(fConjg, f)
    end if

  end function

  !> Lowest nEvals eigenvalues of H via ZHEEV (the Euclidean reference).
  subroutine ReferenceEvalsEuclidean(evalsRef, nEvals)
    use M_Utils_LapackLib

    real(R64), intent(out) :: evalsRef(:)
    integer(I32), intent(in) :: nEvals

    complex(R64), allocatable :: evecsDummy(:, :)
    real(R64), allocatable :: evalsAll(:)

    call LapackLib_DiagonalizeGeneric(evalsAll, evecsDummy, hMat, .false.)
    evalsRef(1:nEvals) = evalsAll(1:nEvals)

  end subroutine

  !> Lowest nEvals eigenvalues of A via ZHEEV of W^(-1/2) H W^(-1/2).
  subroutine ReferenceEvals(evalsRef, nEvals)
    use M_Utils_LapackLib

    real(R64), intent(out) :: evalsRef(:)
    integer(I32), intent(in) :: nEvals

    complex(R64), allocatable :: hSym(:, :), evecsDummy(:, :)
    real(R64), allocatable :: evalsAll(:)
    integer(I32) :: i, j

    allocate (hSym(N_DIM, N_DIM))
    do j = 1, N_DIM
      do i = 1, N_DIM
        hSym(i, j) = hMat(i, j) / sqrt(wts(i) * wts(j))
      end do
    end do
    call LapackLib_DiagonalizeGeneric(evalsAll, evecsDummy, hSym, .false.)
    evalsRef(1:nEvals) = evalsAll(1:nEvals)

  end subroutine

end module

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> LOBPCG backend (diagonalizerList.lobpcg*) on operators of dimension 300
!> that are self-adjoint only in a weighted metric, given to every solver but
!> the Euclidean one as its metricWeights:
!> - A = W⁻¹H without and with a preconditioner (safeguarded
!>   (diag(A) − θ)⁻¹): the lowest 5 eigenvalues match the dense reference
!>   (ZHEEV of W^(-1/2) H W^(-1/2)) to 1e-10, the eigenvectors are
!>   orthonormal in the metric to 1e-12 and their residuals, recomputed from
!>   scratch, are small; convergedQ, nFound = 5; the preconditioner reduces
!>   the iteration count substantially
!> - a repeated solve of an unchanged operator needs no iteration, and after
!>   a small perturbation of H (1e-4) the warm start converges in at most
!>   half the iterations of the cold start
!> - the Euclidean metric (a solver without metricWeights: unit weights): the
!>   Hermitian H itself with the preconditioner (diag(H) − θ)⁻¹; the same
!>   checks, against ZHEEV of H and in the Euclidean inner product
!> - AK (spectrum 1, 2, 2, 2, 3, …) with the preconditioner (AK + 1)⁻¹: the
!>   threefold cluster split by nEvals = 3 (λ3 = λ4), with nGuard 0 and 2
!>   (the block ends inside the cluster with nGuard 0); a start block seeded
!>   (SetStartVectors) with linearly dependent columns (four equal ones and
!>   two that differ by 1e-9 and 1e-7), which the solver completes with
!>   random columns
!> - AK with the nearly exact shift-and-invert preconditioner
!>   (AK − θ + 1e-12)⁻¹: every preconditioned residual is x_j up to 1e-12,
!>   i.e. carries no new direction; the raw residuals replace them
!>   (nRawResiduals > 0) and the solve converges (it used to stop after one
!>   iteration with Ritz values 13…18 instead of 1, 2, 2, 2, 3)
!> - AK + iεW⁻¹S (ε = 1e-6), self-adjoint only up to ε: the residuals reach a
!>   floor far above tol = 1e-10, and the stagnation exit stops the solve
!>   early (stalledQ, convergedQ false, nFound < 5, nIter ≤ 60 of maxIter
!>   1000) with Ritz values still within 1e-4 of the spectrum
program T_DiagonalizerList_Lobpcg
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_DiagonalizerList
  use M_DiagonalizerList_Lobpcg, only: T_DiagonalizerList_E_Lobpcg
  use M_T_DiagonalizerList_Lobpcg
  use testdrive, only: check, error_type

  implicit none

  integer(I32), parameter :: N_EVALS = 5
  integer(I32), parameter :: N_EVALS_CLUSTER = 3
  !> Largest accepted iteration count of the stagnating solve (maxIter 1000)
  integer(I32), parameter :: MAX_ITER_STALL = 60
  type(error_type), allocatable :: error
  type(T_DiagonalizerList_FabricateInput) :: diagonalizerListInput(8)
  character(len=:), allocatable :: jsonFileName
  integer(I32) :: nIterPlain, nIterPrec, nIterSame, nIterWarm, nIterCluster0, nIterCluster2
  integer(I32) :: nIterSeeded, nIterShiftInvert, nIterEuclidean, i, j, iD
  real(R64) :: devStall
  complex(R64), allocatable :: startVectors(:, :)

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  jsonFileName = "test/componentTests/DiagonalizerList/Lobpcg/T_DiagonalizerList_Lobpcg.json"
  call Json_LoadJsonFile(manualFileName_=jsonFileName, relativeToProjectDirQ_=.true.)

  call BuildOperator(offAmp=0.05_R64)
  call BuildKnownSpectrumOperator

  ! 1: A without preconditioner, 2: A with the diagonal preconditioner
  diagonalizerListInput(1) % ApplyMatOnVec => ApplyA
  diagonalizerListInput(1) % dim = N_DIM
  diagonalizerListInput(2) % ApplyMatOnVec => ApplyA
  diagonalizerListInput(2) % ApplyPrecOnVec => ApplyPrec
  diagonalizerListInput(2) % dim = N_DIM
  ! 3, 4: AK, cluster split by nEvals, nGuard 0 / 2; 5: AK, dependent seed
  do iD = 3, 5
    diagonalizerListInput(iD) % ApplyMatOnVec => ApplyAK
    diagonalizerListInput(iD) % ApplyPrecOnVec => ApplyAKFixedPrec
    diagonalizerListInput(iD) % dim = N_DIM
  end do
  ! 6: AK with the nearly exact shift-and-invert preconditioner
  diagonalizerListInput(6) % ApplyMatOnVec => ApplyAK
  diagonalizerListInput(6) % ApplyPrecOnVec => ApplyAKShiftInvert
  diagonalizerListInput(6) % dim = N_DIM
  ! 7: slightly non-self-adjoint AK + iεW⁻¹S (stagnation)
  diagonalizerListInput(7) % ApplyMatOnVec => ApplyAKNonHermitian
  diagonalizerListInput(7) % ApplyPrecOnVec => ApplyAKFixedPrec
  diagonalizerListInput(7) % dim = N_DIM
  ! 8: the Hermitian H in the Euclidean metric, diagonal preconditioner
  diagonalizerListInput(8) % ApplyMatOnVec => ApplyH
  diagonalizerListInput(8) % ApplyPrecOnVec => ApplyPrecH
  diagonalizerListInput(8) % dim = N_DIM
  call DiagonalizerList_Fabricate(diagonalizerListInput)
  ! solvers 1-7 work in the metric with the weights W; 8 gets no
  ! metricWeights (the Euclidean metric: unit weights)
  do iD = 1, 7
    DiagonalizerList(iD) % e % metricWeights = wts
  end do
  call DiagonalizerList_Setup

  !==========================================
  call Say_Section("perform action")
  !==========================================

  print *, "cold start without preconditioner"
  call SolveAndCheck(nIterPlain, 1, ApplyA, ReferenceEvalsOfA(N_EVALS))
  if (allocated(error)) error stop "T_DiagonalizerList_Lobpcg failure (no preconditioner)"

  print *, "cold start with preconditioner"
  call SolveAndCheck(nIterPrec, 2, ApplyA, ReferenceEvalsOfA(N_EVALS))
  if (allocated(error)) error stop "T_DiagonalizerList_Lobpcg failure (preconditioner)"

  print *, "iterations without / with preconditioner: ", nIterPlain, nIterPrec
  call check(error, nIterPrec < nIterPlain / 2, "the preconditioner must reduce the iterations substantially")
  if (allocated(error)) error stop "T_DiagonalizerList_Lobpcg failure (preconditioner efficiency)"

  print *, "warm start, unchanged operator"
  call SolveAndCheck(nIterSame, 2, ApplyA, ReferenceEvalsOfA(N_EVALS))
  call check(error, nIterSame .eq. 0, "warm start from converged eigenvectors must not iterate")
  if (allocated(error)) error stop "T_DiagonalizerList_Lobpcg failure (warm start, unchanged operator)"

  print *, "warm start, perturbed operator"
  call PerturbOperator(amp=1e-4_R64)
  call SolveAndCheck(nIterWarm, 2, ApplyA, ReferenceEvalsOfA(N_EVALS))
  print *, "iterations cold / warm start: ", nIterPrec, nIterWarm
  call check(error, nIterWarm <= nIterPrec / 2, "warm start must converge in very few iterations")
  if (allocated(error)) error stop "T_DiagonalizerList_Lobpcg failure (warm start, perturbed operator)"

  print *, "Euclidean metric (no metricWeights): H with the preconditioner (diag(H) - theta)^-1"
  call SolveAndCheck(nIterEuclidean, 8, ApplyH, ReferenceEvalsOfH(N_EVALS), euclideanQ_=.true.)
  if (allocated(error)) error stop "T_DiagonalizerList_Lobpcg failure (Euclidean metric)"

  print *, "threefold cluster split by nEvals = 3, nGuard 0"
  call SolveAndCheck(nIterCluster0, 3, ApplyAK, dK(1:N_EVALS_CLUSTER))
  if (allocated(error)) error stop "T_DiagonalizerList_Lobpcg failure (cluster split by nEvals, nGuard 0)"

  print *, "threefold cluster split by nEvals = 3, nGuard 2"
  call SolveAndCheck(nIterCluster2, 4, ApplyAK, dK(1:N_EVALS_CLUSTER))
  if (allocated(error)) error stop "T_DiagonalizerList_Lobpcg failure (cluster split by nEvals, nGuard 2)"

  print *, "start block seeded with linearly dependent columns"
  allocate (startVectors(N_DIM, 6))
  do j = 1, 6
    do i = 1, N_DIM
      startVectors(i, j) = cmplx(sin(0.1_R64 * i), cos(0.03_R64 * i), kind=R64)
    end do
  end do
  startVectors(:, 5) = startVectors(:, 5) * (1.0_R64 + 1e-9_R64)
  startVectors(3, 6) = startVectors(3, 6) + 1e-7_R64
  call DiagonalizerList(5) % e % SetStartVectors(startVectors)
  call SolveAndCheck(nIterSeeded, 5, ApplyAK, dK(1:N_EVALS))
  if (allocated(error)) error stop "T_DiagonalizerList_Lobpcg failure (dependent seed)"

  print *, "nearly exact shift-and-invert preconditioner (AK - theta + 1e-12)^-1"
  call SolveAndCheck(nIterShiftInvert, 6, ApplyAK, dK(1:N_EVALS))
  if (allocated(error)) error stop "T_DiagonalizerList_Lobpcg failure (shift-and-invert preconditioner)"
  select type (e => DiagonalizerList(6) % e)
  type is (T_DiagonalizerList_E_Lobpcg)
    call check(error, e % nRawResiduals > 0, "vanishing preconditioned residuals must be replaced by raw residuals")
  end select
  if (allocated(error)) error stop "T_DiagonalizerList_Lobpcg failure (raw residual substitution)"

  print *, "slightly non-self-adjoint operator AK + i eps W^-1 S (eps = 1e-6): stagnation exit"
  call DiagonalizerList(7) % e % Diagonalize(0.0_R64, .true.)
  devStall = maxval(abs(DiagonalizerList(7) % e % evals(1:N_EVALS) - dK(1:N_EVALS)))
  select type (e => DiagonalizerList(7) % e)
  type is (T_DiagonalizerList_E_Lobpcg)
    print *, "  iterations: ", e % nIter, " of maxIter ", e % maxIter, ", matrix-vector products: ", e % nMatVec
    print *, "  convergedQ ", e % convergedQ, ", stalledQ ", e % stalledQ, ", nFound ", e % nFound
    print *, "  relative residuals: ", e % resNorms / max(1.0_R64, abs(dK(1:size(e % resNorms))))
    print *, "  max |Ritz value - eigenvalue of the self-adjoint part|: ", devStall
    call check(error, .not. e % convergedQ .and. e % stalledQ, "a residual floor must end in the stagnation exit")
    if (allocated(error)) error stop "T_DiagonalizerList_Lobpcg failure (stagnation exit)"
    call check(error, e % nIter <= MAX_ITER_STALL, "the stagnation exit must stop the solve early")
    if (allocated(error)) error stop "T_DiagonalizerList_Lobpcg failure (stagnation exit, iterations)"
    call check(error, e % nFound < N_EVALS, "nFound must count only the converged wanted columns")
    if (allocated(error)) error stop "T_DiagonalizerList_Lobpcg failure (stagnation exit, nFound)"
    call check(error, devStall < 1e-4_R64, "the best approximations must be returned")
    if (allocated(error)) error stop "T_DiagonalizerList_Lobpcg failure (stagnation exit, Ritz values)"
  class default
    error stop "T_DiagonalizerList_Lobpcg: wrong backend"
  end select

  print *
  print *, "iterations: cluster nGuard 0 / 2: ", nIterCluster0, nIterCluster2, ", dependent seed: ", nIterSeeded, &
    ", shift-and-invert: ", nIterShiftInvert

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  !> Lowest nEvals eigenvalues of A (dense reference).
  function ReferenceEvalsOfA(nEvals) result(evalsRef)
    integer(I32), intent(in) :: nEvals
    real(R64) :: evalsRef(nEvals)

    call ReferenceEvals(evalsRef, nEvals)

  end function

  !> Lowest nEvals eigenvalues of H (dense Euclidean reference).
  function ReferenceEvalsOfH(nEvals) result(evalsRef)
    integer(I32), intent(in) :: nEvals
    real(R64) :: evalsRef(nEvals)

    call ReferenceEvalsEuclidean(evalsRef, nEvals)

  end function

  !> Diagonalizes with DiagonalizerList(iD) and checks convergedQ, nFound,
  !> the eigenvalues against evalsRef, metric orthonormality and the
  !> residuals recomputed with ApplyOp (in the Euclidean metric if
  !> euclideanQ_, for a solver without metricWeights); returns the iteration
  !> count.
  subroutine SolveAndCheck(nIter, iD, ApplyOp, evalsRef, euclideanQ_)
    integer(I32), intent(out) :: nIter
    integer(I32), intent(in) :: iD
    procedure(I_ApplyMatOnVec) :: ApplyOp
    real(R64), intent(in) :: evalsRef(:)
    logical, intent(in), optional :: euclideanQ_

    real(R64) :: devEvals, devOrtho, maxRes
    complex(R64), allocatable :: av(:)
    complex(R64) :: ovlp
    integer(I32) :: i, j, nE
    logical :: euclideanQ

    euclideanQ = .false.
    if (present(euclideanQ_)) euclideanQ = euclideanQ_

    nE = int(size(evalsRef), I32)
    allocate (av(N_DIM))

    call DiagonalizerList(iD) % e % Diagonalize(0.0_R64, .true.)

    select type (e => DiagonalizerList(iD) % e)
    type is (T_DiagonalizerList_E_Lobpcg)
      nIter = e % nIter
      print *, "  iterations: ", e % nIter, ", matrix-vector products: ", e % nMatVec, &
        ", raw residuals substituted: ", e % nRawResiduals
    class default
      error stop "T_DiagonalizerList_Lobpcg: wrong backend"
    end select

    associate (evals => DiagonalizerList(iD) % e % evals, evecs => DiagonalizerList(iD) % e % evecs)

      call check(error, DiagonalizerList(iD) % e % convergedQ .and. DiagonalizerList(iD) % e % nFound .eq. nE, &
                 "the solve must converge with nFound = nEvals")
      if (allocated(error)) return

      ! eigenvalues against the reference
      devEvals = maxval(abs(evals(1:nE) - evalsRef(:)))
      print *, "  evals:     ", evals(1:nE)
      print *, "  reference: ", evalsRef(:)
      print *, "  max eigenvalue deviation: ", devEvals
      call check(error, devEvals < 1e-10_R64, "eigenvalues must match the reference")
      if (allocated(error)) return

      ! orthonormality in the metric
      devOrtho = 0.0_R64
      do j = 1, nE
        do i = 1, nE
          ovlp = InnerProductOfSolve(evecs(:, i), evecs(:, j), euclideanQ)
          if (i .eq. j) ovlp = ovlp - 1.0_R64
          devOrtho = max(devOrtho, abs(ovlp))
        end do
      end do
      print *, "  max metric orthonormality deviation: ", devOrtho
      call check(error, devOrtho < 1e-12_R64, "eigenvectors must be orthonormal in the metric")
      if (allocated(error)) return

      ! residuals recomputed from scratch
      maxRes = 0.0_R64
      do j = 1, nE
        call ApplyOp(av, evecs(:, j), 0.0_R64)
        av(:) = av(:) - evals(j) * evecs(:, j)
        maxRes = max(maxRes, sqrt(real(InnerProductOfSolve(av, av, euclideanQ), kind=R64)) &
                     / max(1.0_R64, abs(evals(j))))
      end do
      print *, "  max recomputed relative residual: ", maxRes
      call check(error, maxRes < 1e-9_R64, "recomputed residuals must be small")

    end associate

  end subroutine

end program
