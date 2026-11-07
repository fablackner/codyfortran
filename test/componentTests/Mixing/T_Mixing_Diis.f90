!> Mixing module (mixing.diis, mixing.linear) on affine fixed-point maps
!> x ← A x + b with a symmetric A of known spectrum:
!> - Anderson finite termination: DIIS with α = 1 and a history longer than
!>   the dimension finds the fixed point of an n-dimensional affine map within
!>   n + 1 steps, also when the plain iteration diverges (|λ| > 1: here
!>   λ = −1.4 and 1.2) or oscillates (λ = −0.95)
!> - a history of one pair is linear mixing: DIIS with nHistory 1 gives the
!>   iterates of mixing.linear with the same α to round-off
!> - Mixing_Reset empties the history: the next step is the plain damped
!>   step (with an unconverged history of another map stored, which the
!>   growth safeguard does not drop), and the new map terminates within
!>   n + 1 steps
!> - a new vector length resets the mixer (no reset call needed)
program T_Mixing_Diis
  use M_Utils_Types
  use M_Utils_Json
  use M_Mixing
  use testdrive, only: check, error_type
  implicit none

  integer(I32), parameter :: n = 6
  real(R64), parameter :: lambdas(n) = [-1.4_R64, -0.95_R64, 0.3_R64, 0.7_R64, 0.9_R64, 1.2_R64]

  type(error_type), allocatable :: error
  real(R64) :: a(n, n)
  complex(R64) :: b(n), b2(n), xStar(n), xStar2(n)
  complex(R64), allocatable :: xDiis(:), xLin(:), xRaw(:)
  real(R64) :: err, maxDiff
  integer(I32) :: nSteps, k
  logical :: foundQ

  call Json_LoadJsonFile(manualFileName_="test/componentTests/Mixing/T_Mixing_Diis.json", relativeToProjectDirQ_=.true.)

  call BuildMap(a, lambdas)
  b = [(cmplx(1.0_R64 / real(k, R64), 0.5_R64 * real(k, R64), kind=R64), k=1, n)]
  b2 = [(cmplx(-0.5_R64 / real(k, R64), 0.4_R64 * real(n + 1 - k, R64), kind=R64), k=1, n)]
  xStar = FixedPoint(a, b)
  xStar2 = FixedPoint(a, b2)

  !------------------------------------
  ! Anderson finite termination (startThreshold 1e30: no start gate)
  !------------------------------------

  call Mixing_Fabricate
  allocate (xDiis(n), xRaw(n))
  xDiis = (0.0_R64, 0.0_R64)
  call Iterate(nSteps, err, xDiis, a, b, xStar, 3 * n)
  write (*, '(1X, A, I0, A, ES10.2)') "diis, map 1: ", nSteps, " steps, error ", err
  call check(error, nSteps <= n + 1)
  if (allocated(error)) error stop "T_Mixing_Diis failure: Anderson termination within n + 1 steps"

  !------------------------------------
  ! Mixing_Reset before another map: three steps of map 1 (history stored,
  ! residuals of the size of those of map 2), reset, map 2
  !------------------------------------

  call Mixing_Reset
  xDiis = (0.0_R64, 0.0_R64)
  do k = 1, 3
    xRaw = matmul(a, xDiis) + b
    call Mixing_Mix(xDiis, xRaw)
  end do
  call Mixing_Reset
  xRaw = matmul(a, xDiis) + b2
  call Mixing_Mix(xDiis, xRaw)
  maxDiff = maxval(abs(xDiis - xRaw))
  write (*, '(1X, A, ES10.2)') "diis, first step after Mixing_Reset - plain step (alpha 1): ", maxDiff
  call check(error, maxDiff < 1e-14_R64 * maxval(abs(xRaw)))
  if (allocated(error)) error stop "T_Mixing_Diis failure: the step after Mixing_Reset is not the plain step"
  call Iterate(nSteps, err, xDiis, a, b2, xStar2, 3 * n)
  write (*, '(1X, A, I0, A, ES10.2)') "diis, map 2 after Mixing_Reset: ", nSteps + 1, " steps, error ", err
  call check(error, nSteps + 1 <= n + 1)
  if (allocated(error)) error stop "T_Mixing_Diis failure: termination after Mixing_Reset"

  !------------------------------------
  ! a new vector length resets the mixer: the leading n - 1 components of
  ! the map with the last row and column removed
  !------------------------------------

  deallocate (xDiis)
  allocate (xDiis(n - 1))
  xDiis = (0.0_R64, 0.0_R64)
  call Iterate(nSteps, err, xDiis, a(1:n - 1, 1:n - 1), b(1:n - 1), &
               FixedPoint(a(1:n - 1, 1:n - 1), b(1:n - 1)), 3 * n)
  write (*, '(1X, A, I0, A, ES10.2)') "diis, length n - 1: ", nSteps, " steps, error ", err
  call check(error, nSteps <= n)
  if (allocated(error)) error stop "T_Mixing_Diis failure: termination after a change of the vector length"

  !------------------------------------
  ! one stored pair = linear mixing with the same alpha (a contracting map:
  ! the spectrum scaled into (-1, 1))
  !------------------------------------

  call BuildMap(a, lambdas / 1.5_R64)
  call jsonFile % update("mixing.diis.nHistory", 1, foundQ)
  if (.not. foundQ) error stop "T_Mixing_Diis: mixing.diis.nHistory not found"
  call jsonFile % update("mixing.diis.alpha", 0.6_R64, foundQ)
  if (.not. foundQ) error stop "T_Mixing_Diis: mixing.diis.alpha not found"
  call Mixing_Fabricate
  deallocate (xDiis)
  allocate (xDiis(n), xLin(n))
  xDiis = (0.0_R64, 0.0_R64)
  do k = 1, 10
    xRaw = matmul(a, xDiis) + b
    call Mixing_Mix(xDiis, xRaw)
  end do

  call jsonFile % remove("mixing.diis")
  call jsonFile % add("mixing.linear.alpha", 0.6_R64)
  call Mixing_Fabricate
  xLin = (0.0_R64, 0.0_R64)
  do k = 1, 10
    xRaw = matmul(a, xLin) + b
    call Mixing_Mix(xLin, xRaw)
  end do
  maxDiff = maxval(abs(xDiis - xLin))
  write (*, '(1X, A, ES10.2)') "diis with nHistory 1 vs linear, 10 steps: max difference ", maxDiff
  call check(error, maxDiff < 1e-13_R64 * maxval(abs(xLin)))
  if (allocated(error)) error stop "T_Mixing_Diis failure: nHistory 1 differs from linear mixing"

contains

  !> A = Q diag(lambda) Qᵀ with the Householder reflection Q = 1 − 2 v vᵀ of a
  !> fixed unit vector v (dense, symmetric, spectrum lambda).
  subroutine BuildMap(a, lambda)
    real(R64), intent(out) :: a(:, :)
    real(R64), intent(in) :: lambda(:)

    real(R64) :: v(size(lambda)), q(size(lambda), size(lambda))
    integer(I32) :: i, j, m

    m = size(lambda)
    v = [(real(i, R64) + 0.5_R64 * real(mod(i, 3), R64), i=1, m)]
    v = v / norm2(v)
    do j = 1, m
      do i = 1, m
        q(i, j) = -2.0_R64 * v(i) * v(j)
      end do
      q(j, j) = q(j, j) + 1.0_R64
    end do
    do j = 1, m
      do i = 1, m
        a(i, j) = sum(q(i, :) * lambda * q(j, :))
      end do
    end do

  end subroutine

  !> x* = (1 − A)⁻¹ b by Gaussian elimination with partial pivoting.
  function FixedPoint(a, b) result(x)
    real(R64), intent(in) :: a(:, :)
    complex(R64), intent(in) :: b(:)
    complex(R64) :: x(size(b))

    real(R64) :: m(size(b), size(b)), factor
    real(R64), allocatable :: rowTmp(:)
    complex(R64) :: tmp
    integer(I32) :: i, k, p, nn

    nn = size(b)
    m = -a
    do i = 1, nn
      m(i, i) = m(i, i) + 1.0_R64
    end do
    x = b
    do k = 1, nn
      p = k - 1 + maxloc(abs(m(k:, k)), dim=1)
      if (p .ne. k) then
        rowTmp = m(k, :)
        m(k, :) = m(p, :)
        m(p, :) = rowTmp
        tmp = x(k)
        x(k) = x(p)
        x(p) = tmp
      end if
      do i = k + 1, nn
        factor = m(i, k) / m(k, k)
        m(i, k:) = m(i, k:) - factor * m(k, k:)
        x(i) = x(i) - factor * x(k)
      end do
    end do
    do k = nn, 1, -1
      x(k) = (x(k) - sum(m(k, k + 1:) * x(k + 1:))) / m(k, k)
    end do

  end function

  !> Mixing_Mix steps x ← mix(x, A x + b) until |x − x*| < 1e-9 |x*| (at
  !> most maxSteps): nSteps the number of Mixing_Mix calls, err the final
  !> relative error.
  subroutine Iterate(nSteps, err, x, a, b, xStar, maxSteps)
    integer(I32), intent(out) :: nSteps
    real(R64), intent(out) :: err
    complex(R64), intent(inout), contiguous :: x(:)
    real(R64), intent(in) :: a(:, :)
    complex(R64), intent(in) :: b(:)
    complex(R64), intent(in) :: xStar(:)
    integer(I32), intent(in) :: maxSteps

    complex(R64), allocatable :: xNew(:)

    do nSteps = 1, maxSteps
      xNew = matmul(a, x) + b
      call Mixing_Mix(x, xNew)
      err = norm2(abs(x - xStar)) / norm2(abs(xStar))
      if (err < 1e-9_R64) return
    end do
    nSteps = maxSteps + 1

  end subroutine

end program
