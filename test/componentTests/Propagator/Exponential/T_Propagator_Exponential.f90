!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Verifies the exponential integrators (propagator.exponential) on a small
!> model u' = L u + N1(t, u) + N2(u):
!> - L = -i diag(omega), propagated exactly by test bindings of
!>   Method_PropagateExactPart and Method_ApplyExactPartPhi
!> - N1 = -i A(t) M u, a time-dependent Hermitian coupling (IntegratorList(1))
!> - N2 = -i g |u|^2 u, a norm-conserving nonlinearity (IntegratorList(2))
!>
!> For every scheme (krogstad4, adams and adamsLawson with nSteps = 4) the
!> error at t = 2 against a fine classical RK4 solution must decrease with
!> the order of the scheme when the step is halved (4 for Krogstad, 5 for the
!> k = 4 PECE schemes), and the norm must be conserved to the same accuracy.
program T_Propagator_Exponential
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Method
  use M_IntegratorList
  use M_Propagator
  use M_Propagator_Exponential_Adams
  use M_Propagator_Exponential_AdamsLawson
  use testdrive, only: check, error_type

  implicit none

  integer(I32), parameter :: nDim = 4
  real(R64), parameter :: omega(nDim) = [0.0_R64, 1.3_R64, 2.9_R64, 4.1_R64]
  real(R64), parameter :: endTime = 2.0_R64

  type(error_type), allocatable :: error
  type(T_IntegratorList_FabricateInput) :: IntegratorListInput(2)
  character(len=:), allocatable :: jsonFileName
  complex(R64) :: reference(nDim)
  real(R64) :: errCoarse, errFine, normErr

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  jsonFileName = "test/componentTests/Propagator/Exponential/T_Propagator_Exponential.json"
  call Json_LoadJsonFile(manualFileName_=jsonFileName, relativeToProjectDirQ_=.true.)

  Method_PropagateExactPart => PropagateExactPart
  Method_ApplyExactPartPhi => ApplyExactPartPhi

  IntegratorListInput(1) % TimeDerivative => Coupling
  IntegratorListInput(2) % TimeDerivative => Nonlinearity
  call IntegratorList_Fabricate(IntegratorListInput)
  call Propagator_Fabricate()

  reference = ReferenceSolution(20000)

  !==========================================
  call Say_Section("krogstad4")
  !==========================================

  call Errors(errCoarse, errFine, normErr, 50)
  write (*, '(A,3ES12.3)') " krogstad4:   ", errCoarse, errFine, errCoarse / errFine
  call check(error, errCoarse / errFine > 12.0_R64)
  if (allocated(error)) error stop "T_Propagator_Exponential failure: krogstad4 is not of order 4"
  call check(error, errFine < 1e-6_R64 .and. normErr < 1e-6_R64)
  if (allocated(error)) error stop "T_Propagator_Exponential failure: krogstad4 error too large"

  !==========================================
  call Say_Section("adams")
  !==========================================

  call Propagator_Exponential_Adams_Fabricate
  call Errors(errCoarse, errFine, normErr, 50)
  write (*, '(A,3ES12.3)') " adams:       ", errCoarse, errFine, errCoarse / errFine
  call check(error, errCoarse / errFine > 24.0_R64)
  if (allocated(error)) error stop "T_Propagator_Exponential failure: adams is not of order 5"
  call check(error, errFine < 1e-6_R64 .and. normErr < 1e-6_R64)
  if (allocated(error)) error stop "T_Propagator_Exponential failure: adams error too large"

  !==========================================
  call Say_Section("adamsLawson")
  !==========================================

  call Propagator_Exponential_AdamsLawson_Fabricate
  call Errors(errCoarse, errFine, normErr, 50)
  write (*, '(A,3ES12.3)') " adamsLawson: ", errCoarse, errFine, errCoarse / errFine
  call check(error, errCoarse / errFine > 24.0_R64)
  if (allocated(error)) error stop "T_Propagator_Exponential failure: adamsLawson is not of order 5"
  call check(error, errFine < 1e-6_R64 .and. normErr < 1e-6_R64)
  if (allocated(error)) error stop "T_Propagator_Exponential failure: adamsLawson error too large"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  !> Errors at endTime with nSteps and 2 nSteps steps, and the norm error of the latter.
  subroutine Errors(errCoarse, errFine, normErr, nSteps)
    real(R64), intent(out) :: errCoarse, errFine, normErr
    integer(I32), intent(in) :: nSteps

    complex(R64) :: u(nDim)

    u = Propagate(nSteps)
    errCoarse = norm2(abs(u - reference))
    u = Propagate(2 * nSteps)
    errFine = norm2(abs(u - reference))
    normErr = abs(norm2(abs(u)) - norm2(abs(InitialState())))

  end subroutine

  function Propagate(nSteps) result(u)
    integer(I32), intent(in) :: nSteps
    complex(R64) :: u(nDim)

    integer(I32) :: i
    real(R64) :: h

    h = endTime / nSteps
    u = InitialState()
    do i = 0, nSteps - 1
      call Propagator_Propagate(u, i * h, (i + 1) * h)
    end do

  end function

  !> Classical RK4 of the full right-hand side.
  function ReferenceSolution(nSteps) result(u)
    integer(I32), intent(in) :: nSteps
    complex(R64) :: u(nDim)

    complex(R64), dimension(nDim) :: k1, k2, k3, k4
    integer(I32) :: i
    real(R64) :: h, t

    h = endTime / nSteps
    u = InitialState()
    do i = 0, nSteps - 1
      t = i * h
      k1 = FullDerivative(u, t)
      k2 = FullDerivative(u + 0.5_R64 * h * k1, t + 0.5_R64 * h)
      k3 = FullDerivative(u + 0.5_R64 * h * k2, t + 0.5_R64 * h)
      k4 = FullDerivative(u + h * k3, t + h)
      u = u + (h / 6.0_R64) * (k1 + 2.0_R64 * k2 + 2.0_R64 * k3 + k4)
    end do

  end function

  function FullDerivative(u, time) result(du)
    complex(R64), intent(in) :: u(nDim)
    real(R64), intent(in) :: time
    complex(R64) :: du(nDim), d1(nDim), d2(nDim)

    call Coupling(d1, u, time)
    call Nonlinearity(d2, u, time)
    du = cmplx(0.0_R64, -omega, kind=R64) * u + d1 + d2

  end function

  function InitialState() result(u)
    complex(R64) :: u(nDim)

    u = [(0.8_R64, 0.0_R64), (0.0_R64, 0.5_R64), (0.3_R64, 0.1_R64), (-0.1_R64, 0.2_R64)]

  end function

  subroutine Coupling(dState, state, time)
    complex(R64), intent(out), contiguous, target :: dState(:)
    complex(R64), intent(in), contiguous, target :: state(:)
    real(R64), intent(in) :: time

    real(R64) :: a
    integer(I32) :: k

    ! nearest-neighbour Hermitian coupling with a smooth field A(t)
    a = 0.6_R64 * sin(1.7_R64 * time)
    dState(:) = (0.0_R64, 0.0_R64)
    do k = 1, nDim - 1
      dState(k) = dState(k) + state(k + 1)
      dState(k + 1) = dState(k + 1) + state(k)
    end do
    dState(:) = cmplx(0.0_R64, -a, kind=R64) * dState(:)

  end subroutine

  subroutine Nonlinearity(dState, state, time)
    complex(R64), intent(out), contiguous, target :: dState(:)
    complex(R64), intent(in), contiguous, target :: state(:)
    real(R64), intent(in) :: time

    dState(:) = cmplx(0.0_R64, -0.4_R64, kind=R64) * abs(state(:))**2 * state(:)

  end subroutine

  subroutine PropagateExactPart(state, t0, t1, imagTimeQ_)
    complex(R64), intent(inout), contiguous, target :: state(:)
    real(R64), intent(in) :: t0
    real(R64), intent(in) :: t1
    logical, intent(in), optional :: imagTimeQ_

    state(:) = exp(cmplx(0.0_R64, -omega * (t1 - t0), kind=R64)) * state(:)

  end subroutine

  subroutine ApplyExactPartPhi(res, vecs, dt)
    complex(R64), intent(out), contiguous, target :: res(:)
    complex(R64), intent(in), contiguous, target :: vecs(:, :)
    real(R64), intent(in) :: dt

    complex(R64) :: phi(0:size(vecs, 2) - 1)
    integer(I32) :: k

    do k = 1, nDim
      phi = PhiFunctions(cmplx(0.0_R64, -omega(k) * dt, kind=R64), size(vecs, 2) - 1)
      res(k) = sum(phi(:) * vecs(k, :))
    end do

  end subroutine

  !> phi_0(z) ... phi_K(z) by their Taylor series (|z| <= 0.4 here).
  function PhiFunctions(z, nMax) result(res)
    complex(R64), intent(in) :: z
    integer(I32), intent(in) :: nMax
    complex(R64) :: res(0:nMax)

    complex(R64) :: term
    integer(I32) :: m, j

    do m = 0, nMax
      term = cmplx(1.0_R64 / gamma(real(m + 1, kind=R64)), 0.0_R64, kind=R64)
      res(m) = term
      do j = 1, 30
        term = term * z / real(j + m, kind=R64)
        res(m) = res(m) + term
      end do
    end do

  end function

end program
