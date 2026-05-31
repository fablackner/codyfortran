program T_Utils_DerivativeFedvr
  use M_Utils_Types
  use M_Utils_DerivativeFedvr
  use M_Utils_Fedvr
  use testdrive, only: check, error_type, new_unittest
  implicit none

  type(error_type), allocatable :: error

  ! Test the second derivative operator with different boundary conditions
  call test_second_derivative(error)
  if (allocated(error)) error stop 'DerivativeFedvr second derivative test failed'

contains

  subroutine test_second_derivative(error)
    type(error_type), allocatable, intent(out) :: error

    type(T_Fedvr_Ctx) :: fedvrCtx
    type(T_DerivativeFedvr_Ctx) :: derivCtx
    complex(R64), allocatable :: f(:), d2fDx2(:), expected(:)
    real(R64) :: maxErr, relErr, r, xmin, xmax
    real(R64) :: lambda
    integer(I32) :: i, nElements, nLocals
    ! Parameters for test wavefunction f(r) = r * exp(-lambda * r)

    print *, "Testing DerivativeFedvr second derivative..."

    ! Set parameters for the grid
    nElements = 10
    nLocals = 21
    xmin = 0.0_R64
    xmax = 20.0_R64
    lambda = 2.0_R64

    ! Create FEDVR grid with included endpoints
    call Fedvr_CreateCtx(fedvrCtx, xmin, xmax, nElements, nLocals, .true., .false.)

    ! Create derivative context
    call DerivativeFedvr_CreateCtx(derivCtx, fedvrCtx)

    ! Allocate arrays for function, computed second derivative, and expected result
    allocate (f(fedvrCtx % nPoints), d2fDx2(fedvrCtx % nPoints), expected(fedvrCtx % nPoints))

    ! Initialize test wavefunction f(r) = r * exp(-lambda * r) and analytic second derivative
    do i = 1, fedvrCtx % nPoints
      r = fedvrCtx % points(i)
      f(i) = r * exp(-lambda * r)
      expected(i) = ((-2.0_R64 * lambda) + lambda**2 * r) * exp(-lambda * r)
    end do

    ! Compute second derivative using DerivativeFedvr
    call DerivativeFedvr_Do2ndDerivative(d2fDx2, f, derivCtx, fedvrCtx)

    ! Calculate maximum relative error
    maxErr = 0.0_R64
    do i = 1, fedvrCtx % nPoints
      maxErr = max(maxErr, abs(d2fDx2(i) - expected(i)))
    end do
    relErr = maxErr / maxval(abs(expected))

    print *, "Relative error with included endpoints: ", relErr

    ! Check if error is acceptable
    call check(error, relErr < 1.0e-8_R64)

    ! Clean up
    deallocate (f, d2fDx2, expected)
    call DerivativeFedvr_DestroyCtx(derivCtx)
    call Fedvr_DestroyCtx(fedvrCtx)
  end subroutine test_second_derivative

end program T_Utils_DerivativeFedvr
