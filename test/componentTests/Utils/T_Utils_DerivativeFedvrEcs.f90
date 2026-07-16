program T_Utils_DerivativeFedvrEcs
  use M_Utils_Types
  use M_Utils_DerivativeFedvrEcs
  use M_Utils_FedvrEcs
  use testdrive, only: check, error_type
  implicit none

  type(error_type), allocatable :: error

  ! Test contour quadrature and both derivative operators on an ECS grid
  call test_contour_operators(error)
  if (allocated(error)) error stop 'DerivativeFedvrEcs contour operator test failed'

contains

  subroutine test_contour_operators(error)
    type(error_type), allocatable, intent(out) :: error

    type(T_FedvrEcs_Ctx) :: fedvrEcsCtx
    type(T_DerivativeFedvrEcs_Ctx) :: derivCtx
    complex(R64), allocatable :: f(:), dfDx(:), d2fDx2(:), expected1(:), expected2(:)
    complex(R64) :: z, zEnd, integral, exactIntegral
    real(R64) :: relErr1, relErr2, quadErr, xmin, xmax, theta
    real(R64) :: lambda, pi
    integer(I32) :: i, nElements, nLocals, nPoints

    print *, "Testing DerivativeFedvrEcs on an ECS contour..."

    ! Set parameters for the grid; ECS starts at r = 10 with a 2-element ramp
    nElements = 20
    nLocals = 15
    xmin = 0.0_R64
    xmax = 20.0_R64
    lambda = 2.0_R64
    pi = acos(-1.0_R64)
    theta = pi / 6.0_R64

    ! Create FEDVR-ECS grid; exclude r = 0 (test function vanishes there)
    call FedvrEcs_CreateCtx(fedvrEcsCtx, xmin, xmax, nElements, nLocals, .true., .false., &
                            ecsStartElement_=11, ecsAngle_=theta, ecsTransitionElements_=2)
    nPoints = fedvrEcsCtx % nPoints

    ! Create derivative context
    call DerivativeFedvrEcs_CreateCtx(derivCtx, fedvrEcsCtx)

    allocate (f(nPoints), dfDx(nPoints), d2fDx2(nPoints), expected1(nPoints), expected2(nPoints))

    ! Analytic test function f(z) = z * exp(-lambda z) continued onto the contour
    do i = 1, nPoints
      z = fedvrEcsCtx % points(i)
      f(i) = z * exp(-lambda * z)
      expected1(i) = (1.0_R64 - lambda * z) * exp(-lambda * z)
      expected2(i) = (lambda**2 * z - 2.0_R64 * lambda) * exp(-lambda * z)
    end do

    ! Contour integral of an analytic, decaying function equals the real-axis one
    zEnd = fedvrEcsCtx % points(nPoints)
    integral = sum(f * fedvrEcsCtx % weights)
    exactIntegral = (1.0_R64 - exp(-lambda * zEnd) * (1.0_R64 + lambda * zEnd)) / lambda**2
    quadErr = abs(integral - exactIntegral)
    print *, "Contour quadrature error: ", quadErr
    call check(error, quadErr < 1.0e-12_R64)
    if (allocated(error)) return

    ! First derivative along the contour
    call DerivativeFedvrEcs_Do1stDerivative(dfDx, f, derivCtx, fedvrEcsCtx)
    relErr1 = maxval(abs(dfDx - expected1)) / maxval(abs(expected1))
    print *, "Relative error of 1st derivative: ", relErr1
    call check(error, relErr1 < 1.0e-8_R64)
    if (allocated(error)) return

    ! Second derivative along the contour (skip global endpoints)
    call DerivativeFedvrEcs_Do2ndDerivative(d2fDx2, f, derivCtx, fedvrEcsCtx)
    relErr2 = maxval(abs(d2fDx2(2:nPoints - 1) - expected2(2:nPoints - 1))) / maxval(abs(expected2))
    print *, "Relative error of 2nd derivative: ", relErr2
    call check(error, relErr2 < 1.0e-8_R64)
    if (allocated(error)) return

    ! Clean up
    deallocate (f, dfDx, d2fDx2, expected1, expected2)
    call DerivativeFedvrEcs_DestroyCtx(derivCtx)
    call FedvrEcs_DestroyCtx(fedvrEcsCtx)
  end subroutine test_contour_operators

end program T_Utils_DerivativeFedvrEcs
