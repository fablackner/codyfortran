! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> GSL ODEIV2 interface: stepper, control, and evolve wrappers.
!>
!> Exposes a small Fortran-friendly layer over the C API to integrate ODE
!> systems via user-supplied right-hand sides and method selection.
module M_Utils_Odeiv2GslLib
  use M_Utils_Types
  use, intrinsic :: iso_c_binding
  implicit none

  private :: GSL_SUCCESS, GSL_FAILURE
  integer(c_int), parameter :: GSL_SUCCESS = 0
  integer(c_int), parameter :: GSL_FAILURE = -1

  abstract interface
    subroutine I_TimeDerivative(dState, state, time)
      import :: R64
      complex(R64), intent(out), contiguous, target :: dState(:)
      complex(R64), intent(in), contiguous, target :: state(:)
      real(R64), intent(in) :: time
    end subroutine
  end interface

  type, bind(C) :: gsl_odeiv2_Method
    type(c_funptr)     :: function
    type(c_funptr)     :: jacobian
    integer(c_size_t)  :: dimension
    type(c_ptr)        :: params
  end type

  type :: T_Odeiv2GslLib_Ctx
    type(c_ptr) :: driver = c_null_ptr
    type(gsl_odeiv2_Method) :: system
    procedure(I_TimeDerivative), pointer, nopass :: TimeDerivativeSaved => null()
    integer(I32) :: dimension = 0
  end type

  public :: T_Odeiv2GslLib_Ctx
  public :: Odeiv2GslLib_CreateCtx
  public :: Odeiv2GslLib_Integrate
  public :: Odeiv2GslLib_DestroyCtx

  interface
    function gsl_odeiv2_driver_alloc_y_new(sys, stepType, stepSize, absError, relError) bind(C) result(driver)
      import :: c_ptr, c_double
      type(c_ptr) :: driver
      type(c_ptr), value :: sys
      type(c_ptr), value :: stepType
      real(c_double), value :: stepSize
      real(c_double), value :: absError
      real(c_double), value :: relError
    end function
    function gsl_odeiv2_driver_apply_fixed_step(driver, time, stepSize, steps, state) bind(C) result(status)
      import :: c_ptr, c_double, c_size_t, c_int
      integer(c_int) :: status
      type(c_ptr), value :: driver
      type(c_ptr), value :: time
      real(c_double), value :: stepSize
      integer(c_size_t), value :: steps
      type(c_ptr), value :: state
    end function
    subroutine gsl_odeiv2_driver_free(driver) bind(C)
      import :: c_ptr
      type(c_ptr), value :: driver
    end subroutine
  end interface

  type(c_ptr), bind(C) :: gsl_odeiv2_step_rk2
  type(c_ptr), bind(C) :: gsl_odeiv2_step_rk4
  type(c_ptr), bind(C) :: gsl_odeiv2_step_rkf45
  type(c_ptr), bind(C) :: gsl_odeiv2_step_rkck
  type(c_ptr), bind(C) :: gsl_odeiv2_step_rk8pd

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function GslDerivativeWrapper(t, yPtr, dydtPtr, params) bind(C) result(status)
    real(c_double), value :: t
    type(c_ptr), value :: yPtr, dydtPtr, params
    integer(c_int) :: status
    type(T_Odeiv2GslLib_Ctx), pointer :: ctx
    real(R64), pointer :: yReal(:), dydtReal(:)
    complex(R64), allocatable :: yComplex(:), dydtComplex(:)
    integer :: n, i

    call c_f_pointer(params, ctx)
    if (.not. associated(ctx)) then
      status = GSL_FAILURE
      return
    end if
    n = ctx % dimension
    call c_f_pointer(yPtr, yReal, [2 * n])
    call c_f_pointer(dydtPtr, dydtReal, [2 * n])

    allocate (yComplex(n), dydtComplex(n))
    do i = 1, n
      yComplex(i) = cmplx(yReal(i), yReal(n + i), kind=R64)
    end do

    if (associated(ctx % TimeDerivativeSaved)) then
      call ctx % TimeDerivativeSaved(dydtComplex, yComplex, t)
      status = GSL_SUCCESS
      do i = 1, n
        dydtReal(i) = real(dydtComplex(i), kind=R64)
        dydtReal(n + i) = aimag(dydtComplex(i))
      end do
    else
      status = GSL_FAILURE
    end if

    deallocate (yComplex, dydtComplex)
  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Odeiv2GslLib_CreateCtx(ctx, dimension, TimeDerivative, stepSize, absError, relError, stepperType)
    type(T_Odeiv2GslLib_Ctx), intent(out), target :: ctx
    integer(I32), intent(in) :: dimension
    procedure(I_TimeDerivative) :: TimeDerivative
    real(R64), intent(in) :: stepSize, absError, relError
    character(len=*), intent(in) :: stepperType
    type(c_ptr) :: selectedStepper
    real(c_double) :: h, epsAbs, epsRel

    ctx % TimeDerivativeSaved => TimeDerivative
    ctx % dimension = dimension

    select case (trim(adjustl(stepperType)))
    case ("rk2"); selectedStepper = gsl_odeiv2_step_rk2
    case ("rk4"); selectedStepper = gsl_odeiv2_step_rk4
    case ("rkf45"); selectedStepper = gsl_odeiv2_step_rkf45
    case ("rkck"); selectedStepper = gsl_odeiv2_step_rkck
    case ("rk8pd"); selectedStepper = gsl_odeiv2_step_rk8pd
    case default
      error stop "Odeiv2GslLib_CreateCtx: unknown stepper "//trim(stepperType)
    end select

    ctx % system % function = c_funloc(GslDerivativeWrapper)
    ctx % system % jacobian = c_null_funptr
    ctx % system % dimension = int(2 * dimension, c_size_t)
    ctx % system % params = c_loc(ctx)

    h = stepSize
    epsAbs = absError
    epsRel = relError
    ctx % driver = gsl_odeiv2_driver_alloc_y_new(c_loc(ctx % system), selectedStepper, h, epsAbs, epsRel)
    if (.not. c_associated(ctx % driver)) error stop "Odeiv2GslLib_CreateCtx: driver allocation failed"
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Odeiv2GslLib_DestroyCtx(ctx)
    type(T_Odeiv2GslLib_Ctx), intent(inout) :: ctx
    if (c_associated(ctx % driver)) call gsl_odeiv2_driver_free(ctx % driver)
    ctx % driver = c_null_ptr
    ctx % TimeDerivativeSaved => null()
    ctx % dimension = 0
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Odeiv2GslLib_Integrate(status, y, t0, t1, nSteps, ctx)
    integer(I32), intent(out) :: status
    complex(R64), intent(inout), target :: y(:)
    real(R64), intent(in) :: t0, t1
    integer, intent(inout) :: nSteps
    type(T_Odeiv2GslLib_Ctx), intent(in) :: ctx

    real(R64) :: dt
    integer :: n, i
    real(R64), allocatable, target :: yReal(:)
    real(c_double), target :: t_c
    integer(c_int) :: gslStatus

    if (.not. c_associated(ctx % driver)) then
      status = GSL_FAILURE
      return
    end if

    n = size(y)
    if (n .ne. ctx % dimension) error stop "Odeiv2GslLib_Integrate: dimension mismatch"

    if (nSteps <= 0) nSteps = 1
    dt = (t1 - t0) / real(nSteps, R64)

    allocate (yReal(2 * n))
    do i = 1, n
      yReal(i) = real(y(i), kind=R64)
      yReal(n + i) = aimag(y(i))
    end do

    t_c = t0
    gslStatus = gsl_odeiv2_driver_apply_fixed_step(ctx % driver, c_loc(t_c), dt, int(nSteps, c_size_t), c_loc(yReal))

    do i = 1, n
      y(i) = cmplx(yReal(i), yReal(n + i), kind=R64)
    end do

    deallocate (yReal)
    status = int(gslStatus, I32)
  end subroutine

end module
