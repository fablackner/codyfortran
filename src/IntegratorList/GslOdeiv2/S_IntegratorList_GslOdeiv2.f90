! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file S_IntegratorList_GslOdeiv2.f90
!> @brief GSL odeiv2 adaptive stepper integrator implementation.
!>
!> @details
!> Wraps the GNU Scientific Library (GSL) odeiv2 ODE solver framework. GSL
!> provides a suite of adaptive-step integrators with automatic error control.
!>
!> ## Supported Steppers
!>
!> | Type     | Method                                     | Order |
!> |----------|--------------------------------------------|-------|
!> | `rk4`    | Classical Runge–Kutta (fixed step)         | 4     |
!> | `rkf45`  | Runge–Kutta–Fehlberg (adaptive)            | 4(5)  |
!> | `rkck`   | Runge–Kutta–Cash–Karp (adaptive)           | 4(5)  |
!> | `rk8pd`  | Runge–Kutta Prince–Dormand (adaptive)      | 8(9)  |
!>
!> ## Configuration (JSON)
!>
!> | Key           | Type   | Default | Description                |
!> |---------------|--------|---------|----------------------------|
!> | `stepperType` | string | `"rk4"` | GSL stepper algorithm name |
!>
!> ## Error Tolerances
!>
!> Default absolute/relative errors are set internally:
!> - Absolute: 1.0e-6
!> - Relative: 1.0e-10
!>
!> @note The GSL context is lazily initialized on the first call to
!>   `Integrate`, at which point the state dimension becomes known.
submodule(M_IntegratorList_GslOdeiv2) S_IntegratorList_GslOdeiv2

  implicit none

  ! Default error tolerances for GSL adaptive steppers
  real(R64), parameter :: DEFAULT_ABS_ERROR = 1.0e-6_R64
  real(R64), parameter :: DEFAULT_REL_ERROR = 1.0e-10_R64

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine IntegratorList_GslOdeiv2_Allocate(e, path)
    use M_Utils_UnusedVariables
    class(T_IntegratorList_E), allocatable, intent(out) :: e
    character(len=*), intent(in) :: path

    if (.false.) call UnusedVariables_Mark(path)

    allocate (T_IntegratorList_E_GslOdeiv2 :: e)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Fabricate(this)
    use M_Utils_Json
    use M_Utils_Say

    class(T_IntegratorList_E_GslOdeiv2), intent(inout) :: this

    call Say_Fabricate(this % path//".gslOdeiv2")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    this % stepperType = Json_Get("stepperType", "rk4", path_=this % path)
    this % initializedQ = .false.  ! Will be initialized on first call to Integrate

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Setup(this)
    use M_Utils_Say
    use M_Utils_UnusedVariables

    class(T_IntegratorList_E_GslOdeiv2), intent(inout) :: this

    call Say_Setup(this % path//".gslOdeiv2")

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Advance the state using the configured GSL odeiv2 stepper.
  !>
  !> On first invocation, initializes the GSL context with the state dimension
  !> and configured stepper type. Subsequent calls reuse the context.
  !>
  !> @param[inout] this   The GSL integrator instance.
  !> @param[inout] state  Complex state vector (modified in-place).
  !> @param[in]    t0     Starting time.
  !> @param[in]    t1     Target ending time.
  module subroutine Integrate(this, state, t0, t1)
    use M_Utils_Odeiv2GslLib

    class(T_IntegratorList_E_GslOdeiv2), intent(inout) :: this
    complex(R64), intent(inout), contiguous :: state(:)
    real(R64), intent(in) :: t0
    real(R64), intent(in) :: t1

    integer(I32) :: status
    integer :: nSteps

    ! Lazy initialization: create GSL context on first integration call
    if (.not. this % initializedQ) then
      this % dimension = size(state)
      this % step_size = (t1 - t0)

      call Odeiv2GslLib_CreateCtx(this % odeiv2Ctx, this % dimension, this % TimeDerivative, &
                                  this % step_size, DEFAULT_ABS_ERROR, DEFAULT_REL_ERROR, &
                                  trim(this % stepperType))

      this % initializedQ = .true.
    end if

    nSteps = 1
    call Odeiv2GslLib_Integrate(status, state, t0, t1, nSteps, this % odeiv2Ctx)

  end subroutine

end submodule
