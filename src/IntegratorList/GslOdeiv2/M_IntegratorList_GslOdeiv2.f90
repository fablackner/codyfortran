! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> GSL odeiv2 stepper integrator backend.
!>
!> This module provides a concrete `T_IntegratorList_E` that wraps the
!> GNU Scientific Library (GSL) odeiv2 steppers (e.g., rk4, rkf45, rkck, rk8pd)
!> for advancing the state. The right-hand side is supplied via the
!> `TimeDerivative` callback.
module M_IntegratorList_GslOdeiv2
  use M_Utils_Types
  use M_IntegratorList, only: T_IntegratorList_E
  use M_Utils_Odeiv2GslLib, only: T_Odeiv2GslLib_Ctx

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Allocate a GSL odeiv2 integrator element.
    module subroutine IntegratorList_GslOdeiv2_Allocate(e, path)
      class(T_IntegratorList_E), allocatable, intent(out) :: e
      character(len=*), intent(in) :: path
    end subroutine
  end interface

  !=============================================================================
  ! module types
  !=============================================================================

  !> GSL odeiv2 integrator element.
  type, extends(T_IntegratorList_E) :: T_IntegratorList_E_GslOdeiv2
    !> Type of stepper algorithm from GSL to use.
    !> Common values: "rk4" (classical), "rkf45" (Fehlberg),
    !> "rkck" (Cash–Karp), "rk8pd" (Prince–Dormand).
    character(len=64) :: stepperType
    !> System dimension (number of complex components represented).
    integer(I32) :: dimension = 0
    !> Step size for GSL integrator.
    real(R64) :: step_size = 0.0_R64
    !> Whether the GSL context has been initialized.
    logical :: initializedQ = .false.
    !> GSL ODE context (workspaces, stepper and control objects).
    type(T_Odeiv2GslLib_Ctx) :: odeiv2Ctx
  contains
    !> Read stepper type, tolerances and sizes from configuration.
    procedure :: Fabricate
    !> Create and initialize the GSL stepper/control/driver.
    procedure :: Setup
    !> Advance the state using the selected GSL odeiv2 stepper.
    procedure :: Integrate
  end type

  interface
    !> Fabricate the GSL integrator from configuration.
    module subroutine Fabricate(this)
      !> The GSL integrator instance to initialize
      class(T_IntegratorList_E_GslOdeiv2), intent(inout) :: this
    end subroutine
  end interface

  interface
    !> Perform post-fabrication setup (allocate/initialize GSL context).
    module subroutine Setup(this)
      !> The GSL integrator instance to set up
      class(T_IntegratorList_E_GslOdeiv2), intent(inout) :: this
    end subroutine
  end interface

  interface
    !> Integrate from `t0` to `t1` using the configured GSL stepper.
    !> Uses the `TimeDerivative` callback for RHS evaluations.
    module subroutine Integrate(this, state, t0, t1)
      !> The GSL integrator instance
      class(T_IntegratorList_E_GslOdeiv2), intent(inout) :: this
      !> Quantum state to propagate (modified in-place)
      complex(R64), intent(inout), contiguous :: state(:)
      !> Initial time
      real(R64), intent(in) :: t0
      !> Final time
      real(R64), intent(in) :: t1
    end subroutine
  end interface

end module
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
