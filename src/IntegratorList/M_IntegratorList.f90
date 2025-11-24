! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Core integrator framework.
!>
!> This module defines the abstract API that all time-integration backends
!> implement. It provides:
!> - a minimal contract for integrators (`Fabricate`, `Setup`, `Integrate`),
!> - a pointer to the user-supplied time-derivative routine,
!> - lightweight containers for storing heterogeneous integrator elements, and
!> - a module-level fabrication entry point used to build the list of
!>   active integrators from configuration.
!>
!> Concrete implementations live in submodules (e.g. Runge–Kutta variants,
!> GSL ODE steppers, Crank–Nicolson, Expokit, Short Iterative Lanczos) and
!> extend `T_IntegratorList_E`.
module M_IntegratorList
  use M_Utils_Types
  use M_Utils_NoOpProcedures, only: NoOpProcedures_Setup

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  !> Input descriptor passed to `IntegratorList_Fabricate`.
  !> It bundles function pointers and other hooks that integrators need
  !> during construction. At present it carries the time derivative.
  type :: T_IntegratorList_FabricateInput
    !> Procedure that calculates the time derivative d(state)/dt at a given time.
    procedure(I_TimeDerivative), pointer, nopass :: TimeDerivative
  end type

  interface
    !> Builds the list of active integrators from configuration.
    !>
    !> Implementations allocate concrete `T_IntegratorList_E` objects,
    !> attach the provided derivative callback, and store them in the
    !> module-global `integratorList` container for later use.
    module subroutine IntegratorList_Fabricate(input)
      !> Array of construction inputs (at least the time derivative callback).
      type(T_IntegratorList_FabricateInput), intent(in) :: input(:)
    end subroutine
  end interface

  !=============================================================================
  ! type definition
  !=============================================================================

  !> Abstract base for all integrator backends.
  !>
  !> An integrator implements a three-stage lifecycle:
  !> - `Fabricate`: read parameters from configuration and store pointers/paths,
  !> - `Setup`: allocate working memory and initialize algorithmic state, and
  !> - `Integrate`: advance the state from `t0` to `t1` in-place using the
  !>   `TimeDerivative` callback.
  type, abstract :: T_IntegratorList_E
    !> Procedure callback that computes d(state)/dt at a given time.
    procedure(I_TimeDerivative), pointer, nopass :: TimeDerivative

    !> Path to this object's configuration node (used during fabrication).
    character(len=:), allocatable :: path
  contains
    !> Read parameters and capture dependencies from configuration.
    procedure(I_Fabricate), deferred :: Fabricate
    !> Allocate work arrays and initialize algorithmic state.
    procedure(I_Setup), deferred :: Setup
    !> Advance the quantum state from `t0` to `t1` using the chosen method.
    procedure(I_Integrate), deferred :: Integrate
  end type

  abstract interface
    !> Initialize an integrator instance from configuration.
    !>
    !> Concrete types read parameters from the JSON (or equivalent) located at
    !> `this%path` and attach the derivative routine.
    subroutine I_Fabricate(this)
      import :: T_IntegratorList_E
      !> The integrator instance to initialize
      class(T_IntegratorList_E), intent(inout) :: this
    end subroutine
  end interface

  abstract interface
    !> Perform post-fabrication setup.
    !>
    !> Allocate internal buffers and precompute constants needed by the
    !> integration algorithm. Called once before time stepping begins.
    subroutine I_Setup(this)
      import :: T_IntegratorList_E
      !> The integrator instance to set up
      class(T_IntegratorList_E), intent(inout) :: this
    end subroutine
  end interface

  abstract interface
    !> Advance the state by integrating the provided time derivative.
    !>
    !> The concrete method (RK, CN, Krylov/Expokit, GSL stepper, SIL, …)
    !> determines step control, stability and accuracy.
    subroutine I_Integrate(this, state, t0, t1)
      import :: R64, T_IntegratorList_E
      !> The integrator instance
      class(T_IntegratorList_E), intent(inout) :: this
      !> State vector to be propagated (modified in-place)
      complex(R64), intent(inout), contiguous :: state(:)
      !> Starting time for integration
      real(R64), intent(in) :: t0
      !> Target ending time for integration
      real(R64), intent(in) :: t1
    end subroutine
  end interface

  abstract interface
    !> User-supplied time derivative callback.
    !>
    !> Computes d(state)/dt for a given state and time. In quantum dynamics
    !> this is typically -i H |psi(t)>, but the interface is general.
    subroutine I_TimeDerivative(dState, state, time)
      import :: R64
      !> Output array containing the calculated time derivative.
      complex(R64), intent(out), contiguous, target :: dState(:)
      !> Input state at the current time.
      complex(R64), intent(in), contiguous, target  :: state(:)
      !> Current time for the derivative evaluation.
      real(R64), intent(in)             :: time
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Container used to hold heterogeneously-typed integrators.
  type :: T_IntegratorList_Container
    !> Polymorphic instance of a `T_IntegratorList_E` implementation.
    class(T_IntegratorList_E), allocatable :: e
  end type

  !> Global registry of active integrators.
  type(T_IntegratorList_Container), allocatable :: integratorList(:)

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> Optional hook to perform module-level setup across all integrators.
  procedure(I_IntegratorList_Setup), pointer :: IntegratorList_Setup => NoOpProcedures_Setup
  abstract interface
    !> Perform any global integrator setup (no-op by default).
    subroutine I_IntegratorList_Setup
    end subroutine
  end interface

end module

