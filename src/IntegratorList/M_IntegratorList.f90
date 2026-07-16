! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file M_IntegratorList.f90
!> @brief Core ODE integrator framework for quantum state propagation.
!>
!> @details
!> This module defines the abstract API that all time-integration backends
!> implement. It follows the CodyFortranRDM **Interface Module + Submodule**
!> paradigm where:
!>
!> - **Interface modules** (`M_*.f90`) define contracts via abstract types
!>   and deferred procedure pointers.
!> - **Submodules** (`S_*.f90`) provide concrete algorithm implementations.
!>
!> ## Provided Facilities
!>
!> - A minimal contract for integrators (`Fabricate`, `Setup`, `Integrate`)
!> - A pointer to the user-supplied time-derivative routine (`TimeDerivative`)
!> - Lightweight containers for storing heterogeneous integrator elements
!> - A module-level fabrication entry point that builds the integrator list
!>   from JSON configuration at runtime
!>
!> ## Available Backends
!>
!> | Backend         | Method                                         |
!> |-----------------|------------------------------------------------|
!> | `Rk/O1Expl`     | Forward Euler (1st-order explicit)             |
!> | `Rk/O2Expl`     | Midpoint rule (2nd-order explicit)             |
!> | `Rk/O2Impl`     | Implicit midpoint (2nd-order, A-stable)        |
!> | `Rk/O4Expl`     | Classical RK4 (4th-order explicit)             |
!> | `CrankNicolson` | Crank–Nicolson (A-stable, norm-preserving)     |
!> | `Expokit`       | Krylov subspace matrix exponential             |
!> | `GslOdeiv2`     | GSL adaptive steppers (rkf45, rkck, rk8pd, …)  |
!> | `Sil`           | Short Iterative Lanczos                        |
!>
!> ## Lifecycle
!>
!> All integrators follow a three-stage lifecycle:
!> 1. **Fabricate** — Read parameters from JSON, bind procedure pointers.
!> 2. **Setup** — Allocate working memory once state dimensions are known.
!> 3. **Integrate** — Advance the quantum state from `t0` to `t1` in-place.
!>
!> @see AGENTS.md for onboarding documentation and extension guide.
module M_IntegratorList
  use M_Utils_Types
  use M_Utils_NoOpProcedures, only: NoOpProcedures_Setup

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  !> Input descriptor passed to `IntegratorList_Fabricate`.
  !>
  !> Bundles function pointers and other hooks that integrators need during
  !> construction. Each element in the input array corresponds to one
  !> integrator in the JSON `"integratorList"` array.
  !>
  !> @note At present this only carries the time derivative callback, but
  !>   the type is extensible for future hooks (e.g., preconditioning).
  type :: T_IntegratorList_FabricateInput
    !> Procedure that computes d(state)/dt at a given time.
    !> For quantum dynamics this is typically −i Ĥ |ψ⟩.
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

  !> Abstract base for all ODE integrator backends.
  !>
  !> Every integrator implements a **three-stage lifecycle**:
  !>
  !> | Stage       | When Called            | Purpose                                |
  !> |-------------|------------------------|----------------------------------------|
  !> | `Fabricate` | During initialization  | Parse JSON, bind procedure pointers    |
  !> | `Setup`     | After dimensions known | Allocate work arrays, precompute data  |
  !> | `Integrate` | During time loop       | Advance state from `t0` to `t1`        |
  !>
  !> Concrete implementations (RK, CN, Expokit, GSL, SIL) extend this type
  !> and provide algorithm-specific storage and logic.
  type, abstract :: T_IntegratorList_E
    !> Callback that computes d(state)/dt for a given state and time.
    !> Bound during fabrication; invoked by `Integrate`.
    procedure(I_TimeDerivative), pointer, nopass :: TimeDerivative

    !> JSON path to this integrator's configuration node (e.g.,
    !> `"integratorList.rk.o4Expl"`). Stored during fabrication for
    !> parameter retrieval and diagnostic messages.
    character(len=:), allocatable :: path

    !> Optional metric weights defining the inner product used by
    !> metric-aware integrators (currently SIL). When allocated, inner
    !> products are evaluated as sum(conjg(a) * metricWeights * b); when
    !> unallocated, the Euclidean product is used. This is required on grids
    !> whose state representation is not weight-absorbed (e.g., FEDVR),
    !> where Hermitian operators are self-adjoint only in the weighted
    !> metric. Assign directly on the registry element after fabrication and
    !> before the first `Integrate` call (grid weights typically become
    !> available only after `Grid_Setup`).
    real(R64), allocatable :: metricWeights(:)
  contains
    !> @brief Read parameters and capture dependencies from configuration.
    !> @details Called once after allocation. Implementations should read
    !>   JSON keys at `this%path` and store them in instance variables.
    procedure(I_Fabricate), deferred :: Fabricate

    !> @brief Allocate work arrays and initialize algorithmic state.
    !> @details Called once before time stepping begins. At this point
    !>   state dimensions are known.
    procedure(I_Setup), deferred :: Setup

    !> @brief Advance the quantum state from `t0` to `t1` in-place.
    !> @details The numerical scheme (RK, CN, Krylov, …) determines
    !>   step control, stability and accuracy.
    procedure(I_Integrate), deferred :: Integrate
  end type

  abstract interface
    !> Initialize an integrator instance from JSON configuration.
    !>
    !> Concrete types read parameters from the JSON node located at
    !> `this%path` and attach the derivative routine. This stage must
    !> **not** allocate large work arrays—defer that to `Setup`.
    !>
    !> @param[inout] this  The integrator instance to initialize.
    subroutine I_Fabricate(this)
      import :: T_IntegratorList_E
      class(T_IntegratorList_E), intent(inout) :: this
    end subroutine
  end interface

  abstract interface
    !> Perform post-fabrication setup.
    !>
    !> Allocate internal buffers and precompute constants needed by the
    !> integration algorithm. Called once before time stepping begins,
    !> at which point state dimensions are finalized.
    !>
    !> @param[inout] this  The integrator instance to set up.
    subroutine I_Setup(this)
      import :: T_IntegratorList_E
      class(T_IntegratorList_E), intent(inout) :: this
    end subroutine
  end interface

  abstract interface
    !> Advance the state by integrating the time derivative over [t0, t1].
    !>
    !> The concrete method (RK, CN, Krylov/Expokit, GSL stepper, SIL, …)
    !> determines step control, stability and accuracy. The state vector
    !> is modified **in-place**.
    !>
    !> @param[inout] this   The integrator instance.
    !> @param[inout] state  Complex state vector to propagate.
    !> @param[in]    t0     Starting time for integration.
    !> @param[in]    t1     Target ending time for integration.
    subroutine I_Integrate(this, state, t0, t1)
      import :: R64, T_IntegratorList_E
      class(T_IntegratorList_E), intent(inout) :: this
      complex(R64), intent(inout), contiguous :: state(:)
      real(R64), intent(in) :: t0
      real(R64), intent(in) :: t1
    end subroutine
  end interface

  abstract interface
    !> User-supplied time derivative callback (right-hand side of ODE).
    !>
    !> Computes **d(state)/dt** for a given state and time. In quantum
    !> dynamics this is typically:
    !>
    !>     dState = −i Ĥ(t) |ψ(t)⟩
    !>
    !> but the interface is general for any first-order ODE system.
    !>
    !> @param[out] dState  Output array containing the computed derivative.
    !> @param[in]  state   Input state at the current time.
    !> @param[in]  time    Current time for derivative evaluation.
    subroutine I_TimeDerivative(dState, state, time)
      import :: R64
      complex(R64), intent(out), contiguous, target :: dState(:)
      complex(R64), intent(in), contiguous, target  :: state(:)
      real(R64), intent(in)                         :: time
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Polymorphic container for heterogeneous integrator storage.
  !>
  !> Because Fortran arrays must be homogeneous, we wrap the polymorphic
  !> `T_IntegratorList_E` pointer in a thin container type so that
  !> `integratorList(:)` can hold different concrete integrators.
  type :: T_IntegratorList_Container
    !> Polymorphic instance of a concrete `T_IntegratorList_E` backend.
    class(T_IntegratorList_E), allocatable :: e
  end type

  !> Global registry of active integrators.
  !>
  !> Populated by `IntegratorList_Fabricate` from JSON configuration.
  !> The `Propagator` layer indexes into this array to invoke integrators.
  !> Each element holds one polymorphic backend (RK, CN, Expokit, …).
  type(T_IntegratorList_Container), allocatable :: integratorList(:)

  !=============================================================================
  ! module procedure pointers
  !=============================================================================

  !> Optional hook to perform module-level setup across all integrators.
  !>
  !> By default this points to a no-op. After `IntegratorList_Fabricate`
  !> runs, it is re-bound to an internal `Setup` routine that iterates over
  !> `integratorList(:)` and calls `element%Setup()` on each backend.
  procedure(I_IntegratorList_Setup), pointer :: IntegratorList_Setup => NoOpProcedures_Setup

  abstract interface
    !> Global integrator setup interface (iterates over all backends).
    subroutine I_IntegratorList_Setup
    end subroutine
  end interface

end module

