! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Core propagator facade for time evolution of quantum states.
!>
!> This module exposes the public propagator API via procedure pointers that are
!> assigned at runtime by a backend-specific "Fabricate" routine. Concrete
!> backends live in submodules such as `M_Propagator_Single`,
!> `M_Propagator_SplitStep[_Order{2,4}]`, or `M_Propagator_EigenExpansion` and
!> provide a `*_Fabricate` routine that wires the procedure pointers below to
!> their implementations.
!>
!> What this module provides
!> - `Propagator_Setup`: optional initialization hook for allocating work
!>   buffers, reading configuration, and preparing operators. Defaults to a no-op
!>   and is replaced by the chosen backend.
!> - `Propagator_Propagate(state, t0, t1)`: evolves the state from time `t0` to
!>   `t1` in-place according to the backend’s propagation scheme.
!>
!> How to select a backend
!> - Call exactly one backend `*_Fabricate` routine during program startup to
!>   assign the procedure pointers. Calling a different backend’s fabricate later
!>   will rewire the pointers accordingly.
!>
!> Notes
!> - The state vector must be contiguous and use the project-wide real kind
!>   `R64` for its complex components.
!> - Time direction is inferred from `(t1 - t0)` and may be forward or backward.
!> - Backends may manage internal step sizes, work arrays, and error control.
module M_Propagator
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface Propagator_Fabricate
    !> Backend selection entry point (facade contract).
    !>
    !> A concrete backend implements a module procedure with this binding name
    !> and, when called, assigns `Propagator_Setup` and `Propagator_Propagate` to
    !> its own implementations. This routine should be idempotent and safe to
    !> call multiple times; repeated calls simply reassign the procedure
    !> pointers.
    module subroutine Propagator_Fabricate()
    end subroutine
  end interface

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> Pointer to the setup procedure for initializing the propagator system.
  !>
  !> Typical responsibilities of a backend setup include:
  !> - Validating configuration and required operator dimensions
  !> - Allocating workspaces and precomputing constants (e.g., split-step factors)
  !> - Preparing callbacks to Hamiltonian actions or RHS evaluators
  procedure(I_Propagator_Setup), pointer :: Propagator_Setup => NoOpProcedures_Setup
  abstract interface
    !> Initializes the currently selected propagator backend.
    !>
    !> Call this once after fabrication and before the first propagation. It may
    !> be a no-op depending on the backend. Errors should be raised via the
    !> project’s error handling utilities if preconditions are not met.
    subroutine I_Propagator_Setup
    end subroutine
  end interface

  !> Pointer to the main propagation procedure.
  procedure(I_Propagator_Propagate), pointer :: Propagator_Propagate
  abstract interface
    !> In-place time propagation of a quantum state.
    !>
    !> Mathematical contract
    !> \[ |\Psi(t_1)\rangle = \hat{U}(t_1,t_0) |\Psi(t_0)\rangle \]
    !> where the backend defines how the time-evolution operator \(\hat{U}\) is
    !> constructed (e.g., split-operator composition, single-step integrator,
    !> or eigen-expansion).
    !>
    !> Arguments
    !> - `state`: complex, contiguous 1D array representing the state vector;
    !>   updated in-place to time `t1`.
    !> - `t0`: starting time.
    !> - `t1`: target time (may be less than `t0` for backward propagation).
    subroutine I_Propagator_Propagate(state, t0, t1)
      import :: R64
      !> Input/output quantum state. On exit, contains the state at time `t1`.
      complex(R64), intent(inout), contiguous :: state(:)
      !> Starting time for propagation.
      real(R64), intent(in) :: t0
      !> Target ending time for propagation.
      real(R64), intent(in) :: t1
    end subroutine
  end interface

end module

