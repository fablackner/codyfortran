! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> High-level method interface and runtime binding for quantum dynamics.
!>
!> This module defines the public contract that every propagation "method"
!> (e.g., TDHX, Full CI, MCTDHX, TD-2RDM, single-body) must satisfy. Concrete
!> method submodules bind the procedure pointers during their `..._Fabricate`
!> routines based on the chosen configuration, and initialize the packed state
!> vector held here. All methods operate on a contiguous, 1D complex state
!> vector with method-specific internal layout.
module M_Method
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Binds/overrides method procedure pointers at runtime and (optionally)
    !> prepares method-wide defaults based on parsed input configuration.
    !>
    !> Notes
    !> - Should be called once during initialization, after configuration is
    !>   available and before any time stepping.
    !> - Concrete submodules typically rebind `Method_Setup`,
    !>   `Method_TimeDerivative`, and `Method_GetEnergy`.
    module subroutine Method_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Packed quantum state in the chosen method's format.
  !>
  !> The 1D, contiguous complex state contains the full dynamical variables
  !> of the method (e.g., orbitals, CI coefficients, reduced density matrices).
  !> Its exact layout and allocation are determined by the bound `Method_Setup`
  !> routine. Ownership lives in this module; concrete methods should allocate
  !> and resize this array as needed during setup.
  complex(R64), allocatable, target :: Method_state(:)

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> Pointer to the setup procedure for initializing the selected method.
  procedure(I_Method_Setup), pointer :: Method_Setup => NoOpProcedures_Setup
  abstract interface
    !> Initializes the selected method.
    !>
    !> Responsibilities
    !> - Allocate and populate `Method_state` with the method-specific layout.
    !> - Bind any additional, method-specific procedure pointers if required.
    !> - Set initial conditions (ground state, user-provided state, etc.).
    !>
    !> Contract
    !> - Must be idempotent with respect to repeated calls or document otherwise.
    !> - On failure, should raise/stop with an informative message.
    subroutine I_Method_Setup
    end subroutine
  end interface

  !> Pointer to the procedure that computes the TDSE right-hand side.
  procedure(I_Method_TimeDerivative), pointer :: Method_TimeDerivative

  abstract interface
    !> Computes the right-hand side of the time-dependent Schrödinger equation
    !> in units where ħ=1:
    !> \[ \partial_t |\Psi\rangle = -i\,\hat{H}(t)\,|\Psi\rangle. \]
    !>
    !> Semantics
    !> - Overwrites `dState` with `-i * H(state, t)` (no accumulation).
    !> - Must not assume aliasing between `state` and `dState`.
    !> - May use module-level data or cached operator factorizations.
    subroutine I_Method_TimeDerivative(dState, state, time)
      import :: R64
      !> Output: time derivative evaluated at `time` (multiplied by -i).
      complex(R64), intent(out), contiguous, target :: dState(:)
      !> Input: quantum state vector at `time`.
      complex(R64), intent(in), contiguous, target  :: state(:)
      !> Current time (for time-dependent Hamiltonians/operators).
      real(R64), intent(in) :: time
    end subroutine
  end interface

  !> Pointer to the procedure for calculating the energy expectation value.
  procedure(I_Method_GetEnergy), pointer :: Method_GetEnergy
  abstract interface
    !> Returns the energy expectation value of the current quantum state
    !> at time `time`:
    !> \( E(t) = \langle \Psi(t) | \hat{H}(t) | \Psi(t) \rangle \).
    !>
    !> Note: Not every method provides a cheap energy evaluation. Callers
    !> should guard with `associated(Method_GetEnergy)` when optional.
    function I_Method_GetEnergy(time) result(res)
      import :: R64
      !> Energy expectation value at time `time`.
      real(R64) :: res
      !> Current time (for time-dependent Hamiltonians/operators).
      real(R64), intent(in) :: time
    end function
  end interface

end module
