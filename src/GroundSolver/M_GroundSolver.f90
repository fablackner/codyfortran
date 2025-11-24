! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Ground-state solver scaffolding and public interfaces.
!>
!> This module does not contain the numerical algorithms itself. Instead, it
!> declares procedure pointers that are bound at runtime by a corresponding
!> `GroundSolver_Fabricate` routine (provided by a backend/implementation).
!>
!> The bound procedures typically implement one SCF/relaxation iteration that
!> drives an input state towards a stationary (ground) state. Backends may use
!> different physical models or discretizations (e.g. TDHx, Ylm radial grids),
!> but they all conform to the same small interface declared here.
module M_GroundSolver
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface GroundSolver_Fabricate
    !> Bind concrete implementations to the procedure pointers exported by
    !> `M_GroundSolver`.
    !>
    !> Implementations of this routine are responsible for wiring
    !> `GroundSolver_Setup` and `GroundSolver_Approach` to backend
    !> specific procedures (for example a TDHx or Ylm variant). No state is
    !> created here beyond assigning the procedure pointers.
    module subroutine GroundSolver_Fabricate()
    end subroutine
  end interface

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> Pointer to the setup routine for the ground-state solver backend.
  procedure(I_GroundSolver_Setup), pointer :: GroundSolver_Setup => NoOpProcedures_Setup
  abstract interface
    !> Initialize the ground-state solver backend.
    !>
    !> Typical responsibilities include allocating working buffers, reading
    !> configuration, binding model-specific callbacks, and validating input
    !> dimensions. This routine performs no iterations by itself.
    subroutine I_GroundSolver_Setup
    end subroutine
  end interface

  !> Pointer to one ground-state approach/relaxation step.
  procedure(I_GroundSolver_Approach), pointer :: GroundSolver_Approach
  abstract interface
    !> Perform a single approach/relaxation iteration on the state.
    !>
    !> Implementations update the input state in place, typically using a
    !> residual and a mixing scheme to approach a self-consistent solution.
    !> The meaning of the `time` argument is backend-defined: some backends use
    !> it as a physical time parameter, others as an iteration/imaginary-time
    !> step used in the relaxation.
    subroutine I_GroundSolver_Approach(state, alpha, time)
      import :: R64
      !> Input/output state vector. Updated in place by the iteration.
      complex(R64), intent(inout), contiguous, target :: state(:)
      !> Linear/nonlinear mixing parameter used by the backend (0 < alpha <= 1).
      real(R64), intent(in) :: alpha
      !> Backend-defined step parameter (e.g., physical or imaginary time).
      real(R64), intent(in) :: time
    end subroutine
  end interface

end module

