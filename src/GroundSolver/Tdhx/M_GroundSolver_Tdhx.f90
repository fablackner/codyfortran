! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> TDHx backend hooks for the ground-state solver.
!>
!> This module declares the TDHx-specific procedure pointers and factory
!> routine used to wire a concrete implementation. It does not contain
!> numerical kernels; those live in implementation modules that call the
!> corresponding `..._Fabricate` routine to bind the pointers exported here.
module M_GroundSolver_Tdhx
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Bind TDHx-specific callbacks used by the ground-state solver.
    !>
    !> Implementations of this routine assign `GroundSolver_Tdhx_HartreeFockAction`
    !> to the backend function
    module subroutine GroundSolver_Tdhx_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> Pointer to the TDHx Hartree–Fock action used in one iteration.
  procedure(I_GroundSolver_Tdhx_HartreeFockAction), pointer :: GroundSolver_Tdhx_HartreeFockAction
  abstract interface
    !> Compute the TDHx Hartree–Fock action on an orbital.
    !>
    !> Given an input orbital `orb`, compute `dOrb = A_TDHx(orb, time)`, where
    !> `A_TDHx` denotes the effective (model-dependent) action used by the
    !> ground-state iteration. The exact form (e.g., potentials, exchange) is
    !> provided by the backend.
    subroutine I_GroundSolver_Tdhx_HartreeFockAction(dOrb, orb, time)
      import :: R64
      !> Output action applied to the orbital.
      complex(R64), intent(out), contiguous, target :: dOrb(:)
      !> Input orbital/state vector.
      complex(R64), intent(in), contiguous, target :: orb(:)
      !> Backend-defined step parameter (physical or imaginary time).
      real(R64), intent(in) :: time
    end subroutine
  end interface

end module

