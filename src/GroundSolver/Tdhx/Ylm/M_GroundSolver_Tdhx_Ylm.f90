! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Ylm (spherical-harmonics) TDHx hooks for the ground-state solver.
!>
!> This variant exposes interfaces specialized for a radial representation
!> labelled by the angular momentum quantum number `l`. It provides a factory
!> to bind a concrete backend and the signature of the Hartree–Fock radial
!> action applied during one iteration.
module M_GroundSolver_Tdhx_YlmOpt
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Bind Ylm-specific TDHx callbacks for the ground-state solver.
    !>
    !> Implementations assign the radial Hartree–Fock action pointer below to a
    !> backend routine operating on a single `(l)` channel.
    module subroutine GroundSolver_Tdhx_YlmOpt_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> Pointer to the Ylm TDHx Hartree–Fock radial action used per `(l)` channel.
  procedure(I_GroundSolver_Tdhx_YlmOpt_HartreeFockRadialAction), pointer :: GroundSolver_Tdhx_YlmOpt_HartreeFockAction
  abstract interface
    !> Compute the TDHx Hartree–Fock radial action for a given `l`.
    !>
    !> Given an input radial orbital `orbLm` for angular momentum `l`, compute
    !> `dOrbLm = A_TDHx^l(orbLm, time)`, where the operator is defined by the
    !> selected backend and boundary conditions.
    subroutine I_GroundSolver_Tdhx_YlmOpt_HartreeFockRadialAction(dOrbLm, orbLm, l, time)
      import :: I32, R64
      !> Output radial action applied to the orbital for channel `l`.
      complex(R64), intent(out), contiguous, target :: dOrbLm(:)
      !> Input radial orbital/state for channel `l`.
      complex(R64), intent(in), contiguous, target :: orbLm(:)
      !> Angular momentum quantum number (channel index).
      integer(I32), intent(in) :: l
      !> Backend-defined step parameter (physical or imaginary time).
      real(R64), intent(in) :: time
    end subroutine
  end interface

end module

