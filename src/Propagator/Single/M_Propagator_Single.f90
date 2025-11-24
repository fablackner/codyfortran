! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Single-RHS integrator propagator backend.
!>
!> This backend treats time evolution as integration of a single right-hand
!> side (RHS) function of the form \(\dot{\Psi} = f(t, \Psi)\), e.g.
!> \(\dot{\Psi} = -i\,\hat{H}(t)\,\Psi\). It couples a unified time
!> derivative with a chosen time integrator (e.g., RK, exponential integrator,
!> or similar) without operator splitting.
!>
!> Characteristics
!> - Uses a single Hamiltonian (or RHS) action per substep.
!> - Flexible choice of integrator and step control (backend-defined).
!> - Good default when no special structure is exploited.
module M_Propagator_Single
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Fabricate and wire the single-RHS backend.
    !>
    !> Assigns `M_Propagator` procedure pointers to single-integrator
    !> implementations. After calling this routine, invoke
    !> `M_Propagator::Propagator_Setup` once and then use
    !> `M_Propagator::Propagator_Propagate` for time evolution.
    module subroutine Propagator_Single_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module

