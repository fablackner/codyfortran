! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> TD-2RDM with dynamic orbitals (orbital-based).
!>
!> The time evolution advances both the two-body reduced density matrix and the
!> single-particle orbitals, allowing the basis to adapt to the changing state.
!> This can improve accuracy for strongly time-dependent or correlated systems.
!>
!> State layout
!> - The packed state typically stores orbitals and a compressed representation
!>   of the two-body RDM consistent with symmetry and particle statistics.
module M_Method_Mb_OrbBased_Td2rdm
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Binds TD-2RDM (dynamic orbitals) procedure pointers at runtime and sets up
    !> state layout. Typically sets `Method_Setup` and `Method_TimeDerivative`.
    module subroutine Method_Mb_OrbBased_Td2rdm_Fabricate
    end subroutine

    !> Active/weak natural-orbital split only (no-op for nWeakOrbs = 0):
    !> recomputes the weak natural orbitals that the linear split-step part
    !> uses from the current Method_state. The propagator refreshes them at
    !> every evaluation of the non-linear part (validated for splitStep order2
    !> with an RK non-linear integrator); call this after the state was
    !> changed outside the propagator (absorber, purification, a kick) and
    !> before the next Propagator_Propagate, otherwise the next linear step
    !> uses a stale weak space and the split step drops to first order.
    module subroutine Method_Mb_OrbBased_Td2rdm_RefreshWeakOrbs
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Number of lowest natural orbitals per body type that are weak, i.e. get
  !> no Q-space drive (method.mb.orbBased.td2rdm.nWeakOrbs; 0 = off, the
  !> standard TD-2RDM orbital equations).
  integer(I32) :: Method_Mb_OrbBased_Td2rdm_nWeakOrbs = 0

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
