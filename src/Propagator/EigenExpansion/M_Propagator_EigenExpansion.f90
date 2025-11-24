! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Eigen-expansion based propagator backend.
!>
!> This backend advances the state by expanding it in the eigenbasis of the
!> (possibly time-independent) Hamiltonian and applying phase factors in that
!> basis. It is most suitable when eigenpairs are available or inexpensive to
!> obtain for the system size of interest.
!>
!> Characteristics
!> - Works with a unified Hamiltonian action (no operator splitting).
!> - Propagation cost dominated by basis transforms and phase updates.
!> - Accuracy depends on the quality/availability of the eigen-decomposition.
module M_Propagator_EigenExpansion
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Fabricate and wire the eigen-expansion backend.
    !>
    !> Assigns `M_Propagator` procedure pointers to the eigen-expansion
    !> implementations. After calling this routine, use
    !> `M_Propagator::Propagator_Setup` to prepare eigen data if necessary, then
    !> call `M_Propagator::Propagator_Propagate` to advance the state.
    module subroutine Propagator_EigenExpansion_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module

