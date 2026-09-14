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
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
