! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> TD-2RDM with fixed orbitals (orbital-based).
!>
!> Evolves only the two-body reduced density matrix in time while keeping the
!> single-particle orbitals fixed. This reduces computational cost and is suitable
!> when orbital relaxation is negligible or a pre-optimized basis is available.
!>
!> State layout
!> - The packed state contains the two-body RDM only; orbitals are provided
!>   externally or through module-level storage and remain constant in time.
module M_Method_Mb_OrbBased_Td2rdmStaticOrbs
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Binds TD-2RDM (fixed orbitals) procedure pointers at runtime and sets up the
    !> method state. Typically sets `Method_Setup` and `Method_TimeDerivative`.
    module subroutine Method_Mb_OrbBased_Td2rdmStaticOrbs_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
