! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Time-Dependent Hartree–Fock/Exchange (TDHF/TDHX) with time-dependent orbitals.
!>
!> Approximates the state by a single Slater determinant (fermions) or product state
!> (bosons) with time-dependent orbitals. Captures mean-field and exchange effects
!> while neglecting higher-order correlations, offering computational efficiency.
!>
!> State layout
!> - The packed state contains the set of time-dependent orbitals only; no CI
!>   coefficients are stored beyond what is implied by the single-determinant ansatz.
module M_Method_Mb_OrbBased_Tdhx
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Binds TDHF/TDHX procedure pointers at runtime and initializes the method
    !> state. Typically sets `Method_Setup`, `Method_GetEnergy`, and `Method_TimeDerivative`.
    module subroutine Method_Mb_OrbBased_Tdhx_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
