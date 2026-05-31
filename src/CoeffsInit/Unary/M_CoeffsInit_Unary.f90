! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Module M_CoeffsInit_Unary provides a concrete initializer that builds
!> CI coefficients for a unary (single-configuration/product-state)
!> initial condition.
!>
!> Overview
!> --------
!> Constructs a coefficient vector c(1) = 1, c(i>1) = 0, corresponding
!> to the reference Fock configuration |Φ₀⟩ — the first configuration
!> in the enumeration defined by ConfigList.
!>
!> This is the simplest and most common initializer, used whenever
!> the simulation should start from a product state (Slater determinant
!> for fermions, permanents for bosons).
!>
!> Wiring
!> ------
!> The factory `CoeffsInit_Unary_Fabricate` binds `CoeffsInit_Initialize`
!> in M_CoeffsInit to the unary implementation. No setup phase is required.
module M_CoeffsInit_Unary
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Factory that registers the "Unary" initializer.
    !>
    !> Behavior
    !> --------
    !> - Binds `CoeffsInit_Initialize` in `M_CoeffsInit` to produce
    !>   the coefficient vector [1, 0, 0, ..., 0].
    !> - No JSON parameters are required; the `"unary": {}` block
    !>   signals the selection of this initializer.
    module subroutine CoeffsInit_Unary_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

end module
