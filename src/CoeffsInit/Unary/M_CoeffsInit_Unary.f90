! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Module M_CoeffsInit_Unary provides a concrete initializer that builds
!> CI coefficients for a unary (single-configuration/product-state)
!> initial condition.
!>
!> Overview
!> - Constructs a coefficient vector corresponding to a single selected
!>   configuration (or a simple product state) as specified in the input.
!> - This module's factory routine wires the unary implementation into the
!>   generic interface exported by `M_CoeffsInit`.
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
    !> - Reads the input configuration (JSON) that specifies which
    !>   configuration/product-state to populate.
    !> - Assigns `CoeffsInit_Setup` and `CoeffsInit_Initialize` in
    !>   `M_CoeffsInit` to the unary-specific implementations.
    module subroutine CoeffsInit_Unary_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

end module
