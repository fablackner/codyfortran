! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Module M_CoeffsInit_Excited provides a concrete initializer that
!> constructs CI coefficients for excited-state configurations via
!> second-quantized creation and annihilation operators.
!>
!> Overview
!> --------
!> Builds the initial CI state by applying a ladder-operator string
!> â†_{p₁}...â†_{pₙ} â_{q₁}...â_{qₘ} to the reference configuration.
!> This enables preparation of particle-hole excitations, spin-flips,
!> or specific occupation patterns.
!>
!> JSON Parameters
!> ---------------
!>   creates    : array of orbital indices for â† operators
!>   destroys   : array of orbital indices for â  operators
!>   bodyType1  : first body type (default 1)
!>   bodyType2  : second body type (default 2)
!>
!> The same excitation is applied to both body types (e.g., α and β
!> electrons in a spin-unrestricted calculation).
!>
!> Wiring
!> ------
!> The factory `CoeffsInit_Excited_Fabricate` reads the JSON parameters
!> and binds `CoeffsInit_Initialize` in M_CoeffsInit.
module M_CoeffsInit_Excited
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Factory that registers the "Excited" initializer.
    !>
    !> Behavior
    !> --------
    !> - Reads JSON arrays `creates` and `destroys` specifying the
    !>   ladder operators to apply.
    !> - Reads optional `bodyType1` and `bodyType2` (defaults: 1, 2).
    !> - Binds `CoeffsInit_Initialize` to the excited-state builder.
    module subroutine CoeffsInit_Excited_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

end module
