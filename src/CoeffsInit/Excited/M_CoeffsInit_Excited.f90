! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Module M_CoeffsInit_Excited provides a concrete initializer that
!> constructs CI coefficients for an excited state configuration.
!>
!> Overview
!> - Intended to prepare initial amplitudes corresponding to a selected
!>   excited state (e.g., by quantum numbers, index, or mask).
!> - This module's factory routine wires the excited-state implementation
!>   into the generic interface exported by `M_CoeffsInit`.
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
    !> - Reads the input configuration (JSON) for excited-state parameters
    !>   (e.g., state index, spin/symmetry filters, or masks).
    !> - Assigns `CoeffsInit_Setup` and `CoeffsInit_Initialize` in
    !>   `M_CoeffsInit` to the excited-state specific implementations.
    module subroutine CoeffsInit_Excited_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

end module
