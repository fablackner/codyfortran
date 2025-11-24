! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Module M_CoeffsInit_Load provides a concrete initializer that loads
!> CI coefficients from an external file.
!>
!> Overview
!> - Reads persisted amplitudes (format defined by configuration) and
!>   uses them as the initial CI state, optionally normalizing or
!>   validating dimensions.
!> - This module's factory routine wires the loader implementation into
!>   the generic interface exported by `M_CoeffsInit`.
module M_CoeffsInit_Load
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Factory that registers the "Load" initializer.
    !>
    !> Behavior
    !> - Parses input configuration (JSON) for the source file path and
    !>   related options (format, normalization, slicing, etc.).
    !> - Assigns `CoeffsInit_Setup` and `CoeffsInit_Initialize` in
    !>   `M_CoeffsInit` to the loader-specific implementations.
    module subroutine CoeffsInit_Load_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

end module
