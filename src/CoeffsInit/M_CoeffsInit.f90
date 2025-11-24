! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Module M_CoeffsInit defines the public interface for initializing
!> CI coefficients and exposes procedure pointers that are bound at
!> runtime.
!>
!> Overview
!> - `CoeffsInit_Fabricate` selects and wires a concrete initializer
!>   (e.g., unary/single-configuration, excited-state, or load-from-file)
!>   based on the input configuration (JSON) and assigns the procedure
!>   pointers below.
!> - `CoeffsInit_Setup` performs any one-time preparation required by the
!>   chosen initializer (e.g., building masks, reading parameters).
!> - `CoeffsInit_Initialize` fills the CI coefficient vector for the
!>   requested initial state.
!>
!> Notes
!> - This module only declares abstract interfaces and procedure pointers.
!>   The concrete implementations live in specialized modules and are
!>   connected via the `..._Fabricate` routines.
!> - The concrete initializer is selected at runtime from the input file;
!>   if no initializer is configured, a no-op placeholder is used.
module M_CoeffsInit
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Factory that selects a concrete coefficients initializer and
    !> assigns the procedure pointers exported by this module.
    !>
    !> Behavior
    !> - Parses the input configuration (JSON) to determine which
    !>   initializer to use (e.g., Unary, Excited, Load).
    !> - Binds `CoeffsInit_Setup` and `CoeffsInit_Initialize` to the
    !>   chosen implementation.
    !> - May read and cache any initializer-specific parameters.
    module subroutine CoeffsInit_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> Pointer to the setup procedure for the selected initializer.
  !>
  !> Contract
  !> - Must be called once after `CoeffsInit_Fabricate` and before
  !>   `CoeffsInit_Initialize`.
  !> - Performs any one-time preprocessing (e.g., reading parameters,
  !>   constructing masks, allocating caches).
  procedure(I_CoeffsInit_Setup), pointer :: CoeffsInit_Setup => NoOpProcedures_Setup
  abstract interface
    !> Initializes the selected coefficients initializer.
    !>
    !> Effects
    !> - Prepares internal state required by the concrete initializer.
    !> - Should be idempotent; repeated calls have no adverse effects.
    subroutine I_CoeffsInit_Setup
    end subroutine
  end interface

  !> Pointer to the procedure that produces the initial CI coefficients.
  !>
  !> Contract
  !> - Must be called after a successful `CoeffsInit_Setup`.
  !> - Fills the provided coefficient vector with the initial state as
  !>   defined by the chosen initializer.
  procedure(I_CoeffsInit_Initialize), pointer :: CoeffsInit_Initialize
  abstract interface
    !> Fills the CI coefficient vector with the configured initial state.
    !>
    !> @param coeffs [out]
    !>   Complex vector of length equal to the CI space dimension. On
    !>   return it contains the initial amplitudes. The array must be
    !>   contiguous. Implementations may assume the array is correctly
    !>   sized; otherwise they should signal an error.
    subroutine I_CoeffsInit_Initialize(coeffs)
      import :: R64
      !> Output array containing the initialized CI coefficients.
      complex(R64), intent(out), contiguous :: coeffs(:)
    end subroutine
  end interface

end module
