! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Module M_CoeffsInit_Load provides a concrete initializer that loads
!> CI coefficients from an external binary file.
!>
!> Overview
!> --------
!> Reads a pre-computed coefficient vector from disk (`coeffs.in`),
!> enabling simulation restarts, checkpointing, and initialization
!> from externally generated states.
!>
!> File Format
!> -----------
!> The input file must contain raw complex(R64) values in native
!> Fortran binary format (no headers). Use `SaveData` from
!> M_Utils_DataStorage to produce compatible files.
!>
!> Wiring
!> ------
!> The factory `CoeffsInit_Load_Fabricate` binds `CoeffsInit_Initialize`
!> in M_CoeffsInit to the load implementation. No setup phase is required.
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
    !> --------
    !> - Binds `CoeffsInit_Initialize` in `M_CoeffsInit` to read
    !>   coefficients from the binary file `coeffs.in`.
    !> - The file must exist in the working directory at runtime.
    !> - No additional JSON parameters are currently supported.
    module subroutine CoeffsInit_Load_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

end module
