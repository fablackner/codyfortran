! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file M_TwoRdmInit_Load.f90
!> @brief File-backed specialization of the 2-RDM initializer.
!>
!> @details
!> This module provides a 2-RDM initializer that reads tensor blocks from
!> external binary files. This is useful for:
!>
!>   - **Restart/continuation runs**: Resume propagation from saved state
!>   - **Pre-computed initial states**: Load externally prepared 2-RDMs
!>   - **Testing**: Load reference data for validation
!>
!> @par File Format
!> Each body-type block is stored in a separate binary file with naming
!> convention: `rdm2BtXX_YY.in` where XX and YY are zero-padded body type
!> indices (e.g., `rdm2Bt01_02.in` for the alpha-beta block).
!>
!> The binary format uses `storage_size(complex(R64))` for element sizing,
!> consistent with Fortran unformatted I/O via the DataStorage utilities.
!>
!> @par JSON Configuration
!> @code{.json}
!> {
!>   "twoRdmInit": {
!>     "load": { }
!>   }
!> }
!> @endcode
!>
!> No additional parameters are required. Files must exist in the working
!> directory with the expected naming convention.
!>
!> @par Required Files
!>   - `rdm2Bt01_01.in` - alpha-alpha block
!>   - `rdm2Bt01_02.in` - alpha-beta block
!>   - `rdm2Bt02_01.in` - beta-alpha block
!>   - `rdm2Bt02_02.in` - beta-beta block
module M_TwoRdmInit_Load
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> @brief Bind `M_TwoRdmInit` procedure pointers to the file-based
    !>        implementation.
    !>
    !> @details
    !> Currently reads no additional parameters from JSON; the empty
    !> `"twoRdmInit.load"` block simply selects this implementation.
    !> File paths are constructed from body-type indices at load time.
    !>
    !> Binds `TwoRdmInit_InitializeBlock` -> `InitializeBlock`
    module subroutine TwoRdmInit_Load_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

end module
