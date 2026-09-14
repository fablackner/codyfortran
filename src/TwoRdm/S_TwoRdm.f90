! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Implementation submodule for M_TwoRdm.
!>
!> @details
!> Provides the fabrication entry point that reads the JSON configuration
!> and delegates to the appropriate particle-statistics backend (fermionic).
!>
!> Currently supported backends:
!> - `twoRdm.fermionic`: Fermionic 2-RDM with antisymmetric exchange statistics
!>
!> Future extension point for bosonic backends would add:
!> - `twoRdm.bosonic`: Bosonic 2-RDM with symmetric exchange statistics
submodule(M_TwoRdm) S_TwoRdm

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Fabricate the TwoRdm subsystem from JSON configuration.
!>
!> @details
!> Reads the `twoRdm` section of the JSON configuration and wires procedure
!> pointers to the selected backend implementation. This is the top-level
!> entry point for 2-RDM fabrication.
!>
!> JSON configuration example:
!> ```json
!> {
!>   "twoRdm": {
!>     "fermionic": {
!>       "generic": {}
!>     }
!>   }
!> }
!> ```
  module subroutine TwoRdm_Fabricate()
    use M_Utils_Json
    use M_Utils_Say
    use M_TwoRdm_Fermionic

    call Say_Fabricate("twoRdm")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    TwoRdm_SaveTwoRdmState => SaveTwoRdmState

    !------------------------------------
    ! branch on particle statistics
    !------------------------------------

    if (Json_GetExistence("twoRdm.fermionic")) then
      call TwoRdm_Fermionic_Fabricate
    else
      error stop "twoRdm is missing one of: fermionic"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Persist the linearized 2-RDM state to disk.
!>
!> @details
!> Writes the 2-RDM state vector to `twoRdmState.dat` in binary format.
!> The file can be used for checkpointing or post-processing analysis.
!>
!> @param[in] twoRdmState  Linearized 2-RDM vector to save
  subroutine SaveTwoRdmState(twoRdmState)
    use M_Utils_DataStorage

    complex(R64), intent(in), contiguous  :: twoRdmState(:)

    call SaveData(twoRdmState, 'twoRdmState.dat', storage_size(twoRdmState))

  end subroutine

end submodule

