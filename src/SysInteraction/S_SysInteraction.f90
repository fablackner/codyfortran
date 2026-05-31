! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Fabrication submodule for the SysInteraction interface.
!>
!> @details Implements the top-level factory that dispatches to grid-specific
!> interaction back-ends based on the JSON configuration. The dispatch logic
!> mirrors the grid hierarchy:
!>   - `sysInteraction.linear`  → 1D real-space convolutions
!>   - `sysInteraction.lattice` → discrete tight-binding models
!>   - `sysInteraction.ylm`     → spherical-harmonic Coulomb solvers
submodule(M_SysInteraction) S_Interaction

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Top-level fabrication entry point for SysInteraction.
  !>
  !> Reads the JSON path `sysInteraction.*` to determine which grid back-end
  !> to activate, then delegates to the corresponding `*_Fabricate` routine.
  !>
  !> @note Exactly one of `linear`, `lattice`, or `ylm` must be present.
  module subroutine SysInteraction_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Utils_NoOpProcedures
    use M_SysInteraction_Linear
    use M_SysInteraction_Lattice
    use M_SysInteraction_Ylm

    call Say_Fabricate("sysInteraction")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    !------------------------------------
    ! branch: select grid-specific back-end
    !------------------------------------

    if (Json_GetExistence("sysInteraction.linear")) then
      call SysInteraction_Linear_Fabricate

    else if (Json_GetExistence("sysInteraction.lattice")) then
      call SysInteraction_Lattice_Fabricate

    else if (Json_GetExistence("sysInteraction.ylm")) then
      call SysInteraction_Ylm_Fabricate

    else
      error stop "sysInteraction is missing one of: linear, lattice, ylm"

    end if

  end subroutine

end submodule
