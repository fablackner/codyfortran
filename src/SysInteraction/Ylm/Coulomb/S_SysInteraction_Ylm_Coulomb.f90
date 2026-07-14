! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Implementation submodule for Ylm Coulomb interaction.
!>
!> @details Reads the coupling strength and dispatches to the selected radial
!> Poisson solver. Sets `mIndependentQ = .true.` since 1/r is spherically
!> symmetric (the Vₗₘ depends only on l, not m).
submodule(M_SysInteraction_Ylm_Coulomb) S_SysInteraction_Ylm_Coulomb

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Fabricate the Coulomb interaction for Ylm grids.
  !>
  !> Reads `strength` parameter (default 1.0) and selects the radial solver.
  module subroutine SysInteraction_Ylm_Coulomb_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysInteraction
    use M_SysInteraction_Ylm
    use M_SysInteraction_Ylm_Coulomb_StdImpl
    use M_SysInteraction_Ylm_Coulomb_FullEq
    use M_SysInteraction_Ylm_Coulomb_BlockEq
    use M_SysInteraction_Ylm_Coulomb_TwoScan

    call Say_Fabricate("sysInteraction.ylm.coulomb")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysInteraction_Ylm_Coulomb_Strength = Json_Get("sysInteraction.ylm.coulomb.strength", 1.0_R64)
    SysInteraction_Ylm_mIndependentQ = .true.

    ! Real Coulomb kernel: swapped source pairs give the conjugated potential,
    ! but only on a real radial contour (complex scaling breaks the symmetry)
    SysInteraction_conjSymmetricQ = .not. Json_GetExistence("grid.ylm.fedvrEcs")

    !------------------------------------
    ! branch: select radial Poisson solver
    !------------------------------------

    if (Json_GetExistence("sysInteraction.ylm.coulomb.stdImpl")) then
      call SysInteraction_Ylm_Coulomb_StdImpl_Fabricate

    else if (Json_GetExistence("sysInteraction.ylm.coulomb.fullEq")) then
      call SysInteraction_Ylm_Coulomb_FullEq_Fabricate

    else if (Json_GetExistence("sysInteraction.ylm.coulomb.blockEq")) then
      call SysInteraction_Ylm_Coulomb_BlockEq_Fabricate

    else if (Json_GetExistence("sysInteraction.ylm.coulomb.twoScan")) then
      call SysInteraction_Ylm_Coulomb_TwoScan_Fabricate

    else
      error stop "sysInteraction.ylm.coulomb is missing one of: stdImpl, fullEq, blockEq, twoScan"
    end if

  end subroutine

end submodule
