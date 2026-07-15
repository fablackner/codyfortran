! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Implementation submodule for the prolate Coulomb interaction.
submodule(M_SysInteraction_Prolate_Coulomb) S_SysInteraction_Prolate_Coulomb

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Fabricate the Coulomb interaction for prolate grids.
  module subroutine SysInteraction_Prolate_Coulomb_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Grid_Prolate
    use M_SysInteraction
    use M_SysInteraction_Prolate_Coulomb_StdImpl

    call Say_Fabricate("sysInteraction.prolate.coulomb")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysInteraction_Prolate_Coulomb_strength = Json_Get("sysInteraction.prolate.coulomb.strength", 1.0_R64)
    SysInteraction_Prolate_Coulomb_lmax = Json_Get("sysInteraction.prolate.coulomb.lmax", Grid_Prolate_nEta - 1)

    if (SysInteraction_Prolate_Coulomb_lmax > Grid_Prolate_nEta - 1) then
      error stop "sysInteraction.prolate.coulomb.lmax exceeds nEta - 1 (unresolvable by the eta quadrature)"
    end if

    SysInteraction_timeIndependentQ = .true.
    SysInteraction_bodyTypeIndependentQ = .true.

    ! Real Coulomb kernel on a real contour: swapped source pairs give the
    ! conjugated potential
    SysInteraction_conjSymmetricQ = .true.

    !------------------------------------
    ! branch: select xi solver
    !------------------------------------

    if (Json_GetExistence("sysInteraction.prolate.coulomb.stdImpl")) then
      call SysInteraction_Prolate_Coulomb_StdImpl_Fabricate

    else
      error stop "sysInteraction.prolate.coulomb is missing one of: stdImpl"
    end if

  end subroutine

end submodule
