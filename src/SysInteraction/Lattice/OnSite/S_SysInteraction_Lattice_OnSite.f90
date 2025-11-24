! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysInteraction_Lattice_OnSite) S_SysInteraction_Lattice_OnSite

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysInteraction_Lattice_OnSite_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysInteraction
    use M_SysInteraction_Lattice_OnSite_StdImpl

    call Say_Fabricate("sysInteraction.lattice.onSite")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysInteraction_Lattice_OnSite_Strength = Json_Get("sysInteraction.lattice.onSite.strength", 1.0_R64)

    SysInteraction_timeIndependentQ = .true.
    SysInteraction_bodyTypeIndependentQ = .true.

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("sysInteraction.lattice.onSite.stdImpl")) then
      call SysInteraction_Lattice_OnSite_StdImpl_Fabricate

    else
      error stop "sysInteraction.lattice.onSite is missing one of: stdImpl"
    end if

  end subroutine

end submodule
