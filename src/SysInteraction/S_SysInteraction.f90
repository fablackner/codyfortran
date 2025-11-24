! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysInteraction) S_Interaction

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
    ! branch
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
