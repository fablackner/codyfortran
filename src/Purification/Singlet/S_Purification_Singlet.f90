! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Purification_Singlet) S_Purification_Singlet

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Purification_Singlet_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Purification_Singlet_Reduce

    call Say_Fabricate("purification.singlet")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("purification.singlet.reduce")) then
      call Purification_Singlet_Reduce_Fabricate

    else
      error stop "purification.singlet is missing one of: reduce"
    end if

  end subroutine

end submodule

