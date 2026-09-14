! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Purification_SingletTr0) S_Purification_SingletTr0

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Purification_SingletTr0_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Purification_SingletTr0_Reduce

    call Say_Fabricate("purification.singletTr0")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("purification.singletTr0.reduce")) then
      call Purification_SingletTr0_Reduce_Fabricate

    else
      error stop "purification.singletTr0 is missing one of: reduce"
    end if

  end subroutine

end submodule

