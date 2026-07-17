! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysGauge) S_SysGauge

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysGauge_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysGauge_Ylm

    implicit none

    call Say_Fabricate("sysGauge")

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("sysGauge.ylm")) then
      call SysGauge_Ylm_Fabricate

    else
      error stop "sysGauge is missing one of: ylm"
    end if

  end subroutine

end submodule
