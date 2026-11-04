! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysGauge_Fem) S_SysGauge_Fem

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysGauge_Fem_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysGauge_Fem_Length

    call Say_Fabricate("sysGauge.fem")

    !------------------------------------
    ! validate grid requirements
    !------------------------------------

    if (.not. Json_GetExistence("grid.fem")) then
      error stop "grid.fem is required for sysGauge.fem"
    end if

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("sysGauge.fem.length")) then
      call SysGauge_Fem_Length_Fabricate

    else
      error stop "sysGauge.fem is missing one of: length"
    end if

  end subroutine

end submodule
