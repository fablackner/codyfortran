! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysGauge_Prolate) S_SysGauge_Prolate

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysGauge_Prolate_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Grid_Prolate
    use M_SysGauge_Length

    call Say_Fabricate("sysGauge.prolate")

    if (.not. Json_GetExistence("grid.prolate")) then
      error stop "grid.prolate is required for sysGauge.prolate"
    end if

    if (Json_GetExistence("sysGauge.prolate.length")) then
      call SysGauge_Length_Fabricate("sysGauge.prolate.length", Grid_Prolate_FillNodeCoordinate)
      if (any(abs(SysGauge_Length_polarization(1:2)) > 0.0_R64)) &
        error stop "sysGauge.prolate.length.polarization must be along z (x and y couple the m channels)"
    else
      error stop "sysGauge.prolate is missing one of: length"
    end if

  end subroutine

end submodule
