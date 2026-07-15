! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Implementation submodule for prolate orbital initialization dispatch.
submodule(M_OrbsInit_Prolate) S_OrbsInit_Prolate

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine OrbsInit_Prolate_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_OrbsInit_Prolate_Lcao

    call Say_Fabricate("orbsInit.prolate")

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("orbsInit.prolate.lcao")) then
      call OrbsInit_Prolate_Lcao_Fabricate

    else
      error stop "orbsInit.prolate is missing one of: lcao"
    end if

  end subroutine

end submodule
