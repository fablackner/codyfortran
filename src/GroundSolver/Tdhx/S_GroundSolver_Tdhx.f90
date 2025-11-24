! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_GroundSolver_Tdhx) S_GroundSolver_Tdhx

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine GroundSolver_Tdhx_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_GroundSolver_Tdhx_StdImpl
    use M_GroundSolver_Tdhx_YlmOpt

    call Say_Fabricate("groundSolver.tdhx")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("groundSolver.tdhx.stdImpl")) then
      call GroundSolver_Tdhx_StdImpl_Fabricate

    else if (Json_GetExistence("groundSolver.tdhx.ylmOpt")) then
      call GroundSolver_Tdhx_YlmOpt_Fabricate

    else
      error stop "groundSolver.tdhx is missing one of: stdImpl, ylmOpt"

    end if

  end subroutine

end submodule
