! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_GroundSolver) S_GroundSolver

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine GroundSolver_Fabricate()
    use M_Utils_Json
    use M_Utils_Say
    use M_GroundSolver_Tdhx

    call Say_Fabricate("groundSolver")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("groundSolver.tdhx")) then
      call GroundSolver_Tdhx_Fabricate

    else
      error stop "groundSolver is missing one of: tdhx"
    end if

  end subroutine

end submodule
