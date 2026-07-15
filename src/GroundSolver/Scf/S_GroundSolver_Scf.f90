! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Dispatch submodule for SCF implementation selection.
!>
!> @details
!> Reads `groundSolver.scf.*` and dispatches to:
!>   - `stdImpl`: Standard implementation (general grids)
!>   - `ylmOpt`: Ylm-optimized implementation (spherical grids)
submodule(M_GroundSolver_Scf) S_GroundSolver_Scf

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine GroundSolver_Scf_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_GroundSolver_Scf_StdImpl
    use M_GroundSolver_Scf_YlmOpt

    call Say_Fabricate("groundSolver.scf")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("groundSolver.scf.stdImpl")) then
      call GroundSolver_Scf_StdImpl_Fabricate

    else if (Json_GetExistence("groundSolver.scf.ylmOpt")) then
      call GroundSolver_Scf_YlmOpt_Fabricate

    else
      error stop "groundSolver.scf is missing one of: stdImpl, ylmOpt"

    end if

  end subroutine

end submodule
