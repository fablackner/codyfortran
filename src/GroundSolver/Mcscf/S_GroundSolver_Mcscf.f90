! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Dispatch submodule for MCSCF ground-solver implementation selection.
!>
!> @details
!> Reads `groundSolver.mcscf.*` and dispatches to:
!>   - `stdImpl`: Standard implementation (general grids)
submodule(M_GroundSolver_Mcscf) S_GroundSolver_Mcscf

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine GroundSolver_Mcscf_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_GroundSolver_Mcscf_StdImpl

    call Say_Fabricate("groundSolver.mcscf")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("groundSolver.mcscf.stdImpl")) then
      call GroundSolver_Mcscf_StdImpl_Fabricate

    else
      error stop "groundSolver.mcscf is missing one of: stdImpl"

    end if

  end subroutine

end submodule
