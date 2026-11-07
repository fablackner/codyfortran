! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Dispatch submodule for MCSCF ground-solver implementation selection.
!>
!> @details
!> Reads `groundSolver.mcscf.*` and dispatches to:
!>   - `stdImpl`: Standard implementation (general grids)
!>   - `ylmOpt`: Ylm-optimized implementation (spherical grids)
!>   - `gradient`: variational MCSCF by a preconditioned orbital-gradient
!>     iteration on the MCTDHX orbital residual (any grid)
!>
!> The CI matvec is common to all of them: the CI Hamiltonian of
!> Method_Mb_OrbBased_h1/h2, which every backend fills with the current
!> orbitals before its CI diagonalization.
submodule(M_GroundSolver_Mcscf) S_GroundSolver_Mcscf

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine GroundSolver_Mcscf_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Method_Mb_OrbBased
    use M_GroundSolver_Mcscf_StdImpl
    use M_GroundSolver_Mcscf_YlmOpt
    use M_GroundSolver_Mcscf_Gradient

    call Say_Fabricate("groundSolver.mcscf")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    GroundSolver_Mcscf_HamiltonianAction => Method_Mb_OrbBased_ApplyCiHamiltonian

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("groundSolver.mcscf.stdImpl")) then
      call GroundSolver_Mcscf_StdImpl_Fabricate

    else if (Json_GetExistence("groundSolver.mcscf.ylmOpt")) then
      call GroundSolver_Mcscf_YlmOpt_Fabricate

    else if (Json_GetExistence("groundSolver.mcscf.gradient")) then
      call GroundSolver_Mcscf_Gradient_Fabricate

    else
      error stop "groundSolver.mcscf is missing one of: stdImpl, ylmOpt, gradient"

    end if

  end subroutine

end submodule
