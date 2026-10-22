! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Dispatch submodule for ground-state solver backend selection.
!>
!> @details
!> Reads the JSON configuration under `groundSolver.*` and dispatches to the
!> appropriate backend. Currently supported:
!>
!>   - `groundSolver.scf`: Self-consistent field (Hartree–Fock) method
!>   - `groundSolver.mcscf`: Multiconfiguration self-consistent field (MCSCF) method
!>
!> This is the only file that needs modification when adding new top-level
!> ground-state solver backends.
submodule(M_GroundSolver) S_GroundSolver

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine GroundSolver_Fabricate()
    use M_Utils_Json
    use M_Utils_Say
    use M_GroundSolver_Scf
    use M_GroundSolver_Mcscf
    use M_Method_Mb_OrbBased

    call Say_Fabricate("groundSolver")

    ! Declared m (method.mb.orbBased.mOfOrbs): only the Ylm-optimized solvers
    ! keep the orbitals in their m channels (radial Fock operators per l,
    ! alignment per m block); the standard ones diagonalize on the full grid
    ! and mix degenerate orbitals of different m
    if (Method_Mb_OrbBased_definiteMQ) then
      if (.not. (Json_GetExistence("groundSolver.scf.ylmOpt") .or. Json_GetExistence("groundSolver.mcscf.ylmOpt"))) &
        error stop "method.mb.orbBased.mOfOrbs requires groundSolver.scf.ylmOpt or groundSolver.mcscf.ylmOpt"
    end if

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("groundSolver.scf")) then
      call GroundSolver_Scf_Fabricate

    else if (Json_GetExistence("groundSolver.mcscf")) then
      call GroundSolver_Mcscf_Fabricate

    else
      error stop "groundSolver is missing one of: scf, mcscf"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module function GroundSolver_NumberOfSpinUpOrbs() result(nUp)
    use M_Orbs

    integer(I32) :: nUp

    if (Orbs_restrictedQ) then
      nUp = Orbs_nOrbsInState
    else
      nUp = Orbs_nOrbsInState / 2
    end if

  end function

end submodule
