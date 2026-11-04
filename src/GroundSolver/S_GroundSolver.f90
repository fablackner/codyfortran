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
!> For every backend it also binds the orbital preconditioner
!> `GroundSolver_PreconditionAction` (key `groundSolver.precShiftMin`), and
!> it provides `GroundSolver_NoteUnconverged`, the note every backend prints
!> after an unconverged diagonalizer solve.
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
    ! orbital preconditioner (all backends)
    !------------------------------------

    GroundSolver_precShiftMin = Json_Get("precShiftMin", 0.5_R64, path_="groundSolver")
    if (.not. (GroundSolver_precShiftMin > 0.0_R64)) then
      error stop "groundSolver.precShiftMin must be positive"
    end if

    GroundSolver_PreconditionAction => PreconditionAction

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
  module subroutine GroundSolver_NoteUnconverged(iDiag)
    use M_DiagonalizerList

    integer(I32), intent(in) :: iDiag

    if (DiagonalizerList(iDiag) % e % convergedQ) return

    write (*, '(1X, 3A, I0, A, I0, A)') "note: groundSolver: ", DiagonalizerList(iDiag) % e % path, &
      " did not converge (", DiagonalizerList(iDiag) % e % nFound, " of ", DiagonalizerList(iDiag) % e % nEvals, &
      " wanted eigenpairs); continuing with its best approximations"

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

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Orbital preconditioner: dOrb ≈ (F̂ − shift)⁻¹·orb.
  !>
  !> @details
  !> dOrb = (T̂ + max(−shift, GroundSolver_precShiftMin))⁻¹·orb through
  !> `SysKinetic_SolveShifted` (no body type: the mass of body type 1, i.e.
  !> of the spin-restricted orbitals) if the kinetic backend provides it,
  !> else dOrb = orb. The hook is looked up on every call, so the result
  !> does not depend on the order of the Fabricate calls. Bound to
  !> `GroundSolver_PreconditionAction` for every backend.
  subroutine PreconditionAction(dOrb, orb, shift, time)
    use M_SysKinetic
    use M_Utils_UnusedVariables

    complex(R64), intent(out), contiguous, target :: dOrb(:)
    complex(R64), intent(in), contiguous, target :: orb(:)
    real(R64), intent(in) :: shift
    real(R64), intent(in) :: time

    if (.false.) call UnusedVariables_Mark(time)

    if (associated(SysKinetic_SolveShifted)) then
      call SysKinetic_SolveShifted(dOrb, orb, cmplx(max(-shift, GroundSolver_precShiftMin), 0.0_R64, kind=R64))
    else
      dOrb = orb
    end if

  end subroutine

end submodule
