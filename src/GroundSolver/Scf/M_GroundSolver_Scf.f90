! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief SCF (Time-Dependent Hartree–Exchange) backend interface module.
!>
!> @details
!> Provides SCF-specific procedure pointers for the ground-state solver. SCF
!> is a mean-field method that constructs an effective single-particle Fock
!> operator from the Hartree (direct) and Fock (exchange) potentials.
!>
!> This module declares the interface; concrete implementations (StdImpl, YlmOpt)
!> bind the pointers via their `_Fabricate` routines.
!>
!> ## SCF Iteration (SCF)
!>
!> Each iteration builds and diagonalizes the Fock operator:
!>
!>   F̂ = T̂ + V̂_ext + Ĵ − K̂
!>
!> where:
!>   - T̂ = kinetic energy operator
!>   - V̂_ext = external potential (e.g., nuclear Coulomb)
!>   - Ĵ = Hartree (direct) potential: Ĵφ_i = (∑_j ∫|φ_j|²/|r-r'|) φ_i
!>   - K̂ = Fock (exchange) operator: K̂φ_i = ∑_j (∫φ_j*(r')φ_i(r')/|r-r'|) φ_j(r)
!>
!> ## Available Implementations
!>
!> - **stdImpl**: General grid; diagonalizes full nPoints×nPoints Fock matrix
!> - **ylmOpt**: Spherical harmonics; exploits angular symmetry to diagonalize
!>               smaller nRadial×nRadial matrices per l-channel
!>
!> @see M_GroundSolver, M_DiagonalizerList
module M_GroundSolver_Scf
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> @brief Bind SCF-specific callbacks for the ground-state solver.
    !>
    !> @details
    !> Reads `groundSolver.scf.*` from JSON and dispatches to either:
    !>   - `stdImpl`: Standard implementation for general grids
    !>   - `ylmOpt`: Optimized implementation for Ylm (spherical) grids
    !>
    !> Assigns `GroundSolver_Scf_HartreeFockAction` to the selected backend's
    !> Fock operator application routine.
    !>
    !> @pre JSON configuration loaded; `groundSolver.scf` key exists
    !> @post `GroundSolver_Scf_HartreeFockAction` is bound
    module subroutine GroundSolver_Scf_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> @brief Pointer to the SCF Fock operator action (matrix-vector product).
  !>
  !> This callback computes F̂·φ for the iterative diagonalizer (e.g., ARPACK).
  procedure(I_GroundSolver_Scf_HartreeFockAction), pointer :: GroundSolver_Scf_HartreeFockAction
  abstract interface
    !> @brief Apply the SCF Fock operator to an orbital.
    !>
    !> @details
    !> Computes dOrb = F̂ · orb = (T̂ + V̂_ext + Ĵ − K̂) · orb
    !>
    !> This is used as the `ApplyMatOnVec` callback for iterative eigensolvers
    !> (DiagonalizerList). The Hartree potential Ĵ is precomputed and stored;
    !> the exchange term K̂ is computed on-the-fly for each orbital.
    !>
    !> @param[out] dOrb  Result of F̂ applied to orb (same dimension as orb)
    !> @param[in]  orb   Input orbital/trial vector
    !> @param[in]  time  Time parameter (typically 0 for ground state)
    subroutine I_GroundSolver_Scf_HartreeFockAction(dOrb, orb, time)
      import :: R64
      !> Output: Fock operator applied to the orbital, F̂·orb
      complex(R64), intent(out), contiguous, target :: dOrb(:)
      !> Input orbital/trial vector
      complex(R64), intent(in), contiguous, target :: orb(:)
      !> Time parameter (passed through to potential routines)
      real(R64), intent(in) :: time
    end subroutine
  end interface

end module

