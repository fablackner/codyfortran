! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief MCSCF (multiconfiguration self-consistent field) solver interface module.
!>
!> @details
!> Provides MCSCF-specific procedure pointers for the ground-state solver.
!> The MCSCF ground state is found by alternating two diagonalizations:
!>
!>   1. **CI step**: Diagonalize the CI Hamiltonian in the current orbital
!>      basis (coefficient update).
!>   2. **Orbital step**: Diagonalize an effective one-body Fock operator
!>      built from the correlated 1-RDM (orbital update), analogous to the
!>      SCF solver but with RDM-weighted mean-field and exchange terms.
!>
!> This module declares the interface; concrete implementations (StdImpl)
!> bind the pointers via their `_Fabricate` routines.
!>
!> ## MCSCF Iteration
!>
!> With the packed MCTDHX-style state |Ψ⟩ = (c, {φ_j}):
!>
!>   c    ← ground eigenvector of H_{IJ} = ⟨Φ_I|Ĥ[{φ}]|Φ_J⟩
!>   φ_j  ← lowest eigenvectors of F̂[ρ¹] = T̂ + V̂_ext + Ĵ[ρ¹] − K̂[ρ¹]
!>
!> where ρ¹ is the 1-RDM of the updated CI vector. For a single determinant
!> (ρ¹ = identity on the occupied space) this reduces to the Hartree–Fock
!> SCF iteration.
!>
!> ## Required Diagonalizers
!>
!> The program must fabricate two DiagonalizerList entries:
!>   - DiagonalizerList(1): CI space, `dim = Coeffs_nCoeffs`, callback wrapping
!>     `GroundSolver_Mcscf_HamiltonianAction`
!>   - DiagonalizerList(2): orbital space, `dim = Grid_nPoints`, callback
!>     wrapping `GroundSolver_Mcscf_FockAction`
!>
!> ## Available Implementations
!>
!> - **stdImpl**: General grids; full CI diagonalizer plus one orbital-space
!>   Fock diagonalizer
!>
!> @see M_GroundSolver, M_GroundSolver_Scf, M_DiagonalizerList
module M_GroundSolver_Mcscf
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> @brief Bind MCSCF-specific callbacks for the ground-state solver.
    !>
    !> @details
    !> Reads `groundSolver.mcscf.*` from JSON and dispatches to:
    !>   - `stdImpl`: Standard implementation for general grids
    !>
    !> Assigns `GroundSolver_Mcscf_HamiltonianAction` (CI matvec) and
    !> `GroundSolver_Mcscf_FockAction` (orbital matvec) to the selected
    !> backend's routines.
    !>
    !> @pre JSON configuration loaded; `groundSolver.mcscf` key exists
    !> @post `GroundSolver_Mcscf_HamiltonianAction` and
    !>       `GroundSolver_Mcscf_FockAction` are bound
    module subroutine GroundSolver_Mcscf_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> @brief Pointer to the CI Hamiltonian action (matrix-vector product).
  !>
  !> This callback computes Ĥ·c for the CI diagonalizer (DiagonalizerList(1)).
  procedure(I_GroundSolver_Mcscf_HamiltonianAction), pointer :: GroundSolver_Mcscf_HamiltonianAction
  abstract interface
    !> @brief Apply the CI Hamiltonian to a coefficient vector.
    !>
    !> @details
    !> Computes dCoeffs = Ĥ·coeffs where Ĥ is represented by the one- and
    !> two-body Hamiltonian matrices h1/h2 built from the current orbitals
    !> (precomputed and stored by the backend's `Approach` routine).
    !>
    !> @param[out] dCoeffs  Result of Ĥ applied to coeffs (same dimension)
    !> @param[in]  coeffs   Input CI coefficient/trial vector
    !> @param[in]  time     Time parameter (typically 0 for ground state)
    subroutine I_GroundSolver_Mcscf_HamiltonianAction(dCoeffs, coeffs, time)
      import :: R64
      !> Output: CI Hamiltonian applied to the coefficients, Ĥ·coeffs
      complex(R64), intent(out), contiguous, target :: dCoeffs(:)
      !> Input CI coefficient/trial vector
      complex(R64), intent(in), contiguous, target :: coeffs(:)
      !> Time parameter (passed through to potential routines)
      real(R64), intent(in) :: time
    end subroutine
  end interface

  !> @brief Pointer to the RDM-weighted Fock operator action (matrix-vector product).
  !>
  !> This callback computes F̂·φ for the orbital diagonalizer (DiagonalizerList(2)).
  procedure(I_GroundSolver_Mcscf_FockAction), pointer :: GroundSolver_Mcscf_FockAction
  abstract interface
    !> @brief Apply the RDM-weighted Fock operator to an orbital.
    !>
    !> @details
    !> Computes dOrb = F̂·orb = (T̂ + V̂_ext + Ĵ[ρ¹] − K̂[ρ¹])·orb where:
    !>   - Ĵ[ρ¹]: direct potential from the correlated density
    !>            ρ(r) = Σ_{pq} ρ¹_{pq} φ_p*(r)·φ_q(r) (precomputed and stored)
    !>   - K̂[ρ¹]: RDM-weighted exchange, computed on-the-fly as
    !>            K̂·φ = Σ_{pq} ρ¹_{pq}·W(φ_p, φ)·φ_q with
    !>            W(a,b)(r) = ∫a*(r')·b(r')·W(r,r')dr'
    !>
    !> The 1-RDM used here is filled from the updated CI vector by the
    !> backend's `Approach` routine before the diagonalization.
    !>
    !> @param[out] dOrb  Result of F̂ applied to orb (same dimension as orb)
    !> @param[in]  orb   Input orbital/trial vector
    !> @param[in]  time  Time parameter (typically 0 for ground state)
    subroutine I_GroundSolver_Mcscf_FockAction(dOrb, orb, time)
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
