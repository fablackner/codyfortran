! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Ylm-optimized MCSCF implementation interface module.
!>
!> @details
!> Provides the factory routine and interface for a spherical-harmonics (Ylm)
!> optimized MCSCF ground-state solver. The CI step is identical to the
!> standard implementation; the orbital step exploits angular momentum
!> symmetry like the SCF YlmOpt backend:
!>
!> **Key Optimization**: Instead of diagonalizing the RDM-weighted Fock
!> operator F̂[ρ¹] on the full Grid_nPoints×Grid_nPoints space, this backend
!> diagonalizes (lmax+1) smaller nRadial×nRadial radial operators, one per
!> angular momentum channel l. This is valid for spherically symmetric states
!> (e.g., atomic S states) where orbitals factorize as
!> φ_{nlm}(r) = R_nl(r)·Y_lm(θ,φ) and the correlated density is monopole.
!>
!> ## Required Diagonalizers
!>
!> The program must fabricate (lmax+2) DiagonalizerList entries:
!>   - DiagonalizerList(1): CI space, `dim = Coeffs_nCoeffs`, callback wrapping
!>     `GroundSolver_Mcscf_HamiltonianAction`
!>   - DiagonalizerList(l+2), l = 0..lmax: radial orbital space,
!>     `dim = Grid_Ylm_nRadial`, l-specific wrapper callbacks around
!>     `GroundSolver_Mcscf_YlmOpt_FockAction`:
!>
!> ```fortran
!> DiagonalizerListInput(1) % ApplyMatOnVec => GroundSolver_Mcscf_HamiltonianAction
!> DiagonalizerListInput(2) % ApplyMatOnVec => ApplyMatOnVecL0  ! l=0
!> DiagonalizerListInput(3) % ApplyMatOnVec => ApplyMatOnVecL1  ! l=1
!> ! ...
!> ```
!>
!> Each l-channel's `nEvals` must cover the highest radial quantum number
!> nr = n - l - 1 occurring among the orbitals of that channel.
!>
!> @see M_GroundSolver_Mcscf_StdImpl for the general-grid alternative,
!>      M_GroundSolver_Scf_YlmOpt for the single-determinant analogue
module M_GroundSolver_Mcscf_YlmOpt
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> @brief Bind Ylm-optimized MCSCF callbacks for the ground-state solver.
    !>
    !> @details
    !> Assigns the following procedure pointers:
    !>   - `GroundSolver_Setup` → local `Setup`
    !>   - `GroundSolver_Approach` → local `Approach`
    !>   - `GroundSolver_Mcscf_HamiltonianAction` → local `HamiltonianAction`
    !>   - `GroundSolver_Mcscf_YlmOpt_FockAction` → local `FockAction`
    !>
    !> @pre JSON key `groundSolver.mcscf.ylmOpt` exists
    !> @pre Grid is configured as Ylm (spherical harmonics)
    module subroutine GroundSolver_Mcscf_YlmOpt_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> @brief Pointer to the Ylm radial RDM-weighted Fock action for a single l-channel.
  !>
  !> Used as the matvec callback for per-l diagonalizers. The caller must wrap
  !> this in l-specific callbacks (one per diagonalizer).
  procedure(I_GroundSolver_Mcscf_YlmOpt_FockRadialAction), pointer :: GroundSolver_Mcscf_YlmOpt_FockAction
  abstract interface
    !> @brief Apply the radial RDM-weighted Fock operator for channel l.
    !>
    !> @details
    !> Computes dOrbLm = F̂_l[ρ¹] · orbLm where F̂_l[ρ¹] is the RDM-weighted
    !> Fock operator projected onto channel l. This includes:
    !>   - Radial kinetic: T̂_l = -½d²/dr² + l(l+1)/(2r²)
    !>   - Radial external potential: V̂_ext(r)
    !>   - Radial mean-field potential Ĵ[ρ¹]: precomputed monopole (l=0)
    !>     projection of the correlated density's direct potential
    !>   - Exchange K̂[ρ¹]: computed on-the-fly with full angular coupling,
    !>     K̂·φ = Σ_{pq} ρ¹_{pq}·W(φ_p, φ)·φ_q
    !>
    !> The 1-RDM used here is filled from the updated CI vector by the
    !> backend's `Approach` routine before the diagonalization.
    !>
    !> @param[out] dOrbLm  Result of F̂_l applied to orbLm, dimension(nRadial)
    !> @param[in]  orbLm   Input radial orbital for channel l
    !> @param[in]  l       Angular momentum quantum number
    !> @param[in]  time    Time parameter (typically 0 for ground state)
    subroutine I_GroundSolver_Mcscf_YlmOpt_FockRadialAction(dOrbLm, orbLm, l, time)
      import :: I32, R64
      !> Output: radial Fock operator applied to the orbital, F̂_l·orbLm
      complex(R64), intent(out), contiguous, target :: dOrbLm(:)
      !> Input radial orbital for channel l
      complex(R64), intent(in), contiguous, target :: orbLm(:)
      !> Angular momentum quantum number (l = 0, 1, ..., lmax)
      integer(I32), intent(in) :: l
      !> Time parameter (passed through to potential routines)
      real(R64), intent(in) :: time
    end subroutine
  end interface

end module
