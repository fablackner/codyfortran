! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief MCSCF by a preconditioned orbital-gradient iteration (interface).
!>
!> @details
!> Variational MCSCF on the packed MCTDHX state (CI coefficients ⊕
!> orbitals), any grid. Unlike stdImpl/ylmOpt, whose orbital step
!> diagonalizes a Fock operator of the 1-RDM only, the orbital update is
!> driven by the residual of the MCTDHX orbital equations,
!>
!>   g_j = (1 − P) [ h φ_j + Σ_k (ρ⁻¹ Γ W φ)_j ]
!>
!> (ρ⁻¹ the regularized inverse 1-RDM, Γ the 2-RDM), so the fixed point is
!> the stationary point of imaginary-time MCTDHX: the same ground state,
!> without the time-step bias. Per `GroundSolver_Approach`:
!>
!>   1. CI step: c ← lowest eigenvector of H[h1, h2] of the current orbitals
!>      (DiagonalizerList(1), seeded with c, projected on a declared total-M
!>      sector, normalized); its eigenvalue → GroundSolver_Mcscf_Gradient_energy;
!>   2. residual pass with the 1- and 2-RDM of that c;
!>   3. orbital step in the natural-orbital frame of each symmetry block
!>      (Orbs_blockOfOrb), φ̃_j ← φ̃_j − (1 − P) K_j g̃_j with
!>      K_j = `GroundSolver_PreconditionAction`, shift λ_j = Re⟨φ̃_j|F φ̃_j⟩,
!>      rotated back (the stored orbitals keep their gauge), through the
!>      Mixing module (linear or DIIS), then Orbs_Orthonormalize.
!>
!> Why steps 1–3 take this form (weakly occupied orbitals), benchmarks and
!> the configuration: src/GroundSolver/A_GroundSolver.md. Without
!> `SysKinetic_SolveShifted` the preconditioner is the identity (slow). The
!> orbitals keep a declared m (method.mb.orbBased.mOfOrbs): every operator
!> above conserves m.
!>
!> ## JSON Configuration
!>
!> ```json
!> {
!>     "groundSolver": { "mcscf": { "gradient": { } } },
!>     "diagonalizerList": { "lapackCi": { "nEvals": 1 } },
!>     "mixing": { "diis": { "nHistory": 8, "startThreshold": 0.1 } }
!> }
!> ```
!>
!> Requires method.mb.orbBased.mctdhx and DiagonalizerList(1) over the CI
!> space (`dim = Coeffs_nCoeffs`, callback
!> `GroundSolver_Mcscf_HamiltonianAction`); no orbital diagonalizer.
!>
!> @see M_GroundSolver_Mcscf, M_Method_Mb_OrbBased_Mctdhx
module M_GroundSolver_Mcscf_Gradient
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> @brief Bind `GroundSolver_Approach` to the gradient iteration.
    !>
    !> @pre JSON key `groundSolver.mcscf.gradient` exists
    module subroutine GroundSolver_Mcscf_Gradient_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Norm of the orbital energy gradient of the last `GroundSolver_Approach`,
  !> (Σ_j n_j² ‖g̃_j‖²)^½ in the natural-orbital frame of its CI vector (zero
  !> at the MCSCF stationary point), for convergence tests.
  real(R64) :: GroundSolver_Mcscf_Gradient_residualNorm = huge(1.0_R64)

  !> Energy of the last `GroundSolver_Approach`: the CI eigenvalue of the
  !> orbitals its step started from, i.e. the MCSCF energy of the orbitals
  !> that GroundSolver_Mcscf_Gradient_residualNorm refers to (the state after
  !> the call holds the stepped orbitals and that CI vector). Saves the pass of
  !> Method_GetEnergy per iteration.
  real(R64) :: GroundSolver_Mcscf_Gradient_energy = huge(1.0_R64)

end module
