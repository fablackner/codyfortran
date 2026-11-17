! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief MCSCF by a preconditioned orbital-gradient iteration (implementation).
!>
!> @details
!> See M_GroundSolver_Mcscf_Gradient. The packed state layout matches the
!> MCTDHX method: state(1:nC) = CI coefficients, state(nC+1:) = orbitals.
!> Two pair-potential passes per iteration (h2 of the current orbitals, then
!> the correlation operator with the 2-RDM of their CI eigenvector): with the
!> CI vector of the previous orbitals in the residual (one pass) the iteration
!> diverges for weakly occupied orbitals (src/GroundSolver/A_GroundSolver.md).
submodule(M_GroundSolver_Mcscf_Gradient) S_GroundSolver_Mcscf_Gradient

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine GroundSolver_Mcscf_Gradient_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_GroundSolver

    call Say_Fabricate("groundSolver.mcscf.gradient")

    if (.not. Json_GetExistence("method.mb.orbBased.mctdhx")) &
      error stop "groundSolver.mcscf.gradient requires method.mb.orbBased.mctdhx (packed CI ⊕ orbitals state)"

    GroundSolver_Approach => Approach
    GroundSolver_implementsDeclaredMQ = .true.

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief One MCSCF iteration: CI step, residual pass, preconditioned
  !> orbital step (see the module header).
  subroutine Approach(state, time)
    use M_Grid
    use M_Orbs
    use M_Coeffs
    use M_Utils_LapackLib
    use M_Mixing
    use M_GroundSolver
    use M_Method_Mb_OrbBased

    complex(R64), intent(inout), contiguous, target :: state(:)
    real(R64), intent(in) :: time

    complex(R64), contiguous, pointer :: coeffs(:), orbs(:, :), orbsRawFlat(:)
    complex(R64), allocatable, target :: orbsRaw(:, :)
    complex(R64), allocatable :: fOrbs(:, :), natOrbs(:, :), natRot(:, :)
    complex(R64), allocatable :: rdm1(:, :), rdm2(:, :, :, :), invRdm1(:, :), rdm1Spatial(:, :)
    real(R64), allocatable :: lambda(:), natOcc(:)
    integer(I32) :: nC, nG, nOS, j

    nC = Coeffs_nCoeffs
    nG = Grid_nPoints
    nOS = Orbs_nOrbsInState

    coeffs(1:nC) => state(1:nC)
    orbs(1:nG, 1:nOS) => state(nC + 1:)

    !------------------------------------
    ! CI step: the CI eigenvector of the current orbitals (also at the first
    ! call, so that the orbital step sees the correlation of the start
    ! orbitals: with a unary start vector the unoccupied orbitals would get a
    ! bare-h step) and its energy
    !------------------------------------

    call Method_Mb_OrbBased_FillH1(Method_Mb_OrbBased_h1, orbs, time)
    call Method_Mb_OrbBased_FillH2(Method_Mb_OrbBased_h2, orbs, time)
    call CiStep(coeffs, time)

    !------------------------------------
    ! residual pass: F φ = h φ + ρ⁻¹ Γ W φ with the RDMs of that CI vector
    !------------------------------------

    call Coeffs_ApplyH1FillRdm1(coeffs, rdm1_=rdm1)
    call Coeffs_ApplyH2FillRdm2(coeffs, rdm2_=rdm2)
    call Method_Mb_OrbBased_FillInverseRdm1(invRdm1, rdm1)

    allocate (fOrbs(nG, nOS))
    fOrbs(:, :) = 0.0_R64
    call Method_Mb_OrbBased_ApplySingleBodyOp(fOrbs, orbs, time)
    call Method_Mb_OrbBased_ApplyCorrelationOp(fOrbs, orbs, invRdm1, rdm2, time)

    !------------------------------------
    ! natural-orbital frame of the step: rotations among orbitals of one
    ! symmetry block (e.g. 2p and 3p of equal m) leave the state unchanged and
    ! are not driven, so the stored orbitals may be any mixture of natural
    ! orbitals of very different occupation, and the shift λ_j of one
    ! preconditioner per orbital is only meaningful for a natural orbital.
    ! The orbital equations are covariant under such rotations U: the F φ of
    ! the orbitals φ U is (F φ) U.
    !------------------------------------

    call Method_Mb_OrbBased_FillSpatialRdm1(rdm1Spatial, rdm1)
    call LapackLib_DiagonalizeHermitianBlocks(natOcc, natRot, rdm1Spatial, Orbs_blockOfOrb)
    natOrbs = matmul(orbs, natRot)
    fOrbs = matmul(fOrbs, natRot)

    allocate (lambda(nOS))
    do j = 1, nOS
      lambda(j) = real(Grid_InnerProduct(natOrbs(:, j), fOrbs(:, j)), kind=R64)
    end do

    ! Residual g = (1 − P) F φ of the natural orbitals and the norm of the
    ! energy gradient g ρ, (Σ_j n_j² ‖g_j‖²)^½ in this frame
    call Orbs_ProjectOnSubspace(fOrbs, orbs)
    GroundSolver_Mcscf_Gradient_residualNorm = 0.0_R64
    do j = 1, nOS
      GroundSolver_Mcscf_Gradient_residualNorm = GroundSolver_Mcscf_Gradient_residualNorm + &
                                                 natOcc(j)**2 * real(Grid_InnerProduct(fOrbs(:, j), fOrbs(:, j)), kind=R64)
    end do
    GroundSolver_Mcscf_Gradient_residualNorm = sqrt(GroundSolver_Mcscf_Gradient_residualNorm)

    !------------------------------------
    ! preconditioned orbital step, rotated back from the natural-orbital
    ! frame; orbsRaw holds the step, then the stepped orbitals
    !------------------------------------

    allocate (orbsRaw(nG, nOS))
    do j = 1, nOS
      call GroundSolver_PreconditionAction(orbsRaw(:, j), fOrbs(:, j), lambda(j), time)
    end do
    orbsRaw = matmul(orbsRaw, transpose(conjg(natRot)))
    call Orbs_ProjectOnSubspace(orbsRaw, orbs)
    orbsRaw(:, :) = orbs(:, :) - orbsRaw(:, :)

    !------------------------------------
    ! mix and orthonormalize the orbitals
    !------------------------------------

    orbsRawFlat(1:nG * nOS) => orbsRaw
    call Mixing_Mix(state(nC + 1:nC + nG * nOS), orbsRawFlat)
    call Orbs_Orthonormalize(orbs)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief coeffs ← lowest eigenvector of H[h1, h2] (DiagonalizerList(1),
  !> seeded with coeffs), projected on a declared total-M sector (the
  !> diagonalizer works on full-length vectors, whose eigenvectors lie in the
  !> sector only up to round-off) and normalized;
  !> GroundSolver_Mcscf_Gradient_energy ← its eigenvalue.
  subroutine CiStep(coeffs, time)
    use M_Coeffs
    use M_DiagonalizerList
    use M_GroundSolver

    complex(R64), intent(inout), contiguous :: coeffs(:)
    real(R64), intent(in) :: time

    call DiagonalizerList(1) % e % SetStartVectors(reshape(coeffs, [size(coeffs), 1]))
    call DiagonalizerList(1) % e % Diagonalize(time, .true.)
    call GroundSolver_NoteUnconverged(1)

    coeffs(:) = DiagonalizerList(1) % e % evecs(:, 1)
    call Coeffs_ProjectOnSector(coeffs)
    call Coeffs_Normalize(coeffs)
    GroundSolver_Mcscf_Gradient_energy = DiagonalizerList(1) % e % evals(1)

  end subroutine

end submodule
