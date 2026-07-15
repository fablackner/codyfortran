! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Standard SCF implementation submodule.
!>
!> @details
!> Implements the SCF iteration for general grids (non-Ylm). The algorithm:
!>
!>   1. **Setup**: Allocate working arrays for potentials and temporaries
!>   2. **Approach** (per iteration):
!>      a. Fill external potential V̂_ext(t)
!>      b. Build Hartree potential: Ĵ = 2·∑_j ∫|φ_j|²/|r-r'| (factor 2 for spin)
!>      c. Diagonalize Fock operator F̂ = T̂ + V̂_ext + Ĵ − K̂ via ARPACK
!>      d. Mix: orbs_new = (1-α)·orbs_old + α·evecs
!>      e. Copy spin-up → spin-down (restricted)
!>      f. Gram–Schmidt orthonormalize
!>
!> ## Working Arrays
!>
!> - `hartreePotential`: Precomputed direct (Hartree) potential on grid
!> - `externalPotential`: External potential V̂_ext (recomputed if time-dependent)
!> - `interactionPotential`: Temporary for exchange term computation
!> - `src`: Interaction source density ρ_ij(r) = φ_i*(r)·φ_j(r)
!> - `dOrbTmp`: Temporary for operator action results
!>
!> @note The exchange term K̂ is computed on-the-fly in `HartreeFockAction`
!>       since it depends on which orbital is being acted upon.
submodule(M_GroundSolver_Scf_StdImpl) S_GroundSolver_Scf_StdImpl

  implicit none

  !-----------------------------------------------------------------------------
  ! Module-level working arrays (allocated in Setup, used in Approach/HartreeFockAction)
  !-----------------------------------------------------------------------------

  !> Temporary for operator action results, dimension(Grid_nPoints)
  complex(R64), allocatable :: dOrbTmp(:)

  !> Precomputed Hartree (direct) potential on grid, dimension(potSize)
  complex(R64), allocatable :: hartreePotential(:)

  !> External potential V̂_ext, dimension(Grid_nPoints) - recomputed each Approach
  complex(R64), allocatable :: externalPotential(:)

  !> Temporary for exchange potential computation, dimension(potSize)
  complex(R64), allocatable :: interactionPotential(:)

  !> Interaction source density ρ_ij(r) = φ_i*(r)·φ_j(r), dimension(potSize)
  complex(R64), allocatable :: src(:)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine GroundSolver_Scf_StdImpl_Fabricate
    use M_Utils_Say
    use M_Utils_Json
    use M_GroundSolver
    use M_GroundSolver_Scf

    call Say_Fabricate("groundSolver.scf.stdImpl")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    GroundSolver_Setup => Setup
    GroundSolver_Approach => Approach
    GroundSolver_Scf_HartreeFockAction => HartreeFockAction

  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Allocate working arrays for the standard SCF implementation.
  !>
  !> @details
  !> Allocates potentials and temporaries based on the configured grid size.
  !> For Ylm grids, the potential size is (2·lmax_pot+1)² × nRadial to accommodate
  !> all angular momentum channels of the interaction.
  subroutine Setup
    use M_Utils_Say
    use M_Grid
    use M_Grid_Ylm
    use M_SysInteraction_Ylm

    integer(I32) :: lmaxPot, potSize

    call Say_Setup("groundSolver.scf.stdImpl")

    lmaxPot = SysInteraction_Ylm_lmax
    potSize = (2 * lmaxPot + 1)**2 * Grid_Ylm_nRadial

    allocate (hartreePotential(potSize))
    allocate (dOrbTmp(Grid_nPoints))

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Apply the SCF Fock operator to an orbital: dOrb = F̂·orb.
  !>
  !> @details
  !> Computes dOrb = (T̂ + V̂_ext + Ĵ − K̂)·orb where:
  !>   - T̂·orb: Kinetic energy via `SysKinetic_MultiplyWithKineticOp`
  !>   - V̂_ext·orb: External potential (precomputed in Approach)
  !>   - Ĵ·orb: Hartree potential (precomputed in Approach)
  !>   - K̂·orb: Exchange, computed on-the-fly as
  !>            K̂·φ = ∑_j W(φ_j, φ)·φ_j where W(a,b)(r) = ∫a*(r')b(r')/|r-r'|dr'
  !>
  !> The exchange sum runs over occupied orbitals (j = 1 to nOrbs/2 for restricted).
  subroutine HartreeFockAction(dOrb, orb, time)
    use M_SysKinetic
    use M_SysPotential
    use M_SysInteraction
    use M_Method_Mb_OrbBased
    use M_Orbs

    complex(R64), intent(out), contiguous, target :: dOrb(:)
    complex(R64), intent(in), contiguous, target :: orb(:)
    real(R64), intent(in) :: time

    integer(I32) :: j

    dOrb = 0.0_R64

    ! Kinetic energy: T̂·orb
    call SysKinetic_MultiplyWithKineticOp(dOrbTmp, orb, time)
    dOrb = dOrb + dOrbTmp

    ! External potential: V̂_ext·orb (externalPotential filled in Approach)
    call SysPotential_MultiplyWithExternalPotential(dOrbTmp, externalPotential, orb)
    dOrb = dOrb + dOrbTmp

    ! Hartree (direct) term: Ĵ·orb (hartreePotential filled in Approach)
    call SysInteraction_MultiplyWithInteractionPotential(dOrbTmp, hartreePotential, orb)
    dOrb = dOrb + dOrbTmp

    ! Exchange term: −K̂·orb = −∑_j W(φ_j, orb)·φ_j
    ! Loop over occupied orbitals (spin-up only for restricted; factor 2 in Hartree)
    do j = 1, Orbs_nOrbsInState / 2
      call SysInteraction_FillInteractionSrc(src, Orbs_orbs(:, j), orb(:))
      call SysInteraction_FillInteractionPotential(interactionPotential, src, time)
      call SysInteraction_MultiplyWithInteractionPotential(dOrbTmp, interactionPotential, Orbs_orbs(:, j))
      dOrb = dOrb - dOrbTmp
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Perform one SCF iteration: build Fock, diagonalize, mix, orthonormalize.
  !>
  !> @details
  !> Algorithm:
  !>   1. Fill external potential V̂_ext for current time
  !>   2. Build Hartree potential: Ĵ = 2·∑_j ∫|φ_j|²/|r-r'| (factor 2 for spin degeneracy)
  !>   3. Diagonalize Fock operator via DiagonalizerList(1) (uses HartreeFockAction callback)
  !>   4. Mix eigenvectors with old orbitals: φ_new = (1-α)·φ_old + α·evec
  !>   5. Copy spin-up orbitals to spin-down (restricted calculation)
  !>   6. Gram–Schmidt orthonormalize the mixed orbitals
  !>
  !> @param[inout] state  Packed orbital state, dimension(Grid_nPoints × nOrbsInState)
  !> @param[in]    alpha  Mixing parameter (0 < α ≤ 1)
  !> @param[in]    time   Time for potential evaluation (typically 0)
  subroutine Approach(state, alpha, time)
    use M_Grid
    use M_Grid_Ylm
    use M_DiagonalizerList
    use M_Orbs
    use M_SysPotential
    use M_SysInteraction
    use M_Method_Mb_OrbBased

    complex(R64), intent(inout), contiguous, target :: state(:)
    real(R64), intent(in) :: alpha
    real(R64), intent(in) :: time

    complex(R64), contiguous, pointer :: orbs(:, :)
    integer(I32) :: nG, nOS, j

    nG = Grid_nPoints
    nOS = Orbs_nOrbsInState
    orbs(1:nG, 1:nOS) => state(1:)

    ! Step 1: Fill external potential for current time
    call SysPotential_FillExternalPotential(externalPotential, time)

    ! Step 2: Build Hartree potential from occupied orbitals
    ! Factor 2 accounts for spin degeneracy (restricted: spin-up = spin-down)
    hartreePotential = 0.0_R64
    do j = 1, nOS / 2
      call SysInteraction_FillInteractionSrc(src, orbs(:, j), orbs(:, j))
      call SysInteraction_FillInteractionPotential(interactionPotential, 2.0_R64 * src, time)
      hartreePotential = hartreePotential + interactionPotential
    end do

    ! Step 3: Diagonalize Fock operator (HartreeFockAction is the matvec callback)
    call DiagonalizerList(1) % e % Diagonalize(time, .true.)

    ! Step 4-5: Mix eigenvectors with old orbitals and copy spin-up → spin-down
    do j = 1, nOS / 2
      orbs(:, j) = (1.0_R64 - alpha) * orbs(:, j) + alpha * DiagonalizerList(1) % e % evecs(:, j)
      orbs(:, nOS / 2 + j) = orbs(:, j)
    end do

    ! Step 6: Gram–Schmidt orthonormalize
    call Orbs_Orthonormalize(orbs)

  end subroutine

end submodule
