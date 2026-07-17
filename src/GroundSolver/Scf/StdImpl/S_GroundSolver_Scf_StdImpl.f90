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
!>      b. Build raw Hartree potential: Ĵ = 2·∑_j ∫|φ_j|²/|r-r'| (factor 2 for
!>         spin); for `mixTarget="potential"` mix it into the accepted
!>         potential (Mixing_Mix), else use it directly
!>      c. Diagonalize Fock operator F̂ = T̂ + V̂_ext + Ĵ − K̂ via ARPACK
!>      d. Gauge-align eigenvectors with old orbitals (Orbs_AlignOnReference)
!>      e. Update orbitals: for `mixTarget="orbitals"` mix the aligned
!>         eigenvectors into the old orbitals (Mixing_Mix); for
!>         `mixTarget="potential"` replace them outright
!>      f. Copy spin-up → spin-down (restricted)
!>      g. Gram–Schmidt orthonormalize
!>
!> ## JSON Configuration
!>
!> @code{.json}
!> "groundSolver": { "scf": { "stdImpl": { "mixTarget": "orbitals" } } }
!> @endcode
!>
!> | Parameter   | Type   | Default    | Description                            |
!> |-------------|--------|------------|----------------------------------------|
!> | `mixTarget` | string | "orbitals" | What the mixer damps: "orbitals" mixes |
!> |             |        |            | gauge-aligned eigenvectors into the    |
!> |             |        |            | old orbitals; "potential" mixes the    |
!> |             |        |            | Hartree potential (the Poisson solve   |
!> |             |        |            | is linear in ρ) and replaces orbitals  |
!> |             |        |            | outright; "density" purifies the mixed |
!> |             |        |            | 1-RDM (1−α)·D_old + α·D_new to its top |
!> |             |        |            | natural orbitals (Orbs_MixOccupiedSpace|
!> |             |        |            | with α from mixing.linear.alpha)       |
!>
!> The mixing *algorithm* (linear, DIIS) is configured separately under the
!> top-level `mixing` block; see M_Mixing.
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

  !> Accepted Hartree (direct) potential on grid, dimension(potSize).
  !> For mixTarget="potential" this carries the mixed potential across
  !> iterations; for mixTarget="orbitals" it is rebuilt every Approach.
  complex(R64), allocatable :: hartreePotentialMixed(:)

  !> Raw Hartree potential built from the current orbitals, dimension(potSize)
  complex(R64), allocatable :: hartreePotentialRaw(:)

  !> External potential V̂_ext, dimension(Grid_nPoints) - recomputed each Approach
  complex(R64), allocatable :: externalPotential(:)

  !> Temporary for exchange potential computation, dimension(potSize)
  complex(R64), allocatable :: interactionPotential(:)

  !> Interaction source density ρ_ij(r) = φ_i*(r)·φ_j(r), dimension(potSize)
  complex(R64), allocatable :: src(:)

  !> Gauge-aligned new orbitals for mixing, dimension(Grid_nPoints, nOrbsInState/2)
  complex(R64), allocatable, target :: orbsRaw(:, :)

  !-----------------------------------------------------------------------------
  ! Mixing configuration
  !-----------------------------------------------------------------------------

  !> What quantity the mixer acts on: "orbitals" (gauge-aligned eigenvectors)
  !> or "potential" (Hartree potential; equivalent to density mixing since the
  !> interaction potential is linear in the source density) or "density"
  character(len=:), allocatable :: mixTarget

  !> Linear mixing alpha for density mixing.
  real(R64) :: alphaDensity = 1.0_R64

  !> Temporary array for density mixing.
  complex(R64), allocatable :: orbsTmp(:, :)
  real(R64), allocatable :: lambdaDiscarded(:)

  !> Whether hartreePotentialMixed already holds an accepted (mixed) potential
  logical :: hartreePotentialMixedInitializedQ = .false.

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine GroundSolver_Scf_StdImpl_Fabricate
    use M_Utils_Say
    use M_Utils_Json
    use M_GroundSolver
    use M_GroundSolver_Scf

    call Say_Fabricate("groundSolver.scf.stdImpl")

    !------------------------------------
    ! read configuration parameters
    !------------------------------------

    mixTarget = Json_Get("mixTarget", "orbitals", path_="groundSolver.scf.stdImpl")
    if (mixTarget .ne. "orbitals" .and. mixTarget .ne. "potential" .and. mixTarget .ne. "density") then
      error stop "groundSolver.scf.stdImpl.mixTarget must be one of: orbitals, potential, density"
    end if

    if (mixTarget .eq. "density") then
      alphaDensity = Json_Get("alpha", 1.0_R64, path_="mixing.linear")
    end if

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
    use M_GroundSolver
    use M_Orbs

    integer(I32) :: lmaxPot, potSize, nUp

    call Say_Setup("groundSolver.scf.stdImpl")

    lmaxPot = SysInteraction_Ylm_lmax
    potSize = (2 * lmaxPot + 1)**2 * Grid_Ylm_nRadial

    nUp = GroundSolver_NumberOfSpinUpOrbs()

    allocate (hartreePotentialMixed(potSize))
    allocate (hartreePotentialRaw(potSize))
    allocate (dOrbTmp(Grid_nPoints))
    allocate (orbsRaw(Grid_nPoints, nUp))
    allocate (orbsTmp(Grid_nPoints, nUp))
    allocate (lambdaDiscarded(nUp))

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
    use M_GroundSolver
    use M_Orbs

    complex(R64), intent(out), contiguous, target :: dOrb(:)
    complex(R64), intent(in), contiguous, target :: orb(:)
    real(R64), intent(in) :: time

    integer(I32) :: j, nUp

    nUp = GroundSolver_NumberOfSpinUpOrbs()

    dOrb = 0.0_R64

    ! Kinetic energy: T̂·orb
    call SysKinetic_MultiplyWithKineticOp(dOrbTmp, orb, time)
    dOrb = dOrb + dOrbTmp

    ! External potential: V̂_ext·orb (externalPotential filled in Approach)
    call SysPotential_MultiplyWithExternalPotential(dOrbTmp, externalPotential, orb)
    dOrb = dOrb + dOrbTmp

    ! Hartree (direct) term: Ĵ·orb (hartreePotentialMixed filled in Approach)
    call SysInteraction_MultiplyWithInteractionPotential(dOrbTmp, hartreePotentialMixed, orb)
    dOrb = dOrb + dOrbTmp

    ! Exchange term: −K̂·orb = −∑_j W(φ_j, orb)·φ_j
    ! Loop over occupied orbitals (spin-up only; factor 2 in Hartree)
    do j = 1, nUp
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
  !>   4. Gauge-align eigenvectors with old orbitals
  !>   5. Update orbitals per mixTarget: mix aligned eigenvectors ("orbitals")
  !>      or replace outright ("potential", damping applied to Ĵ in step 2)
  !>   6. Copy spin-up orbitals to spin-down (restricted calculation)
  !>   7. Gram–Schmidt orthonormalize the updated orbitals
  !>
  !> @param[inout] state  Packed orbital state, dimension(Grid_nPoints × nOrbsInState)
  !> @param[in]    time   Time for potential evaluation (typically 0)
  subroutine Approach(state, time)
    use M_Mixing
    use M_Grid
    use M_Grid_Ylm
    use M_DiagonalizerList
    use M_GroundSolver
    use M_Orbs
    use M_SysPotential
    use M_SysInteraction
    use M_Method_Mb_OrbBased

    complex(R64), intent(inout), contiguous, target :: state(:)
    real(R64), intent(in) :: time

    complex(R64), contiguous, pointer :: orbs(:, :)
    complex(R64), contiguous, pointer :: orbsRawFlat(:)
    integer(I32) :: nG, nOS, nUp, j

    nG = Grid_nPoints
    nOS = Orbs_nOrbsInState
    nUp = GroundSolver_NumberOfSpinUpOrbs()
    orbs(1:nG, 1:nOS) => state(1:)

    ! Step 1: Fill external potential for current time
    call SysPotential_FillExternalPotential(externalPotential, time)

    ! Step 2: Build raw Hartree potential from occupied orbitals
    ! Factor 2 accounts for spin degeneracy (restricted: spin-up = spin-down)
    hartreePotentialRaw = 0.0_R64
    do j = 1, nUp
      call SysInteraction_FillInteractionSrc(src, orbs(:, j), orbs(:, j))
      call SysInteraction_FillInteractionPotential(interactionPotential, 2.0_R64 * src, time)
      hartreePotentialRaw = hartreePotentialRaw + interactionPotential
    end do

    if (mixTarget .eq. "potential") then
      ! Damping acts on the (gauge-invariant) Hartree potential: mix the raw
      ! potential into the accepted one; the very first potential is taken as is
      if (hartreePotentialMixedInitializedQ) then
        call Mixing_Mix(hartreePotentialMixed, hartreePotentialRaw)
      else
        hartreePotentialMixed = hartreePotentialRaw
        hartreePotentialMixedInitializedQ = .true.
      end if
    else
      hartreePotentialMixed = hartreePotentialRaw
    end if

    ! Step 3: Diagonalize Fock operator (HartreeFockAction is the matvec callback)
    call DiagonalizerList(1) % e % Diagonalize(time, .true.)

    ! Step 4: Gauge-align the eigenvectors with the old orbitals (removes the
    ! arbitrary per-vector phase and the arbitrary rotation within degenerate
    ! shells)
    orbsRaw = DiagonalizerList(1) % e % evecs(:, 1:nUp)
    ! ARPACK normalizes in the Euclidean metric; alignment and mixing require
    ! orthonormality with respect to the grid metric (Grid_InnerProduct)
    call Grid_Orthonormalize(orbsRaw)
    call Orbs_AlignOnReference(orbsRaw, orbs(:, 1:nUp))

    ! Step 5: Update the spin-up orbitals and copy spin-up → spin-down
    if (mixTarget .eq. "orbitals") then
      orbsRawFlat(1:nG * nUp) => orbsRaw
      call Mixing_Mix(state(1:nG * nUp), orbsRawFlat)
    else if (mixTarget .eq. "density") then
      if (hartreePotentialMixedInitializedQ) then
        call Orbs_MixOccupiedSpace(orbs(:, 1:nUp), orbsRaw, alphaDensity, orbsTmp, lambdaDiscarded)
        orbs(:, 1:nUp) = orbsTmp
      else
        orbs(:, 1:nUp) = orbsRaw
        hartreePotentialMixedInitializedQ = .true.
      end if
    else
      ! mixTarget="potential": pure replacement, damping happened in Step 2
      do j = 1, nUp
        orbs(:, j) = orbsRaw(:, j)
      end do
    end if

    ! Restricted: the state holds only the shared spatial set, no copy needed
    if (.not. Orbs_restrictedQ) then
      do j = 1, nUp
        orbs(:, nUp + j) = orbs(:, j)
      end do
    end if

    ! Step 6: Gram–Schmidt orthonormalize
    call Orbs_Orthonormalize(orbs)

  end subroutine

end submodule
