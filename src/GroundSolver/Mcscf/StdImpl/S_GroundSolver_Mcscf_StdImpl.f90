! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Standard MCSCF implementation submodule.
!>
!> @details
!> Implements the MCSCF iteration with two diagonalizations per step. The
!> packed state layout matches the MCTDHX method:
!> state(1:nC) = CI coefficients, state(nC+1:) = orbitals (column-major).
!>
!>   1. **Setup**: Allocate working arrays
!>   2. **Approach** (per iteration):
!>      a. Fill external potential V̂_ext(t) and build h1/h2 in the current
!>         orbital basis
!>      b. Diagonalize the CI Hamiltonian via DiagonalizerList(1)
!>         (HamiltonianAction is the matvec callback)
!>      c. Mix the phase-aligned CI ground eigenvector with the old
!>         coefficients: c_new = (1-α)·c_old + α·evec, then normalize
!>      d. Fill the 1-RDM from the updated coefficients and build the direct
!>         (Hartree) potential from the correlated density
!>         ρ(r) = Σ_{pq} ρ¹_{pq} φ_p*(r)·φ_q(r)
!>      e. Diagonalize the RDM-weighted Fock operator F̂[ρ¹] via
!>         DiagonalizerList(2) (FockAction is the matvec callback)
!>      f. Gauge-align the Fock eigenvectors with the old orbitals
!>         (Orbs_AlignOnReference), mix, copy spin-up → spin-down
!>         (restricted), Gram–Schmidt orthonormalize
!>
!> ## Working Arrays
!>
!> - `h1`, `h2`: One-/two-body Hamiltonian matrices in the current orbital
!>   basis; rebuilt each `Approach` and shared with `HamiltonianAction`
!> - `rdm1Up`: Spin-up block of the correlated 1-RDM, used by `FockAction`
!> - `hartreePotential`: Direct potential from the correlated density
!> - `externalPotential`: External potential V̂_ext (recomputed each Approach)
!> - `exchangePotentials`: Per-orbital exchange potentials W(φ_p, orb)
!>
!> @note The CI diagonalization requires DiagonalizerList(1) with
!>       `dim = Coeffs_nCoeffs`; the orbital diagonalization requires
!>       DiagonalizerList(2) with `dim = Grid_nPoints` and
!>       `nEvals = Orbs_nOrbsInState/2`.
!> @note Restricted calculations only: as in the SCF stdImpl, the spin-up
!>       orbitals are relaxed and copied to the spin-down block.
submodule(M_GroundSolver_Mcscf_StdImpl) S_GroundSolver_Mcscf_StdImpl

  implicit none

  !-----------------------------------------------------------------------------
  ! Module-level working arrays (allocated in Setup/Approach, used in callbacks)
  !-----------------------------------------------------------------------------

  !> One-body Hamiltonian matrix in the current orbital basis, filled in Approach
  complex(R64), allocatable :: h1(:, :)

  !> Two-body Hamiltonian matrix in the current orbital basis, filled in Approach
  complex(R64), allocatable :: h2(:, :, :, :)

  !> Spin-up block of the correlated 1-RDM, filled in Approach
  complex(R64), allocatable :: rdm1Up(:, :)

  !> Accepted direct (Hartree) potential from the correlated density
  complex(R64), allocatable :: hartreePotentialMixed(:)

  !> Raw direct potential from the correlated density
  complex(R64), allocatable :: hartreePotentialRaw(:)

  !> External potential V̂_ext, dimension(Grid_nPoints) - recomputed each Approach
  complex(R64), allocatable :: externalPotential(:)

  !> Per-orbital exchange potentials W(φ_p, orb), dimension(potSize, nOrbsInState/2)
  complex(R64), allocatable :: exchangePotentials(:, :)

  !> RDM-weighted sum of exchange potentials, dimension(potSize)
  complex(R64), allocatable :: weightedPotential(:)

  !> Accumulated RDM-weighted interaction source for the direct term
  complex(R64), allocatable :: weightedSrc(:)

  !> Temporary for pairwise interaction sources
  complex(R64), allocatable :: src(:)

  !> Temporary for exchange potential computation
  complex(R64), allocatable :: interactionPotential(:)

  !> Temporary for operator action results, dimension(Grid_nPoints)
  complex(R64), allocatable :: dOrbTmp(:)

  !> Gauge-aligned new orbitals for mixing, dimension(Grid_nPoints, nOrbsInState/2)
  complex(R64), allocatable, target :: orbsRaw(:, :)

  !> CI mixing parameter. The CI coefficients are mixed linearly with phase
  !> alignment: an eigenvector is not a fixed-point iterate in a linear space,
  !> so DIIS residual extrapolation does not apply to it cleanly.
  real(R64) :: alphaCi = 1.0_R64

  !-----------------------------------------------------------------------------
  ! Mixing configuration
  !-----------------------------------------------------------------------------

  !> What quantity the mixer acts on: "orbitals" or "potential"
  character(len=:), allocatable :: mixTarget

  !> Whether hartreePotentialMixed already holds an accepted (mixed) potential
  logical :: hartreePotentialMixedInitializedQ = .false.

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine GroundSolver_Mcscf_StdImpl_Fabricate
    use M_Utils_Say
    use M_Utils_Json
    use M_GroundSolver
    use M_GroundSolver_Mcscf

    call Say_Fabricate("groundSolver.mcscf.stdImpl")

    !------------------------------------
    ! read configuration parameters
    !------------------------------------

    alphaCi = Json_Get("alphaCi", 1.0_R64, path_="groundSolver.mcscf.stdImpl")

    mixTarget = Json_Get("mixTarget", "orbitals", path_="groundSolver.mcscf.stdImpl")
    if (mixTarget .ne. "orbitals" .and. mixTarget .ne. "potential") then
      error stop "groundSolver.mcscf.stdImpl.mixTarget must be one of: orbitals, potential"
    end if

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    GroundSolver_Setup => Setup
    GroundSolver_Approach => Approach
    GroundSolver_Mcscf_HamiltonianAction => HamiltonianAction
    GroundSolver_Mcscf_FockAction => FockAction

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Allocate working arrays for the standard MCSCF implementation.
  subroutine Setup
    use M_Utils_Say
    use M_Grid
    use M_Orbs

    call Say_Setup("groundSolver.mcscf.stdImpl")

    allocate (dOrbTmp(Grid_nPoints))
    allocate (rdm1Up(Orbs_nOrbsInState / 2, Orbs_nOrbsInState / 2))
    allocate (orbsRaw(Grid_nPoints, Orbs_nOrbsInState / 2))

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Apply the CI Hamiltonian to a coefficient vector: dCoeffs = Ĥ·coeffs.
  !>
  !> @details
  !> Uses the h1/h2 matrices prepared in `Approach` for the current orbitals.
  !> This is the matvec callback for DiagonalizerList(1) over the CI space.
  subroutine HamiltonianAction(dCoeffs, coeffs, time)
    use M_Coeffs

    complex(R64), intent(out), contiguous, target :: dCoeffs(:)
    complex(R64), intent(in), contiguous, target :: coeffs(:)
    real(R64), intent(in) :: time

    dCoeffs = 0.0_R64

    call Coeffs_ApplyH1FillRdm1(coeffs, dCoeffs_=dCoeffs, h1_=h1)
    call Coeffs_ApplyH2FillRdm2(coeffs, dCoeffs_=dCoeffs, h2_=h2)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Apply the RDM-weighted Fock operator to an orbital: dOrb = F̂[ρ¹]·orb.
  !>
  !> @details
  !> Computes dOrb = (T̂ + V̂_ext + Ĵ[ρ¹] − K̂[ρ¹])·orb where:
  !>   - T̂·orb: Kinetic energy via `SysKinetic_MultiplyWithKineticOp`
  !>   - V̂_ext·orb: External potential (precomputed in Approach)
  !>   - Ĵ[ρ¹]·orb: Direct potential from the correlated density
  !>     (precomputed in Approach)
  !>   - K̂[ρ¹]·orb: RDM-weighted exchange, computed on-the-fly as
  !>     K̂·φ = Σ_{pq} ρ¹_{pq}·W(φ_p, φ)·φ_q
  !>
  !> For ρ¹ = identity on the occupied space this reduces to the Hartree–Fock
  !> Fock operator of the SCF stdImpl.
  !>
  !> This is the matvec callback for DiagonalizerList(2) over the orbital space.
  subroutine FockAction(dOrb, orb, time)
    use M_SysKinetic
    use M_SysPotential
    use M_SysInteraction
    use M_Orbs

    complex(R64), intent(out), contiguous, target :: dOrb(:)
    complex(R64), intent(in), contiguous, target :: orb(:)
    real(R64), intent(in) :: time

    integer(I32) :: p, q, nUp

    nUp = Orbs_nOrbsInState / 2

    dOrb = 0.0_R64

    ! Kinetic energy: T̂·orb
    call SysKinetic_MultiplyWithKineticOp(dOrbTmp, orb, time)
    dOrb = dOrb + dOrbTmp

    ! External potential: V̂_ext·orb (externalPotential filled in Approach)
    call SysPotential_MultiplyWithExternalPotential(dOrbTmp, externalPotential, orb)
    dOrb = dOrb + dOrbTmp

    ! Direct (Hartree) term: Ĵ[ρ¹]·orb (hartreePotentialMixed filled in Approach)
    call SysInteraction_MultiplyWithInteractionPotential(dOrbTmp, hartreePotentialMixed, orb)
    dOrb = dOrb + dOrbTmp

    ! Exchange term: −K̂[ρ¹]·orb = −Σ_{pq} ρ¹_{pq}·W(φ_p, orb)·φ_q
    ! Spin-up block only (restricted; the trial vector carries no spin index)
    do p = 1, nUp
      call SysInteraction_FillInteractionSrc(src, Orbs_orbs(:, p), orb(:))
      call SysInteraction_FillInteractionPotential(interactionPotential, src, time)
      exchangePotentials(:, p) = interactionPotential
    end do

    do q = 1, nUp
      weightedPotential = 0.0_R64
      do p = 1, nUp
        weightedPotential = weightedPotential + rdm1Up(p, q) * exchangePotentials(:, p)
      end do
      call SysInteraction_MultiplyWithInteractionPotential(dOrbTmp, weightedPotential, Orbs_orbs(:, q))
      dOrb = dOrb - dOrbTmp
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Perform one MCSCF iteration: CI diagonalization + Fock diagonalization.
  !>
  !> @details
  !> Algorithm:
  !>   1. Fill external potential and h1/h2 in the current orbital basis
  !>   2. Diagonalize the CI Hamiltonian via DiagonalizerList(1)
  !>   3. Mix the phase-aligned CI ground eigenvector with the old coefficients
  !>      and normalize
  !>   4. Fill the 1-RDM from the updated coefficients and build the direct
  !>      potential from the correlated density
  !>   5. Diagonalize the RDM-weighted Fock operator via DiagonalizerList(2)
  !>   6. Gauge-align the Fock eigenvectors with the old orbitals, mix, copy
  !>      spin-up → spin-down, Gram–Schmidt orthonormalize
  !>
  !> @param[inout] state  Packed MCTDHX-style state: coefficients ⊕ orbitals
  !> @param[in]    time   Time for potential evaluation (typically 0)
  subroutine Approach(state, time)
    use M_Mixing
    use M_Grid
    use M_DiagonalizerList
    use M_Coeffs
    use M_Orbs
    use M_SysPotential
    use M_SysInteraction
    use M_Method_Mb_OrbBased

    complex(R64), intent(inout), contiguous, target :: state(:)
    real(R64), intent(in) :: time

    complex(R64), contiguous, pointer :: coeffs(:)
    complex(R64), contiguous, pointer :: orbs(:, :)
    complex(R64), contiguous, pointer :: orbsRawFlat(:)
    complex(R64), allocatable :: rdm1(:, :)
    complex(R64) :: overlap, phase
    integer(I32) :: nC, nG, nOS, p, q, j

    nC = Coeffs_nCoeffs
    nG = Grid_nPoints
    nOS = Orbs_nOrbsInState

    coeffs(1:nC) => state(1:nC)
    orbs(1:nG, 1:nOS) => state(nC + 1:)

    ! Step 1: External potential and Hamiltonian matrices in the current basis
    call SysPotential_FillExternalPotential(externalPotential, time)
    call Method_Mb_OrbBased_FillH1(h1, orbs, time)
    call Method_Mb_OrbBased_FillH2(h2, orbs, time)

    ! Step 2: Diagonalize CI Hamiltonian (HamiltonianAction is the matvec callback)
    call DiagonalizerList(1) % e % Diagonalize(time, .true.)

    ! Step 3: Mix phase-aligned CI ground eigenvector with old coefficients
    ! The eigenvector's global phase is arbitrary; align it so the mixing
    ! interpolates instead of cancelling.
    overlap = dot_product(DiagonalizerList(1) % e % evecs(:, 1), coeffs)
    phase = (1.0_R64, 0.0_R64)
    if (0.0_R64 < abs(overlap)) phase = overlap / abs(overlap)

    coeffs = (1.0_R64 - alphaCi) * coeffs + alphaCi * phase * DiagonalizerList(1) % e % evecs(:, 1)
    call Coeffs_Normalize(coeffs)

    ! Step 4: 1-RDM from updated coefficients; direct potential from the
    ! correlated density ρ(r) = Σ_{pq} ρ¹_{pq}·φ_p*(r)·φ_q(r) (all body types)
    call Coeffs_ApplyH1FillRdm1(coeffs, rdm1_=rdm1)
    rdm1Up = rdm1(1:nOS / 2, 1:nOS / 2)

    do p = 1, nOS
      do q = 1, nOS
        if (abs(rdm1(p, q)) .eq. 0.0_R64) cycle
        call SysInteraction_FillInteractionSrc(src, orbs(:, p), orbs(:, q))
        if (.not. allocated(weightedSrc)) then
          allocate (weightedSrc(size(src)))
          weightedSrc = 0.0_R64
        end if
        weightedSrc = weightedSrc + rdm1(p, q) * src
      end do
    end do

    call SysInteraction_FillInteractionPotential(hartreePotentialRaw, weightedSrc, time)
    weightedSrc = 0.0_R64

    if (mixTarget .eq. "potential") then
      if (hartreePotentialMixedInitializedQ) then
        call Mixing_Mix(hartreePotentialMixed, hartreePotentialRaw)
      else
        hartreePotentialMixed = hartreePotentialRaw
        hartreePotentialMixedInitializedQ = .true.
      end if
    else
      hartreePotentialMixed = hartreePotentialRaw
    end if

    if (.not. allocated(exchangePotentials)) then
      allocate (exchangePotentials(size(hartreePotentialMixed), nOS / 2))
      allocate (weightedPotential(size(hartreePotentialMixed)))
    end if

    ! Step 5: Diagonalize RDM-weighted Fock operator (FockAction is the matvec callback)
    call DiagonalizerList(2) % e % Diagonalize(time, .true.)

    ! Step 6: Gauge-align the Fock eigenvectors with the old orbitals (removes
    ! the arbitrary per-vector phase and the arbitrary rotation within
    ! degenerate shells), then mix and copy spin-up → spin-down (restricted)
    orbsRaw = DiagonalizerList(2) % e % evecs(:, 1:nOS / 2)
    ! ARPACK normalizes in the Euclidean metric; alignment and mixing require
    ! orthonormality with respect to the grid metric (Grid_InnerProduct)
    call Grid_Orthonormalize(orbsRaw)
    call Orbs_AlignOnReference(orbsRaw, orbs(:, 1:nOS / 2))

    if (mixTarget .eq. "orbitals") then
      orbsRawFlat(1:nG * (nOS / 2)) => orbsRaw
      call Mixing_Mix(state(nC + 1:nC + nG * (nOS / 2)), orbsRawFlat)
    else
      do j = 1, nOS / 2
        orbs(:, j) = orbsRaw(:, j)
      end do
    end if

    do j = 1, nOS / 2
      orbs(:, nOS / 2 + j) = orbs(:, j)
    end do

    ! Gram–Schmidt orthonormalize
    call Orbs_Orthonormalize(orbs)

  end subroutine

end submodule
