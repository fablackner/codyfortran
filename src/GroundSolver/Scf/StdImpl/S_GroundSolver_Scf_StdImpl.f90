! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Standard SCF implementation submodule.
!>
!> @details
!> Implements the SCF iteration on any grid back-end (the Ne3d tests run it
!> on Ylm grids; `ylmOpt` is the faster Ylm-specialized variant). The algorithm:
!>
!>   1. **Setup**: Allocate working arrays for potentials and temporaries
!>   2. **Approach** (per iteration):
!>      a. Fill external potential V̂_ext(t)
!>      b. Build raw Hartree potential Ĵ = ∫ρ(r')/|r-r'| from ONE solve of the
!>         total source ρ = 2·∑_j |φ_j|² (factor 2 for spin; the interaction
!>         potential is linear in the source; SysInteraction_FillDirectPotential
!>         with factor_ = 2); for `mixTarget="potential"`
!>         mix it into the accepted potential (Mixing_Mix), else use it directly
!>      c. With `aceQ`: build the adaptively compressed exchange operator
!>         (`BuildAce`) of the current occupied orbitals
!>      d. Diagonalize Fock operator F̂ = T̂ + V̂_ext + Ĵ − K̂ via
!>         DiagonalizerList(1) (a diagonalizer without previous eigenvectors
!>         gets the current occupied orbitals as start vectors,
!>         SetStartVectors); if it did not converge (convergedQ), a one-line
!>         note is printed and the SCF continues with its best approximations
!>         (GroundSolver_NoteUnconverged)
!>      e. Gauge-align eigenvectors with old orbitals (Orbs_AlignOnReference)
!>      f. Update orbitals: for `mixTarget="orbitals"` mix the aligned
!>         eigenvectors into the old orbitals (Mixing_Mix); for
!>         `mixTarget="potential"` replace them outright
!>      g. Copy spin-up → spin-down (restricted)
!>      h. Gram–Schmidt orthonormalize
!>
!> ## Adaptively Compressed Exchange (ACE, Lin 2016)
!>
!> With ξ_j = K̂φ_j = ∑_k W(φ_k, φ_j) φ_k for the nUp occupied orbitals
!> (nUp(nUp+1)/2 interaction-potential solves when
!> SysInteraction_conjSymmetricQ, else nUp²), M_ij = ⟨φ_i|ξ_j⟩ = L Lᴴ
!> (Cholesky) and ζ = ξ L⁻ᴴ, the operator K̂_ACE = ∑_k |ζ_k⟩⟨ζ_k| equals K̂
!> on span{φ_j} and K̂ − K̂_ACE is positive semidefinite, so the SCF fixed
!> point (lowest nUp eigenvectors) is unchanged while the Fock action costs
!> nUp inner products instead of nUp solves per matrix–vector product. The
!> construction (SysInteraction_FillExchangeAce) and the application
!> together with Ĵ in one pass (SysInteraction_AddMeanField) are the
!> grid-agnostic mean-field helpers of SysInteraction, shared with the frozen
!> mean-field reference of the Cayley exact part of Method.
!>
!> ## Preconditioner
!>
!> The preconditioner of a lobpcg Fock diagonalizer is the backend-independent
!> `GroundSolver_PreconditionAction` (M_GroundSolver, bound by
!> GroundSolver_Fabricate): (T̂ + max(−θ, σ_min))⁻¹ when the kinetic backend
!> provides `SysKinetic_SolveShifted`, else the identity.
!>
!> ## Eigensolvers and the exact exchange
!>
!> LOBPCG assumes a Fock operator self-adjoint in the grid metric. The
!> exact exchange is so only up to the discretization of the interaction:
!> on grid.fem up to the multipole truncation of the Poisson boundary values
!> (|⟨φ|K̂g⟩ − conj⟨g|K̂φ⟩| ≈ 5e-4 relative for multipoleLmax 6), enough to
!> stall LOBPCG at relative residual floors of 1e-7 (its stagnation exit
!> then stops each solve early, unconverged). K̂_ACE = ∑_k |ζ_k⟩⟨ζ_k| is
!> Hermitian by construction and costs nUp inner products instead of nUp
!> Poisson solves per Fock action (fresh H2 SCF of T_Fem_03: 5.5 s with ACE,
!> 21 s with the exact exchange, 61 s before the stagnation exit): when a
!> Fock solve with aceQ = false does not converge (the generic convergedQ of
!> the diagonalizer), Approach prints a one-time hint recommending
!> aceQ = true.
!>
!> ## JSON Configuration
!>
!> @code{.json}
!> "groundSolver": { "scf": { "stdImpl": { "mixTarget": "orbitals", "aceQ": false } } }
!> @endcode
!>
!> | Parameter      | Type   | Default    | Description                            |
!> |----------------|--------|------------|----------------------------------------|
!> | `mixTarget`    | string | "orbitals" | What the mixer damps: "orbitals" mixes |
!> |                |        |            | gauge-aligned eigenvectors into the    |
!> |                |        |            | old orbitals; "potential" mixes the    |
!> |                |        |            | Hartree potential (the Poisson solve   |
!> |                |        |            | is linear in ρ) and replaces orbitals  |
!> |                |        |            | outright; "density" purifies the mixed |
!> |                |        |            | 1-RDM (1−α)·D_old + α·D_new to its top |
!> |                |        |            | natural orbitals (Orbs_MixOccupiedSpace|
!> |                |        |            | with α from mixing.linear.alpha)       |
!> | `aceQ`         | bool   | false      | Exchange via ACE (built once per       |
!> |                |        |            | Approach) instead of nUp solves per    |
!> |                |        |            | Fock action                            |
!>
!> The mixing *algorithm* (linear, DIIS) is configured separately under the
!> top-level `mixing` block; see M_Mixing.
!>
!> ## Working Arrays
!>
!> - `hartreePotentialMixed/Raw`: direct (Hartree) potential on grid
!> - `externalPotential`: External potential V̂_ext (recomputed if time-dependent)
!> - `interactionPotential`: Temporary for exchange term computation
!> - `src`: Interaction source ρ_ij(r) = φ_i*(r)·φ_j(r) of the exact exchange
!> - `dOrbTmp`: Temporary for operator action results
!> - `aceZeta`: ACE vectors ζ_k (aceQ only)
!>
!> @note Without `aceQ` the exchange term K̂ is computed on-the-fly in
!>       `HartreeFockAction` since it depends on which orbital is being acted
!>       upon.
submodule(M_GroundSolver_Scf_StdImpl) S_GroundSolver_Scf_StdImpl

  implicit none

  !-----------------------------------------------------------------------------
  ! Module-level working arrays (allocated in Setup, used in Approach/HartreeFockAction)
  !-----------------------------------------------------------------------------

  !> Temporary for operator action results, dimension(Grid_nPoints)
  complex(R64), allocatable :: dOrbTmp(:)

  !> Accepted Hartree (direct) potential on grid, sized from
  !> SysInteraction_FillInteractionPotential on first Approach.
  !> For mixTarget="potential" this carries the mixed potential across
  !> iterations; for mixTarget="orbitals" it is rebuilt every Approach.
  complex(R64), allocatable :: hartreePotentialMixed(:)

  !> Raw Hartree potential built from the current orbitals, sized from
  !> SysInteraction_FillInteractionPotential on first Approach
  complex(R64), allocatable :: hartreePotentialRaw(:)

  !> External potential V̂_ext, dimension(Grid_nPoints) - recomputed each Approach
  complex(R64), allocatable :: externalPotential(:)

  !> Temporary for exchange potential computation, sized by
  !> SysInteraction_FillInteractionPotential
  complex(R64), allocatable :: interactionPotential(:)

  !> Interaction source density ρ_ij(r) = φ_i*(r)·φ_j(r) of the exact
  !> exchange, sized by SysInteraction_FillInteractionSrc
  complex(R64), allocatable :: src(:)

  !> Gauge-aligned new orbitals for mixing, dimension(Grid_nPoints, nOrbsInState/2)
  complex(R64), allocatable, target :: orbsRaw(:, :)

  !-----------------------------------------------------------------------------
  ! Adaptively compressed exchange (ACE)
  !-----------------------------------------------------------------------------

  !> Apply the exchange through the ACE operator ∑_k |ζ_k⟩⟨ζ_k|
  logical :: aceQ = .false.

  !> ACE vectors ζ_k, dimension(Grid_nPoints, nUp), rebuilt in every Approach
  complex(R64), allocatable :: aceZeta(:, :)

  !> Whether aceZeta holds the ACE operator (of the orbitals of the last
  !> Approach); before that, HartreeFockAction applies the exact exchange
  logical :: aceBuiltQ = .false.

  !> Whether the aceQ hint after an unconverged Fock solve with the exact
  !> exchange (aceQ = false) was printed (once per run)
  logical :: aceHintPrintedQ = .false.

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

    aceQ = Json_Get("aceQ", .false., path_="groundSolver.scf.stdImpl")

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
  !> Allocates temporaries based on the configured grid size. Potential
  !> buffers are sized lazily in Approach from the interaction layout returned
  !> by SysInteraction (as in Mcscf StdImpl).
  subroutine Setup
    use M_Utils_Say
    use M_Grid
    use M_GroundSolver
    use M_Orbs

    integer(I32) :: nUp

    call Say_Setup("groundSolver.scf.stdImpl")

    nUp = GroundSolver_NumberOfSpinUpOrbs()

    allocate (dOrbTmp(Grid_nPoints))
    allocate (orbsRaw(Grid_nPoints, nUp))
    allocate (orbsTmp(Grid_nPoints, nUp))
    allocate (lambdaDiscarded(nUp))
    if (aceQ) allocate (aceZeta(Grid_nPoints, nUp))

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
  !>            or, with `aceQ`, as K̂_ACE·φ = ∑_k ζ_k ⟨ζ_k|φ⟩ (built in Approach
  !>            from the orbitals of that iteration; exact exchange before the
  !>            first Approach), applied together with Ĵ in one pass
  !>            (SysInteraction_AddMeanField)
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

    if (aceQ .and. aceBuiltQ) then
      ! Hartree (direct) term and exchange via ACE in one pass:
      ! (Ĵ − K̂_ACE)·orb = Ĵ·orb − ∑_k ζ_k ⟨ζ_k|orb⟩ (hartreePotentialMixed
      ! and aceZeta of the last Approach; exact on its occupied space)
      call SysInteraction_AddMeanField(dOrb, dOrbTmp, orb, hartreePotentialMixed, aceZeta)
    else
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
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Build the adaptively compressed exchange (ACE) operator.
  !>
  !> @details
  !> For the occupied orbitals φ_1..φ_nUp (columns of `orbs`), closed shell
  !> (exchange within one spin, every orbital occupied once per spin):
  !> SysInteraction_FillExchangeAce with the weighted orbitals = orbs, i.e.
  !>   ξ_j = K̂φ_j = ∑_k W(φ_k, φ_j) φ_k, with the pair potentials
  !>   W(φ_k, φ_j) solved for k ≤ j only when SysInteraction_conjSymmetricQ
  !>   (W(φ_j, φ_k) = conj(W(φ_k, φ_j)) via ConjugateInteractionPotential),
  !>   else for all pairs;
  !>   M_ij = ⟨φ_i|ξ_j⟩ (Hermitian positive definite, Hermitized), M = L Lᴴ;
  !>   ζ = ξ L⁻ᴴ (stored in `aceZeta`),
  !> so that ∑_k ζ_k ⟨ζ_k|φ⟩ = ξ M⁻¹ ⟨ξ|φ⟩ reproduces K̂φ for φ in span{φ_j}.
  !> If M is numerically not positive definite (e.g. a slightly non-symmetric
  !> discretized interaction kernel), the pseudo-inverse of M on its positive
  !> eigen-subspace is used and a warning is printed.
  !>
  !> @param[in] orbs  Occupied orbitals (Grid_nPoints x nUp)
  !> @param[in] time  Time for the interaction potential
  subroutine BuildAce(orbs, time)
    use M_SysInteraction

    complex(R64), intent(in), contiguous :: orbs(:, :)
    real(R64), intent(in) :: time

    logical :: positiveDefiniteQ

    call SysInteraction_FillExchangeAce(aceZeta, positiveDefiniteQ, orbs, orbs, time)
    if (.not. positiveDefiniteQ) &
      write (*, *) "Warning: groundSolver.scf.stdImpl.aceQ: exchange matrix <phi_i|K phi_j> not positive definite, ", &
      "using its positive eigen-subspace"
    aceBuiltQ = .true.

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Perform one SCF iteration: build Fock, diagonalize, mix, orthonormalize.
  !>
  !> @details
  !> Algorithm:
  !>   1. Fill external potential V̂_ext for current time
  !>   2. Build Hartree potential Ĵ from one solve of the total source
  !>      2·∑_j |φ_j|² (factor 2 for spin degeneracy); with aceQ build the ACE
  !>      exchange operator of the current orbitals
  !>   3. Diagonalize Fock operator via DiagonalizerList(1) (uses HartreeFockAction
  !>      callback); a diagonalizer without previous eigenvectors gets the
  !>      current occupied orbitals as start vectors (SetStartVectors)
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

    ! Step 2: Build raw Hartree potential from occupied orbitals with ONE
    ! solve of the total source (the potential is linear in the source).
    ! Factor 2 accounts for spin degeneracy (restricted: spin-up = spin-down)
    call SysInteraction_FillDirectPotential(hartreePotentialRaw, orbs(:, 1:nUp), orbs(:, 1:nUp), time, factor_=2.0_R64)

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

    ! Exchange via ACE: compress K̂ of the current orbitals once per iteration
    if (aceQ) call BuildAce(orbs(:, 1:nUp), time)

    ! Step 3: Diagonalize Fock operator (HartreeFockAction is the matvec callback).
    ! A diagonalizer without previous eigenvectors gets the current occupied
    ! orbitals as start vectors of its first solve (backends that cannot
    ! start from given vectors ignore them)
    if (.not. allocated(DiagonalizerList(1) % e % evecs) .and. DiagonalizerList(1) % e % dim .eq. nG) then
      call DiagonalizerList(1) % e % SetStartVectors(orbs(:, 1:nUp))
    end if
    call DiagonalizerList(1) % e % Diagonalize(time, .true.)
    call GroundSolver_NoteUnconverged(1)
    ! An unconverged solve with the exact exchange: it is self-adjoint in the
    ! metric only up to its discretization (grid.fem: the multipole truncation
    ! of the Poisson boundary values), which can stall a metric-based solver
    ! (LOBPCG) at a residual floor; ACE is self-adjoint by construction
    if (.not. aceQ .and. .not. DiagonalizerList(1) % e % convergedQ .and. .not. aceHintPrintedQ) then
      write (*, '(1X, 3A)') "hint: groundSolver.scf.stdImpl: ", DiagonalizerList(1) % e % path, &
        " did not converge with aceQ = false: the exact exchange is Hermitian only up to its discretization"
      write (*, '(1X, A)') "  (grid.fem: multipole truncation of the Poisson boundary values), which can stall an "// &
        "eigensolver at a residual floor; aceQ = true is recommended (Hermitian by construction, and much "// &
        "cheaper). Printed once."
      aceHintPrintedQ = .true.
    end if

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
