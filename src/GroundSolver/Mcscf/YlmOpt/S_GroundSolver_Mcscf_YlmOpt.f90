! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Ylm-optimized MCSCF implementation submodule.
!>
!> @details
!> Implements the MCSCF iteration exploiting spherical symmetry. The CI step
!> is identical to the standard implementation; the orbital step diagonalizes
!> smaller nRadial×nRadial radial RDM-weighted Fock operators per l-channel
!> instead of one Grid_nPoints×Grid_nPoints matrix (cf. SCF YlmOpt). The
!> packed state layout matches the MCTDHX method:
!> state(1:nC) = CI coefficients, state(nC+1:) = orbitals (column-major).
!>
!> ## Algorithm
!>
!> **Setup**: Allocate radial potentials and full-grid temporaries.
!>
!> **Approach** (per iteration):
!>   1. Build h1/h2 in the current orbital basis
!>   2. Diagonalize the CI Hamiltonian via DiagonalizerList(1)
!>      (HamiltonianAction is the matvec callback)
!>   3. Mix the phase-aligned CI ground eigenvector with the old
!>      coefficients: c_new = (1-αCI)·c_old + αCI·evec, then normalize
!>   4. Fill the 1-RDM from the updated coefficients, accumulate the
!>      RDM-weighted interaction source of the correlated density
!>      ρ(r) = Σ_{pq} ρ¹_{pq} φ_p*(r)·φ_q(r), extract its monopole (l=0)
!>      component, and solve for the radial Hartree potential
!>   5. Combine with the radial external potential and diagonalize each
!>      l-channel's radial Fock operator via DiagonalizerList(l+2)
!>      (l-wrappers around FockAction are the matvec callbacks)
!>   6. Map eigenvectors back using hydrogen-like quantum numbers (n,l,m),
!>      gauge-align with the old orbitals, mix, copy spin-up → spin-down
!>      (restricted), Gram–Schmidt orthonormalize
!>
!> **FockAction** (per l-channel):
!>   - Radial kinetic + centrifugal: SysKinetic_Ylm_MultiplyWithRadialKineticOp
!>   - Mean-field potential: Gaunt-weighted monopole (external + Hartree[ρ¹])
!>   - Exchange: RDM-weighted with full angular coupling via
!>     Grid_Ylm_{Set,Get}LmComponent
!>
!> ## Working Arrays
!>
!> - `h1`, `h2`: One-/two-body Hamiltonian matrices in the current orbital
!>   basis; rebuilt each `Approach` and shared with `HamiltonianAction`
!> - `rdm1Up`: Spin-up block of the correlated 1-RDM, used by `FockAction`
!> - `potLm`: Combined radial potential (external + Hartree[ρ¹], l=0)
!> - `weightedSrc`: RDM-weighted interaction source of the correlated density
!> - `exchangePotentials`: Per-orbital exchange potentials W(φ_p, orb)
!>
!> @note Requires OrbsInit_Ylm_HydrogenLike to track (n,l,m) quantum numbers.
!> @note The monopole truncation of the direct potential is exact for
!>       spherically symmetric (S) states.
!> @note Restricted calculations only: the spin-up orbitals are relaxed and
!>       copied to the spin-down block.
submodule(M_GroundSolver_Mcscf_YlmOpt) S_GroundSolver_Mcscf_YlmOpt

  implicit none

  !-----------------------------------------------------------------------------
  ! Module-level working arrays (allocated in Setup, used in Approach/callbacks)
  !-----------------------------------------------------------------------------

  !> One-body Hamiltonian matrix in the current orbital basis, filled in Approach
  complex(R64), allocatable :: h1(:, :)

  !> Two-body Hamiltonian matrix in the current orbital basis, filled in Approach
  complex(R64), allocatable :: h2(:, :, :, :)

  !> Spin-up block of the correlated 1-RDM, filled in Approach
  complex(R64), allocatable :: rdm1Up(:, :)

  !> Temporary for radial operator action, dimension(Grid_Ylm_nRadial)
  complex(R64), allocatable :: dOrbLmTmp(:)

  !> Combined radial potential (external + Hartree[ρ¹]), dimension(Grid_Ylm_nRadial)
  complex(R64), allocatable :: potLm(:)

  !> Accepted radial Hartree (direct) potential, dimension(Grid_Ylm_nRadial)
  complex(R64), allocatable :: hartreePotentialMixed(:)

  !> Raw radial Hartree potential from the correlated density, dimension(Grid_Ylm_nRadial)
  complex(R64), allocatable :: hartreePotentialRaw(:)

  !> Radial source (l=0 component of the correlated density), dimension(Grid_Ylm_nRadial)
  complex(R64), allocatable :: srcLm(:)

  !> Temporary for pairwise interaction sources, dimension(potSize)
  complex(R64), allocatable :: src(:)

  !> Accumulated RDM-weighted interaction source for the direct term, dimension(potSize)
  complex(R64), allocatable :: weightedSrc(:)

  !> Temporary for full-grid interaction potential, dimension(potSize)
  complex(R64), allocatable :: interactionPotential(:)

  !> Per-orbital exchange potentials W(φ_p, orb), dimension(potSize, nOrbsInState/2)
  complex(R64), allocatable :: exchangePotentials(:, :)

  !> RDM-weighted sum of exchange potentials, dimension(potSize)
  complex(R64), allocatable :: weightedPotential(:)

  !> Full-grid orbital temporary for exchange, dimension(Grid_nPoints)
  complex(R64), allocatable :: orb(:)

  !> Full-grid action temporary for exchange, dimension(Grid_nPoints)
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
  module subroutine GroundSolver_Mcscf_YlmOpt_Fabricate
    use M_Utils_Say
    use M_Utils_Json
    use M_GroundSolver
    use M_GroundSolver_Mcscf

    call Say_Fabricate("groundSolver.mcscf.ylmOpt")

    !------------------------------------
    ! read configuration parameters
    !------------------------------------

    alphaCi = Json_Get("alphaCi", 1.0_R64, path_="groundSolver.mcscf.ylmOpt")

    mixTarget = Json_Get("mixTarget", "orbitals", path_="groundSolver.mcscf.ylmOpt")
    if (mixTarget .ne. "orbitals" .and. mixTarget .ne. "potential") then
      error stop "groundSolver.mcscf.ylmOpt.mixTarget must be one of: orbitals, potential"
    end if

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    GroundSolver_Setup => Setup
    GroundSolver_Approach => Approach
    GroundSolver_Mcscf_HamiltonianAction => HamiltonianAction
    GroundSolver_Mcscf_YlmOpt_FockAction => FockAction

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Allocate working arrays for the Ylm-optimized MCSCF implementation.
  !>
  !> @details
  !> Allocates radial potentials and full-grid temporaries. The potential size
  !> is (2·lmax_pot+1)² × nRadial to accommodate all angular momentum channels.
  subroutine Setup
    use M_Utils_Say
    use M_Grid
    use M_Grid_Ylm
    use M_SysInteraction_Ylm
    use M_Orbs

    integer(I32) :: lmaxPot, potSize

    call Say_Setup("groundSolver.mcscf.ylmOpt")

    lmaxPot = SysInteraction_Ylm_lmax
    potSize = (2 * lmaxPot + 1)**2 * Grid_Ylm_nRadial

    allocate (src(potSize))
    allocate (weightedSrc(potSize))
    allocate (interactionPotential(potSize))
    allocate (exchangePotentials(potSize, Orbs_nOrbsInState / 2))
    allocate (weightedPotential(potSize))
    allocate (dOrbLmTmp(Grid_Ylm_nRadial))
    allocate (potLm(Grid_Ylm_nRadial))
    allocate (hartreePotentialMixed(Grid_Ylm_nRadial))
    allocate (hartreePotentialRaw(Grid_Ylm_nRadial))
    allocate (srcLm(Grid_Ylm_nRadial))
    allocate (orb(Grid_nPoints))
    allocate (dOrbTmp(Grid_nPoints))
    allocate (orbsRaw(Grid_nPoints, Orbs_nOrbsInState / 2))
    allocate (rdm1Up(Orbs_nOrbsInState / 2, Orbs_nOrbsInState / 2))

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
  !> @brief Apply the radial RDM-weighted Fock operator for channel l:
  !>        dOrbLm = F̂_l[ρ¹]·orbLm.
  !>
  !> @details
  !> Computes:
  !>   1. Radial kinetic + centrifugal: T̂_l·orbLm via SysKinetic_Ylm
  !>   2. Mean-field potential: G(l,0;l,0;0,0)·potLm·orbLm where G is the
  !>      Gaunt coefficient and potLm the combined external + Hartree[ρ¹]
  !>      radial potential (precomputed in Approach)
  !>   3. Exchange: RDM-weighted with full angular coupling,
  !>      K̂·φ = Σ_{pq} ρ¹_{pq}·W(φ_p, φ)·φ_q; computed by expanding the trial
  !>      orbital to the full grid and extracting the l-component afterwards
  !>
  !> For ρ¹ = identity on the occupied space this reduces to the SCF YlmOpt
  !> Hartree–Fock action.
  !>
  !> This is the matvec callback for DiagonalizerList(l+2) over the radial space.
  subroutine FockAction(dOrbLm, orbLm, l, time)
    use M_Utils_SphericalHarmonics
    use M_SysKinetic_Ylm
    use M_SysInteraction
    use M_Orbs
    use M_Grid_Ylm

    complex(R64), intent(out), contiguous, target :: dOrbLm(:)
    complex(R64), intent(in), contiguous, target :: orbLm(:)
    integer(I32), intent(in) :: l
    real(R64), intent(in) :: time

    real(R64) :: gVal
    integer(I32) :: p, q, nUp

    nUp = Orbs_nOrbsInState / 2

    dOrbLm = 0.0_R64

    ! Gaunt coefficient for monopole coupling: <Y_l0|Y_00|Y_l0>
    gVal = SphericalHarmonics_GauntCoefficient( &
           l1=0, m1=0, &
           l2=l, m2=0, &
           l3=l, m3=0)

    ! Radial kinetic + centrifugal: T̂_l·orbLm
    call SysKinetic_Ylm_MultiplyWithRadialKineticOp(dOrbLmTmp, orbLm, l, 0, time)
    dOrbLm = dOrbLm + dOrbLmTmp

    ! Mean-field potential (Gaunt-weighted monopole): potLm filled in Approach
    dOrbLm = dOrbLm + gVal * potLm * orbLm

    ! Exchange: expand to full grid, compute RDM-weighted exchange, extract
    ! l-component. Spin-up block only (restricted; the trial vector carries
    ! no spin index): −K̂[ρ¹]·orb = −Σ_{pq} ρ¹_{pq}·W(φ_p, orb)·φ_q
    orb = 0.0_R64
    call Grid_Ylm_SetLmComponent(orb, l, 0, orbLm)

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
      call Grid_Ylm_GetLmComponent(dOrbLmTmp, l, 0, dOrbTmp)
      dOrbLm = dOrbLm - dOrbLmTmp
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Perform one Ylm-optimized MCSCF iteration: CI diagonalization +
  !>        per-l radial Fock diagonalizations.
  !>
  !> @details
  !> Algorithm:
  !>   1. Build h1/h2 in the current orbital basis
  !>   2. Diagonalize the CI Hamiltonian via DiagonalizerList(1)
  !>   3. Mix the phase-aligned CI ground eigenvector with the old coefficients
  !>      and normalize
  !>   4. Fill the 1-RDM from the updated coefficients, build the monopole
  !>      radial Hartree potential from the correlated density
  !>      ρ(r) = Σ_{pq} ρ¹_{pq}·φ_p*(r)·φ_q(r) (all body types)
  !>   5. Diagonalize each l-channel's radial Fock operator via
  !>      DiagonalizerList(l+2)
  !>   6. Map eigenvectors to orbitals using (n,l,m) quantum numbers,
  !>      gauge-align with the old orbitals, mix, copy spin-up → spin-down,
  !>      Gram–Schmidt orthonormalize
  !>
  !> The orbital-to-eigenvector mapping uses the radial quantum number
  !> nr = n - l - 1 to select the correct eigenvector from the l-channel
  !> diagonalization.
  !>
  !> @param[inout] state  Packed MCTDHX-style state: coefficients ⊕ orbitals
  !> @param[in]    time   Time for potential evaluation (typically 0)
  subroutine Approach(state, time)
    use M_Mixing
    use M_Grid
    use M_Grid_Ylm
    use M_DiagonalizerList
    use M_Coeffs
    use M_Orbs
    use M_SysInteraction
    use M_SysInteraction_Ylm
    use M_SysPotential_Ylm
    use M_OrbsInit_Ylm_HydrogenLike
    use M_Method_Mb_OrbBased

    complex(R64), intent(inout), contiguous, target :: state(:)
    real(R64), intent(in) :: time

    complex(R64), contiguous, pointer :: coeffs(:)
    complex(R64), contiguous, pointer :: orbs(:, :)
    complex(R64), contiguous, pointer :: orbsRawFlat(:)
    complex(R64), allocatable :: rdm1(:, :)
    complex(R64) :: overlap, phase
    integer(I32) :: nC, nG, nOS, p, q, j
    integer(I32) :: n, l, m, nr

    nC = Coeffs_nCoeffs
    nG = Grid_nPoints
    nOS = Orbs_nOrbsInState

    coeffs(1:nC) => state(1:nC)
    orbs(1:nG, 1:nOS) => state(nC + 1:)

    ! Step 1: Hamiltonian matrices in the current basis
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

    ! Step 4: 1-RDM from updated coefficients; monopole radial Hartree
    ! potential from the correlated density
    ! ρ(r) = Σ_{pq} ρ¹_{pq}·φ_p*(r)·φ_q(r) (all body types)
    call Coeffs_ApplyH1FillRdm1(coeffs, rdm1_=rdm1)
    rdm1Up = rdm1(1:nOS / 2, 1:nOS / 2)

    weightedSrc = 0.0_R64
    do p = 1, nOS
      do q = 1, nOS
        if (abs(rdm1(p, q)) .eq. 0.0_R64) cycle
        call SysInteraction_FillInteractionSrc(src, orbs(:, p), orbs(:, q))
        weightedSrc = weightedSrc + rdm1(p, q) * src
      end do
    end do

    call Grid_Ylm_GetLmComponent(srcLm, 0, 0, weightedSrc)
    call SysInteraction_Ylm_FillInteractionPotentialRadial(hartreePotentialRaw, srcLm, 0, 0, time)

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

    ! Combined radial potential: external (l=0, m=0 component) + Hartree[ρ¹]
    call SysPotential_Ylm_FillExternalPotentialRadial(potLm, 0, 0, time)
    potLm = potLm + hartreePotentialMixed

    ! Step 5: Diagonalize each l-channel's radial Fock operator (l-wrappers
    ! around FockAction are the matvec callbacks)
    do l = 0, Grid_Ylm_lmax
      call DiagonalizerList(l + 2) % e % Diagonalize(time, .true.)
    end do

    ! Step 6: Map eigenvectors to orbitals using (n,l,m)
    do j = 1, nOS / 2
      n = OrbsInit_Ylm_HydrogenLike_n(j)
      l = OrbsInit_Ylm_HydrogenLike_l(j)
      m = OrbsInit_Ylm_HydrogenLike_m(j)
      nr = n - l - 1  ! Radial quantum number (node count)

      orbsRaw(:, j) = 0.0_R64
      call Grid_Ylm_SetLmComponent(orbsRaw(:, j), l, m, DiagonalizerList(l + 2) % e % evecs(:, nr + 1))
    end do

    ! The radial eigenvectors are normalized in the Euclidean metric;
    ! alignment and mixing require orthonormality with respect to the grid
    ! metric (Grid_InnerProduct), then gauge-align the new orbitals with the
    ! old ones (removes the arbitrary phase of the radial eigenvectors),
    ! then mix and copy spin-up → spin-down (restricted)
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
