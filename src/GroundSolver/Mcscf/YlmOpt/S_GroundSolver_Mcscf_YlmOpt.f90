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
!>      (GroundSolver_Mcscf_HamiltonianAction is the matvec callback)
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

  !> Accumulated RDM-weighted interaction source for the direct term, dimension(potSize)
  complex(R64), allocatable :: weightedSrc(:)

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

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine GroundSolver_Mcscf_YlmOpt_Fabricate
    use M_Utils_Say
    use M_Utils_Json
    use M_GroundSolver

    call Say_Fabricate("groundSolver.mcscf.ylmOpt")

    !------------------------------------
    ! read configuration parameters
    !------------------------------------

    alphaCi = Json_Get("alphaCi", 1.0_R64, path_="groundSolver.mcscf.ylmOpt")
    if (.not. (alphaCi > 0.0_R64 .and. alphaCi <= 1.0_R64)) error stop "groundSolver.mcscf.ylmOpt.alphaCi must be in (0, 1]"

    mixTarget = Json_Get("mixTarget", "orbitals", path_="groundSolver.mcscf.ylmOpt")
    if (mixTarget .ne. "orbitals" .and. mixTarget .ne. "potential") then
      error stop "groundSolver.mcscf.ylmOpt.mixTarget must be one of: orbitals, potential"
    end if

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    GroundSolver_Setup => Setup
    GroundSolver_Approach => Approach
    GroundSolver_Mcscf_YlmOpt_FockAction => FockAction
    GroundSolver_implementsDeclaredMQ = .true.

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Allocate working arrays for the Ylm-optimized MCSCF implementation.
  !>
  !> @details
  !> Allocates radial potentials and full-grid temporaries. The potential size
  !> holds all (l,m) channels of the interaction expansion × nRadial.
  subroutine Setup
    use M_Utils_Json
    use M_Utils_Say
    use M_Grid
    use M_Grid_Ylm
    use M_SysInteraction_Ylm
    use M_OrbsInit_Ylm_HydrogenLike
    use M_GroundSolver
    integer(I32) :: potSize, nUp

    call Say_Setup("groundSolver.mcscf.ylmOpt")

    ! The eigenvectors are mapped to the orbitals by the (n, l, m) of their
    ! hydrogen-like start
    if (.not. Json_GetExistence("orbsInit.ylm.hydrogenLike") .or. Json_GetExistence("orbsInit.load")) &
      error stop "groundSolver.mcscf.ylmOpt takes the (n, l, m) of every orbital from orbsInit.ylm.hydrogenLike &
                 &(without orbsInit.load)"
    if (.not. allocated(OrbsInit_Ylm_HydrogenLike_n)) &
      error stop "groundSolver.mcscf.ylmOpt: call OrbsInit_Fabricate before GroundSolver_Setup"
    if (size(OrbsInit_Ylm_HydrogenLike_n) < GroundSolver_NumberOfSpinUpOrbs()) &
      error stop "groundSolver.mcscf.ylmOpt: orbsInit.ylm.hydrogenLike needs (n, l, m) for every spin-up orbital"

    potSize = SysInteraction_Ylm_nPoints

    nUp = GroundSolver_NumberOfSpinUpOrbs()

    if (allocated(weightedSrc)) deallocate (weightedSrc)
    allocate (weightedSrc(potSize))
    if (allocated(dOrbLmTmp)) deallocate (dOrbLmTmp)
    allocate (dOrbLmTmp(Grid_Ylm_nRadial))
    if (allocated(potLm)) deallocate (potLm)
    allocate (potLm(Grid_Ylm_nRadial))
    if (allocated(hartreePotentialRaw)) deallocate (hartreePotentialRaw)
    allocate (hartreePotentialRaw(Grid_Ylm_nRadial))
    if (allocated(srcLm)) deallocate (srcLm)
    allocate (srcLm(Grid_Ylm_nRadial))
    if (allocated(orb)) deallocate (orb)
    allocate (orb(Grid_nPoints))
    if (allocated(dOrbTmp)) deallocate (dOrbTmp)
    allocate (dOrbTmp(Grid_nPoints))
    if (allocated(orbsRaw)) deallocate (orbsRaw)
    allocate (orbsRaw(Grid_nPoints, nUp))
    if (allocated(rdm1Up)) deallocate (rdm1Up)
    allocate (rdm1Up(nUp, nUp))

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
    use M_GroundSolver_Mcscf
    use M_SysKinetic_Ylm
    use M_Grid_Ylm

    complex(R64), intent(out), contiguous, target :: dOrbLm(:)
    complex(R64), intent(in), contiguous, target :: orbLm(:)
    integer(I32), intent(in) :: l
    real(R64), intent(in) :: time

    real(R64) :: gVal
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

    dOrbTmp = 0.0_R64
    call GroundSolver_Mcscf_SubtractExchange(dOrbTmp, orb, rdm1Up, time)
    call Grid_Ylm_GetLmComponent(dOrbLmTmp, l, 0, dOrbTmp)
    dOrbLm = dOrbLm + dOrbLmTmp

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
    use M_Grid
    use M_GroundSolver_Mcscf
    use M_Grid_Ylm
    use M_DiagonalizerList
    use M_GroundSolver
    use M_Coeffs
    use M_Orbs
    use M_SysInteraction_Ylm
    use M_SysPotential_Ylm
    use M_OrbsInit_Ylm_HydrogenLike
    use M_Method_Mb_OrbBased

    complex(R64), intent(inout), contiguous, target :: state(:)
    real(R64), intent(in) :: time

    complex(R64), contiguous, pointer :: coeffs(:)
    complex(R64), contiguous, pointer :: orbs(:, :)
    complex(R64), allocatable :: rdm1(:, :)
    integer(I32) :: nC, nG, nOS, nUp, j
    integer(I32) :: n, l, m, nr

    nC = Coeffs_nCoeffs
    nG = Grid_nPoints
    nOS = Orbs_nOrbsInState
    nUp = GroundSolver_NumberOfSpinUpOrbs()

    coeffs(1:nC) => state(1:nC)
    orbs(1:nG, 1:nOS) => state(nC + 1:)

    ! Step 1: Hamiltonian matrices in the current basis
    call Method_Mb_OrbBased_FillH1(Method_Mb_OrbBased_h1, orbs, time)
    call Method_Mb_OrbBased_FillH2(Method_Mb_OrbBased_h2, orbs, time)

    ! Steps 2 and 3: diagonalize the CI Hamiltonian, mix the phase-aligned
    ! ground eigenvector into the coefficients
    call GroundSolver_Mcscf_CiStep(coeffs, alphaCi, time)

    ! Step 4: 1-RDM from updated coefficients; monopole radial Hartree
    ! potential from the correlated density
    ! ρ(r) = Σ_{pq} ρ¹_{pq}·φ_p*(r)·φ_q(r) (all body types)
    call Coeffs_ApplyH1FillRdm1(coeffs, rdm1_=rdm1)
    rdm1Up = rdm1(1:nUp, 1:nUp)

    call GroundSolver_Mcscf_FillCorrelatedSource(weightedSrc, rdm1, orbs)

    call Grid_Ylm_GetLmComponent(srcLm, 0, 0, weightedSrc, mmax_=SysInteraction_Ylm_mmax)
    call SysInteraction_Ylm_FillInteractionPotentialRadial(hartreePotentialRaw, srcLm, 0, 0, time)

    if (mixTarget .eq. "potential") then
      call GroundSolver_MixPotential(hartreePotentialMixed, hartreePotentialRaw)
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
      call GroundSolver_NoteUnconverged(l + 2)
    end do

    ! Step 6: Map eigenvectors to orbitals using (n,l,m)
    do j = 1, nUp
      n = OrbsInit_Ylm_HydrogenLike_n(j)
      l = OrbsInit_Ylm_HydrogenLike_l(j)
      m = OrbsInit_Ylm_HydrogenLike_m(j)
      nr = n - l - 1  ! Radial quantum number (node count)

      orbsRaw(:, j) = 0.0_R64
      call Grid_Ylm_SetLmComponent(orbsRaw(:, j), l, m, DiagonalizerList(l + 2) % e % evecs(:, nr + 1))
    end do

    ! The radial eigenvectors are normalized in the Euclidean metric:
    ! orthonormalize in the grid metric and gauge-align on the old orbitals
    ! (per block of equal declared m: never mixes orbitals of different m),
    ! then mix or replace, copy spin-up → spin-down and orthonormalize
    call GroundSolver_AlignOrbitals(orbsRaw, orbs(:, 1:nUp), blockOf_=Orbs_blockOfOrb(1:nUp))
    call GroundSolver_StoreOrbitals(state(nC + 1:nC + nG * nOS), orbsRaw, mixQ=mixTarget .eq. "orbitals")

  end subroutine

end submodule
