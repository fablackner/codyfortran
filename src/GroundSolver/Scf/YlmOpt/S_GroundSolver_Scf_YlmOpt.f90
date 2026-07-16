! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Ylm-optimized SCF implementation submodule.
!>
!> @details
!> Implements the SCF iteration exploiting spherical symmetry. For atoms with
!> a central potential, orbitals factorize as φ_{nlm}(r) = R_nl(r)·Y_lm(θ,φ).
!> This allows diagonalizing smaller nRadial×nRadial matrices per l-channel
!> instead of a single large Grid_nPoints×Grid_nPoints matrix.
!>
!> ## Algorithm
!>
!> **Setup**: Allocate working arrays for radial potentials and sources.
!>
!> **Approach** (per iteration):
!>   1. Extract l=0 radial component of external potential
!>   2. Build monopole Hartree potential from |φ_j|² (l=0 component)
!>   3. Diagonalize per-l radial Fock operators via DiagonalizerList(l+1)
!>   4. Map eigenvectors back using hydrogen-like quantum numbers (n,l,m)
!>   5. Gauge-align with old orbitals, mix, copy spin-up → spin-down,
!>      orthonormalize
!>
!> **HartreeFockAction** (per l-channel):
!>   - Radial kinetic + centrifugal: SysKinetic_Ylm_MultiplyWithRadialKineticOp
!>   - Mean-field potential: Gaunt-weighted monopole
!>   - Exchange: full angular coupling via Grid_Ylm_{Set,Get}LmComponent
!>
!> ## Working Arrays
!>
!> - `potLm`: Combined external + Hartree radial potential (l=0)
!> - `src`: Full angular interaction source density
!> - `orb`, `dOrbTmp`: Full-grid temporaries for exchange computation
!>
!> @note Requires OrbsInit_Ylm_HydrogenLike to track (n,l,m) quantum numbers.
submodule(M_GroundSolver_Scf_YlmOpt) S_GroundSolver_Scf_YlmOpt

  implicit none

  !-----------------------------------------------------------------------------
  ! Module-level working arrays (allocated in Setup, used in Approach/HartreeFockAction)
  !-----------------------------------------------------------------------------

  !> Temporary for radial operator action, dimension(Grid_Ylm_nRadial)
  complex(R64), allocatable :: dOrbLmTmp(:)

  !> Combined radial potential (external + Hartree), dimension(Grid_Ylm_nRadial)
  complex(R64), allocatable :: potLm(:)

  !> Accepted Hartree (direct) potential on grid, dimension(Grid_Ylm_nRadial)
  complex(R64), allocatable :: hartreePotentialMixed(:)

  !> Raw Hartree potential built from the current orbitals, dimension(Grid_Ylm_nRadial)
  complex(R64), allocatable :: hartreePotentialRaw(:)

  !> Full angular interaction source density, dimension(potSize)
  complex(R64), allocatable :: src(:)

  !> Temporary source for summation, dimension(potSize)
  complex(R64), allocatable :: srcTmp(:)

  !> Radial source (l=0 component), dimension(Grid_Ylm_nRadial)
  complex(R64), allocatable :: srcLm(:)

  !> Temporary for full-grid interaction potential, dimension(potSize)
  complex(R64), allocatable :: interactionPotential(:)

  !> Full-grid orbital temporary for exchange, dimension(Grid_nPoints)
  complex(R64), allocatable :: orb(:)

  !> Full-grid action temporary for exchange, dimension(Grid_nPoints)
  complex(R64), allocatable :: dOrbTmp(:)

  !> Gauge-aligned new orbitals for mixing, dimension(Grid_nPoints, nOrbsInState/2)
  complex(R64), allocatable, target :: orbsRaw(:, :)

  !-----------------------------------------------------------------------------
  ! Mixing configuration
  !-----------------------------------------------------------------------------

  !> What quantity the mixer acts on: "orbitals" or "potential"
  character(len=:), allocatable :: mixTarget

  !> Whether hartreePotentialMixed already holds an accepted (mixed) potential
  logical :: hartreePotentialMixedInitializedQ = .false.

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine GroundSolver_Scf_YlmOpt_Fabricate
    use M_Utils_Say
    use M_Utils_Json
    use M_GroundSolver
    use M_GroundSolver_Scf

    call Say_Fabricate("groundSolver.scf.ylmOpt")

    !------------------------------------
    ! read configuration parameters
    !------------------------------------

    mixTarget = Json_Get("mixTarget", "orbitals", path_="groundSolver.scf.ylmOpt")
    if (mixTarget .ne. "orbitals" .and. mixTarget .ne. "potential") then
      error stop "groundSolver.scf.ylmOpt.mixTarget must be one of: orbitals, potential"
    end if

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    GroundSolver_Setup => Setup
    GroundSolver_Approach => Approach
    GroundSolver_Scf_YlmOpt_HartreeFockAction => HartreeFockAction

  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Allocate working arrays for the Ylm-optimized SCF implementation.
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

    call Say_Setup("groundSolver.scf.ylmOpt")

    lmaxPot = SysInteraction_Ylm_lmax
    potSize = (2 * lmaxPot + 1)**2 * Grid_Ylm_nRadial

    allocate (src(potSize))
    allocate (srcTmp(potSize))
    allocate (interactionPotential(potSize))
    allocate (dOrbLmTmp(Grid_Ylm_nRadial))
    allocate (potLm(Grid_Ylm_nRadial))
    allocate (hartreePotentialMixed(Grid_Ylm_nRadial))
    allocate (hartreePotentialRaw(Grid_Ylm_nRadial))
    allocate (srcLm(Grid_Ylm_nRadial))
    allocate (orb(Grid_nPoints))
    allocate (dOrbTmp(Grid_nPoints))
    allocate (orbsRaw(Grid_nPoints, Orbs_nOrbsInState / 2))

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Apply the radial Fock operator for channel l: dOrbLm = F̂_l·orbLm.
  !>
  !> @details
  !> Computes:
  !>   1. Radial kinetic + centrifugal: T̂_l·orbLm via SysKinetic_Ylm
  !>   2. Mean-field potential: G(l,0;l,0;0,0)·potLm·orbLm where G is the Gaunt coefficient
  !>   3. Exchange: full angular coupling to all occupied orbitals
  !>
  !> The Gaunt coefficient G(l,0;l,0;0,0) = <Y_l0|Y_00|Y_l0> couples the monopole
  !> potential to the l-channel. Exchange is computed by expanding the trial
  !> orbital to the full grid, computing the exchange there, and extracting the
  !> l-component.
  subroutine HartreeFockAction(dOrbLm, orbLm, l, time)
    use M_Utils_SphericalHarmonics
    use M_SysKinetic_Ylm
    use M_SysPotential_Ylm
    use M_SysInteraction
    use M_Method_Mb_OrbBased
    use M_Orbs
    use M_Grid
    use M_Grid_Ylm

    complex(R64), intent(out), contiguous, target :: dOrbLm(:)
    complex(R64), intent(in), contiguous, target :: orbLm(:)
    integer(I32), intent(in) :: l
    real(R64), intent(in) :: time

    real(R64) :: gVal
    integer(I32) :: j

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

    ! Exchange: expand to full grid, compute exchange, extract l-component
    orb = 0.0_R64
    call Grid_Ylm_SetLmComponent(orb, l, 0, orbLm)
    do j = 1, Orbs_nOrbsInState / 2
      call SysInteraction_FillInteractionSrc(src, Orbs_orbs(:, j), orb(:))
      call SysInteraction_FillInteractionPotential(interactionPotential, src, time)
      call SysInteraction_MultiplyWithInteractionPotential(dOrbTmp, interactionPotential, Orbs_orbs(:, j))
      call Grid_Ylm_GetLmComponent(dOrbLmTmp, l, 0, dOrbTmp)
      dOrbLm = dOrbLm - dOrbLmTmp
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Perform one Ylm-optimized SCF iteration.
  !>
  !> @details
  !> Algorithm:
  !>   1. Extract radial external potential V̂_ext(r) for l=0,m=0
  !>   2. Build monopole Hartree potential from occupied orbitals
  !>   3. Diagonalize each l-channel's radial Fock operator
  !>   4. Map eigenvectors to orbitals using (n,l,m) quantum numbers
  !>   5. Gauge-align with old orbitals, mix, copy spin-up → spin-down
  !>   6. Gram–Schmidt orthonormalize
  !>
  !> The orbital-to-eigenvector mapping uses the radial quantum number nr = n - l - 1
  !> to select the correct eigenvector from the l-channel diagonalization.
  !>
  !> @param[inout] state  Packed orbital state
  !> @param[in]    time   Time for potential evaluation
  subroutine Approach(state, time)
    use M_Mixing
    use M_Grid
    use M_Grid_Ylm
    use M_DiagonalizerList
    use M_Orbs
    use M_SysInteraction
    use M_SysInteraction_Ylm
    use M_SysPotential_Ylm
    use M_OrbsInit_Ylm_HydrogenLike

    complex(R64), intent(inout), contiguous, target :: state(:)
    real(R64), intent(in) :: time

    complex(R64), contiguous, pointer :: orbs(:, :)
    complex(R64), contiguous, pointer :: orbsRawFlat(:)
    integer(I32) :: nG, nOS, j
    integer(I32) :: n, l, m, nr

    nG = Grid_nPoints
    nOS = Orbs_nOrbsInState
    orbs(1:nG, 1:nOS) => state(1:)

    ! Step 1: Extract radial external potential (l=0, m=0 component)
    call SysPotential_Ylm_FillExternalPotentialRadial(potLm, 0, 0, time)

    ! Step 2: Build monopole (l=0) Hartree potential from occupied orbitals
    src = 0.0_R64
    do j = 1, nOS / 2
      call SysInteraction_FillInteractionSrc(srcTmp, orbs(:, j), orbs(:, j))
      src = src + 2.0_R64 * srcTmp  ! Factor 2 for spin degeneracy
    end do
    call Grid_Ylm_GetLmComponent(srcLm, 0, 0, src)
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

    potLm = potLm + hartreePotentialMixed

    ! Step 3: Diagonalize each l-channel's radial Fock operator
    do l = 0, Grid_Ylm_lmax
      call DiagonalizerList(l + 1) % e % Diagonalize(time, .true.)
    end do

    ! Step 4: Map eigenvectors to orbitals using (n,l,m)
    do j = 1, nOS / 2
      n = OrbsInit_Ylm_HydrogenLike_n(j)
      l = OrbsInit_Ylm_HydrogenLike_l(j)
      m = OrbsInit_Ylm_HydrogenLike_m(j)
      nr = n - l - 1  ! Radial quantum number (node count)

      orbsRaw(:, j) = 0.0_R64
      call Grid_Ylm_SetLmComponent(orbsRaw(:, j), l, m, DiagonalizerList(l + 1) % e % evecs(:, nr + 1))
    end do

    ! Step 5: The radial eigenvectors are normalized in the Euclidean metric;
    ! alignment and mixing require orthonormality with respect to the grid
    ! metric (Grid_InnerProduct), then gauge-align the new orbitals with the
    ! old ones (removes the arbitrary phase of the radial eigenvectors),
    ! then mix and copy spins
    call Grid_Orthonormalize(orbsRaw)
    call Orbs_AlignOnReference(orbsRaw, orbs(:, 1:nOS / 2))

    if (mixTarget .eq. "orbitals") then
      orbsRawFlat(1:nG * (nOS / 2)) => orbsRaw
      call Mixing_Mix(state(1:nG * (nOS / 2)), orbsRawFlat)
    else
      do j = 1, nOS / 2
        orbs(:, j) = orbsRaw(:, j)
      end do
    end if

    do j = 1, nOS / 2
      orbs(:, nOS / 2 + j) = orbs(:, j)  ! Copy spin-up → spin-down
    end do

    ! Step 6: Gram–Schmidt orthonormalize
    call Orbs_Orthonormalize(orbs)

  end subroutine

end submodule
