! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Method_Mb_OrbBased) S_Method_Mb_OrbBased
  !-----------------------------------------------------------------------------
  ! Orbital-based many-body method implementation.
  !
  ! Provides:
  !   - Orbital index mappings per body type
  !   - Default implementations for operator application (kinetic, potential,
  !     interaction, mean-field, correlation)
  !   - Hamiltonian matrix construction (h1, h2)
  !   - Energy calculation via RDM contraction
  !
  ! Delegates to specific methods: TDCI, TDHX, MCTDHX, TD-2RDM, TD-2RDM-StaticOrbs
  !-----------------------------------------------------------------------------

  implicit none

  !=============================================================================
  ! local data
  !=============================================================================

  !> Target group (m, body type) of the definite-m node path: the columns
  !> (pair colPair, orbital colOrb) whose product V_pair * orbital carries
  !> this m, the orbitals of the group, the offset of its columns in the h2
  !> columns of all groups, and the weights (columns x orbitals) of the
  !> correlation sum (refilled at every call).
  type :: T_TargetGroup
    integer(I32), allocatable :: colPair(:)
    integer(I32), allocatable :: colOrb(:)
    integer(I32), allocatable :: orbs(:)
    integer(I32) :: colOffset = 0
    complex(R64), allocatable :: weights(:, :)
  end type

  !> Whether the definite-m node path is bound (Method_Mb_OrbBased_Fabricate).
  logical :: ylmNodesQ = .false.

  !> Skeleton of the definite-m node path (SetupYlmNodes, built once from
  !> Method_Mb_OrbBased_mOfOrb in Method_Mb_OrbBased_Setup): ordered pairs
  !> (k2, l2) of equal body type within the |m| cap of the potentials, every
  !> unordered pair (k2 >= l2) followed by its swap, with M = m_k2 - m_l2 and
  !> the index of the stored potential (negative for the swap, which uses the
  !> conjugate); the target groups; the quadrature weights of the full
  !> measure on the node points.
  integer(I32), allocatable :: nodePairK2(:), nodePairL2(:), nodePairM(:), nodePairPot(:)
  integer(I32) :: nodeNPots = 0
  type(T_TargetGroup), allocatable :: nodeGroups(:)
  real(R64), allocatable :: nodeWeights(:)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Mb_OrbBased_Fabricate
    !---------------------------------------------------------------------------
    ! Initializes orbital-based infrastructure and dispatches to concrete method.
    !
    ! JSON keys read:
    !   - method.mb.orbBased.nOrbs:                 Orbitals per body type
    !   - method.mb.orbBased.regularizationParameter: Numerical stability (default 1e-10)
    !   - method.mb.orbBased.exactFieldFreeQ:       Exact field-free one-body
    !                                               propagation in split steps (default false)
    !   - method.mb.orbBased.exactFieldFreeGpuQ:    cuBLAS variant of the exact
    !                                               propagation (default false)
    !   - method.mb.orbBased.exactFieldFreeCut:     Spectral cut energy of h0
    !                                               (default 0 = off; mctdhx only)
    !   - method.mb.orbBased.mOfOrbs:               Declared m of every spin
    !                                               orbital (Ylm grids, optional)
    !
    ! Binds default operator-application procedures, then branches to the
    ! selected method (TDCI, TDHX, MCTDHX, TD-2RDM, or TD-2RDM-StaticOrbs).
    !---------------------------------------------------------------------------
    use M_Utils_Json
    use M_Utils_Say
    use M_Utils_NoOpProcedures
    use M_Grid_Ylm
    use M_Method
    use M_Method_Mb
    use M_Method_Mb_OrbBased_Tdci
    use M_Method_Mb_OrbBased_Tdhx
    use M_Method_Mb_OrbBased_Mctdhx
    use M_Method_Mb_OrbBased_Td2rdm
    use M_Method_Mb_OrbBased_Td2rdmStaticOrbs

    integer(I32) :: i, ibt, index, counter

    call Say_Fabricate("method.mb.orbBased")

    !------------------------------------
    ! parse orbital configuration
    !------------------------------------

    allocate (Method_Mb_OrbBased_nOrbs(Method_Mb_nBodyTypes))
    Method_Mb_OrbBased_nOrbs = Json_Get("method.mb.orbBased.nOrbs", [(0, i=1, Method_Mb_nBodyTypes)])

    Method_Mb_OrbBased_nOrbsSum = sum(Method_Mb_OrbBased_nOrbs)

    !------------------------------------
    ! build orbital index mappings
    !------------------------------------

    allocate (Method_Mb_OrbBased_bodyTypeOfOrb(Method_Mb_OrbBased_nOrbsSum))
    allocate (Method_Mb_OrbBased_nOrbsStart(Method_Mb_nBodyTypes))
    allocate (Method_Mb_OrbBased_nOrbsEnd(Method_Mb_nBodyTypes))

    counter = 0
    do ibt = 1, Method_Mb_nBodyTypes
      Method_Mb_OrbBased_nOrbsStart(ibt) = counter + 1
      do index = 1, Method_Mb_OrbBased_nOrbs(ibt)
        counter = counter + 1
        Method_Mb_OrbBased_bodyTypeOfOrb(counter) = ibt
      end do
      Method_Mb_OrbBased_nOrbsEnd(ibt) = counter
    end do

    Method_Mb_OrbBased_regularizationParameter = Json_Get("method.mb.orbBased.regularizationParameter", 1e-10_R64)

    !------------------------------------
    ! declared m of the orbitals (Ylm grids)
    !------------------------------------

    Method_Mb_OrbBased_definiteMQ = Json_GetExistence("method.mb.orbBased.mOfOrbs")
    if (Method_Mb_OrbBased_definiteMQ) then
      if (.not. Json_GetExistence("grid.ylm")) error stop "method.mb.orbBased.mOfOrbs requires a Ylm grid"
      Method_Mb_OrbBased_mOfOrb = Json_Get("method.mb.orbBased.mOfOrbs", [(0, i=1, Method_Mb_OrbBased_nOrbsSum)])
      if (size(Method_Mb_OrbBased_mOfOrb) .ne. Method_Mb_OrbBased_nOrbsSum) &
        error stop "method.mb.orbBased.mOfOrbs needs one m per spin orbital (sum of method.mb.orbBased.nOrbs)"
      if (any(abs(Method_Mb_OrbBased_mOfOrb) > Grid_Ylm_mmax)) &
        error stop "method.mb.orbBased.mOfOrbs: |m| exceeds grid.ylm.mmax"
    end if

    !------------------------------------
    ! bind default procedure pointers
    !------------------------------------

    Method_GetEnergy => GetEnergy
    Method_Mb_OrbBased_FillH1 => FillH1
    Method_Mb_OrbBased_FillH2 => FillH2
    Method_Mb_OrbBased_ApplyKineticOp => ApplyKineticOp
    Method_Mb_OrbBased_ApplyCorrelationOp => ApplyCorrelationOp
    Method_Mb_OrbBased_ApplyHartreeFockOp => ApplyHartreeFockOp
    Method_Mb_OrbBased_ApplySingleBodyOp => ApplySingleBodyOp

    ! Declared m: correlation operator and h2 by the definite-m node quadrature
    ylmNodesQ = Method_Mb_OrbBased_definiteMQ .and. YlmNodesFastPathQ()
    if (ylmNodesQ) then
      Method_Mb_OrbBased_ApplyCorrelationOp => ApplyCorrelationOpYlmNodes
      Method_Mb_OrbBased_FillH2 => FillH2YlmNodes
      write (*, '(A)') " method.mb.orbBased.mOfOrbs: correlation operator and h2 by the definite-m node quadrature"
    else if (Method_Mb_OrbBased_definiteMQ) then
      write (*, '(A)') " method.mb.orbBased.mOfOrbs: no node quadrature (needs a non-ECS grid and a &
                       &body-type-independent Ylm interaction), generic correlation operator"
    else if (Json_GetExistence("grid.ylm")) then
      write (*, '(A)') " method.mb.orbBased.mOfOrbs not declared: generic correlation operator and h2"
    end if

    ! Conditionally bind potential/gauge/interaction if configured
    if (Json_GetExistence("sysPotential")) then
      Method_Mb_OrbBased_ApplyPotentialOp => ApplyPotentialOp
    else
      Method_Mb_OrbBased_ApplyPotentialOp => NoOpProcedures_ApplyPotentialOp
    end if

    if (Json_GetExistence("sysGauge")) then
      Method_Mb_OrbBased_ApplyGaugeOp => ApplyGaugeOp
    else
      Method_Mb_OrbBased_ApplyGaugeOp => NoOpProcedures_ApplyPotentialOp
    end if

    if (Json_GetExistence("sysInteraction")) then
      Method_Mb_OrbBased_ApplyInteractionOp => ApplyInteractionOp
    else
      Method_Mb_OrbBased_ApplyInteractionOp => NoOpProcedures_ApplyInteractionOp
    end if

    ! Exact field-free one-body propagation: the linear split-step integrator
    ! then only propagates the gauge coupling; the methods with a split-step
    ! linear part bind Method_PropagateExactPart
    Method_Mb_OrbBased_exactFieldFreeQ = Json_Get("method.mb.orbBased.exactFieldFreeQ", .false.)
    if (Method_Mb_OrbBased_exactFieldFreeQ) then
      Method_Mb_OrbBased_ApplyLinOp => ApplyGaugeOpLin
    else
      Method_Mb_OrbBased_ApplyLinOp => ApplySingleBodyOpLin
    end if

    ! Spectral cut of h0 (needs the eigenvectors of the exact part); the
    ! orbital derivatives of MCTDHX are projected onto the kept eigenvectors
    Method_Mb_OrbBased_exactFieldFreeCut = Json_Get("method.mb.orbBased.exactFieldFreeCut", 0.0_R64)
    if (Method_Mb_OrbBased_exactFieldFreeCut .ne. 0.0_R64) then
      if (.not. Method_Mb_OrbBased_exactFieldFreeQ) &
        error stop "method.mb.orbBased.exactFieldFreeCut requires method.mb.orbBased.exactFieldFreeQ"
      if (.not. Json_GetExistence("method.mb.orbBased.mctdhx")) &
        error stop "method.mb.orbBased.exactFieldFreeCut is only implemented for method.mb.orbBased.mctdhx"
      if (Method_Mb_OrbBased_exactFieldFreeCut < 0.0_R64) &
        error stop "method.mb.orbBased.exactFieldFreeCut must be positive"
      Method_Mb_OrbBased_ApplySpectralCut => Method_Mb_OrbBased_ProjectBelowFieldFreeCut
    else
      Method_Mb_OrbBased_ApplySpectralCut => NoOpProcedures_ApplyAbsorber
    end if

    !------------------------------------
    ! branch: method selection
    !------------------------------------

    if (Json_GetExistence("method.mb.orbBased.tdci")) then
      call Method_Mb_OrbBased_Tdci_Fabricate

    else if (Json_GetExistence("method.mb.orbBased.tdhx")) then
      call Method_Mb_OrbBased_Tdhx_Fabricate

    else if (Json_GetExistence("method.mb.orbBased.mctdhx")) then
      call Method_Mb_OrbBased_Mctdhx_Fabricate

    else if (Json_GetExistence("method.mb.orbBased.td2rdm")) then
      call Method_Mb_OrbBased_Td2rdm_Fabricate

    else if (Json_GetExistence("method.mb.orbBased.td2rdmStaticOrbs")) then
      call Method_Mb_OrbBased_Td2rdmStaticOrbs_Fabricate

    else
      error stop "method.mb.orbBased is missing one of: tdci, tdhx, mctdhx, td2rdm, td2rdmStaticOrbs"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Mb_OrbBased_Setup

    if (ylmNodesQ) call SetupYlmNodes

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Mb_OrbBased_FillSpatialRdm1(rdm1Spatial, rdm1)
    use M_Orbs

    complex(R64), intent(out), allocatable :: rdm1Spatial(:, :)
    complex(R64), intent(in), contiguous   :: rdm1(:, :)

    integer(I32) :: nO, nOS

    nOS = Orbs_nOrbsInState
    nO = Method_Mb_OrbBased_nOrbsSum

    allocate (rdm1Spatial(nOS, nOS))

    if (nO .eq. nOS) then
      ! Unrestricted or single-species: direct copy
      rdm1Spatial(:, :) = rdm1(1:nOS, 1:nOS)
    else if (nO .eq. 2 * nOS) then
      ! Restricted: sum over spin blocks
      rdm1Spatial(:, :) = rdm1(1:nOS, 1:nOS) + rdm1(nOS + 1:, nOS + 1:)
    else
      error stop "invalid nO to nOS relation in Method_Mb_OrbBased_FillSpatialRdm1"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Mb_OrbBased_FillInverseRdm1(invRdm1, rdm1)
    !---------------------------------------------------------------------------
    ! Regularized inverse of the spatial 1-RDM (spin-summed for restricted):
    !   M = U n/(n^2 + eps) U^dagger  (natural orbitals U, occupations n),
    ! eps = Method_Mb_OrbBased_regularizationParameter.
    !---------------------------------------------------------------------------
    use M_Utils_RdmDiagonalize
    use M_Orbs

    complex(R64), intent(out), allocatable :: invRdm1(:, :)
    complex(R64), intent(in), contiguous   :: rdm1(:, :)

    real(R64), allocatable :: natocc(:)
    complex(R64), allocatable :: natorb(:, :)
    complex(R64), allocatable :: rdm1Spatial(:, :)

    integer(I32) :: nOS
    integer(I32) :: j1, k1
    real(R64) :: reg

    nOS = Orbs_nOrbsInState

    call Method_Mb_OrbBased_FillSpatialRdm1(rdm1Spatial, rdm1)
    allocate (invRdm1(nOS, nOS))

    !------------------------------------
    ! Diagonalize 1-RDM → natural orbitals
    !------------------------------------

    call RdmDiagonalize_Rdm1(natocc, natorb, rdm1Spatial, blockOf_=Orbs_blockOfOrb)

    reg = Method_Mb_OrbBased_regularizationParameter

    ! The transform to natural orbitals, regularized inversion, and
    ! back-transform combine into a single matrix acting on index 3:
    !   M_{k₁,j₁} = Σ_j U_{k₁,j} · n_j/(n_j² + ε) · U†_{j,j₁}
    do j1 = 1, nOS
      do k1 = 1, nOS
        invRdm1(k1, j1) = sum(natorb(k1, :) * (natocc(:) / (natocc(:)**2 + reg)) * conjg(natorb(j1, :)))
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Mb_OrbBased_SaveRdm2Blocks(rdm2)
    use M_Utils_DataStorage
    use M_Method_Mb

    complex(R64), intent(in), contiguous :: rdm2(:, :, :, :)

    complex(R64), allocatable :: rdm2Bt(:, :, :, :)
    character(len=256) :: filename
    integer(I32)      :: ibt1, ibt2

    do ibt1 = 1, Method_Mb_nBodyTypes
      do ibt2 = 1, Method_Mb_nBodyTypes

        write (filename, '("rdm2Bt",I2.2,"_",I2.2,".in")') ibt1, ibt2

        rdm2Bt = rdm2(Method_Mb_OrbBased_nOrbsStart(ibt1):Method_Mb_OrbBased_nOrbsEnd(ibt1), &
                      Method_Mb_OrbBased_nOrbsStart(ibt2):Method_Mb_OrbBased_nOrbsEnd(ibt2), &
                      Method_Mb_OrbBased_nOrbsStart(ibt1):Method_Mb_OrbBased_nOrbsEnd(ibt1), &
                      Method_Mb_OrbBased_nOrbsStart(ibt2):Method_Mb_OrbBased_nOrbsEnd(ibt2))

        call SaveData(rdm2Bt, trim(filename), storage_size(rdm2Bt))

      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function YlmNodesFastPathQ() result(res)
    !---------------------------------------------------------------------------
    ! Whether the definite-m node quadrature (ApplyCorrelationOpYlmNodes,
    ! FillH2YlmNodes) applies to orbitals of declared m: Hermitian (non-ECS)
    ! Ylm grid, body-type-independent Ylm interaction with potentials within
    ! the angular quadrature tables of the grid. Requires Grid and
    ! SysInteraction to be fabricated.
    !---------------------------------------------------------------------------
    use M_Utils_Json
    use M_Grid_Ylm
    use M_SysInteraction
    use M_SysInteraction_Ylm

    logical :: res

    res = .false.
    if (.not. Json_GetExistence("grid.ylm")) return
    if (Json_GetExistence("grid.ylm.fedvrEcs")) return
    if (.not. Json_GetExistence("sysInteraction.ylm")) return
    if (.not. SysInteraction_bodyTypeIndependentQ) return
    if (SysInteraction_Ylm_lmax > 2 * Grid_Ylm_lmax + 2) return
    if (SysInteraction_Ylm_mmax > min(2 * Grid_Ylm_mmax, 2 * Grid_Ylm_lmax + 2)) return

    res = .true.

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function GetEnergy(time) result(res)
    !---------------------------------------------------------------------------
    ! Computes the total energy via RDM contraction:
    !   E = Tr[ρ¹ h¹] + ½ Tr[ρ² h²]
    !
    ! Uses the current state to extract RDMs and compute Hamiltonian matrices,
    ! then contracts them via RdmObservables_Energy.
    !---------------------------------------------------------------------------
    use M_Utils_RdmObservables
    use M_Method
    use M_Orbs

    real(R64)                :: res
    real(R64), intent(in)    :: time

    ! Extract RDMs from current state
    call Method_Mb_OrbBased_FillRdm1(Method_Mb_OrbBased_rdm1, Method_state)
    call Method_Mb_OrbBased_FillRdm2(Method_Mb_OrbBased_rdm2, Method_state)

    ! Build Hamiltonian matrices in orbital basis
    call Method_Mb_OrbBased_FillH1(Method_Mb_OrbBased_h1, Orbs_orbs, time)
    call Method_Mb_OrbBased_FillH2(Method_Mb_OrbBased_h2, Orbs_orbs, time)

    ! Contract: E = Tr[ρ¹ h¹] + ½ Tr[ρ² h²]
    res = RdmObservables_Energy(Method_Mb_OrbBased_rdm1, &
                                Method_Mb_OrbBased_h1, &
                                Method_Mb_OrbBased_rdm2, &
                                Method_Mb_OrbBased_h2)

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillH1(h1, orbs, time)
    !---------------------------------------------------------------------------
    ! Builds the one-body Hamiltonian matrix in the orbital basis:
    !   h¹_{ij} = ⟨φ_i| T̂ + V̂_ext + Ĥ_A |φ_j⟩
    ! where Ĥ_A is the velocity-gauge coupling (zero unless sysGauge is set).
    !
    ! Only computes elements where both orbitals belong to the same body type.
    ! For restricted calculations, duplicates the computed block for spin-down.
    !---------------------------------------------------------------------------
    use M_Method
    use M_Orbs
    use M_Grid
    use M_SysKinetic
    use M_SysPotential

    complex(R64), intent(out), allocatable :: h1(:, :)
    complex(R64), intent(in), contiguous                :: orbs(:, :)
    real(R64), intent(in)                  :: time

    integer(I32) :: nOS, nO, nG
    integer(I32) :: i1, j1, i1bt, j1bt
    complex(R64), allocatable :: orbTmps(:, :)

    nG = Grid_nPoints
    nOS = Orbs_nOrbsInState
    nO = Method_Mb_OrbBased_nOrbsSum

    if (.not. allocated(h1)) allocate (h1(nO, nO))
    h1(:, :) = 0.0_R64

    ! Apply (T̂ + V̂ + Ĥ_A) to all orbitals
    allocate (orbTmps(nG, nOS))
    orbTmps(:, :) = 0.0_R64

    call Method_Mb_OrbBased_ApplyKineticOp(orbTmps, orbs, time)
    call Method_Mb_OrbBased_ApplyPotentialOp(orbTmps, orbs, time)
    call Method_Mb_OrbBased_ApplyGaugeOp(orbTmps, orbs, time)

    ! Compute matrix elements via inner products
    !$omp parallel do default(shared) private(i1, j1, i1bt, j1bt) collapse(2)
    do j1 = 1, nOS
      do i1 = 1, nOS

        i1bt = Method_Mb_OrbBased_bodyTypeOfOrb(i1)
        j1bt = Method_Mb_OrbBased_bodyTypeOfOrb(j1)

        ! Skip cross-body-type elements (vanish by orthogonality)
        if (i1bt .ne. j1bt) cycle

        h1(i1, j1) = Grid_InnerProduct(orbs(:, i1), orbTmps(:, j1))

      end do
    end do
    !$omp end parallel do

    ! For restricted calculations, duplicate for second spin block
    if (Orbs_restrictedQ) then
      h1(nOS + 1:, nOS + 1:) = h1(:nOS, :nOS)
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillH2(h2, orbs, time)
    !---------------------------------------------------------------------------
    ! Builds the two-body Hamiltonian tensor in the orbital basis:
    !   h²_{i₁i₂j₁j₂} = ½ ∫∫ φ*_{i₁}(r₁) φ*_{i₂}(r₂) W(r₁,r₂) φ_{j₁}(r₁) φ_{j₂}(r₂) dr₁dr₂
    !
    ! The factor ½ accounts for double-counting in the Hamiltonian.
    ! Only computes elements where orbitals within each pair share body type.
    ! For restricted calculations, duplicates blocks for spin combinations.
    !---------------------------------------------------------------------------
    use M_Method
    use M_SysInteraction
    use M_Orbs
    use M_Grid

    complex(R64), intent(out), allocatable :: h2(:, :, :, :)
    complex(R64), intent(in), contiguous   :: orbs(:, :)
    real(R64), intent(in)                  :: time

    integer(I32) :: nG, nOS, nO
    integer(I32) :: i2, j2
    integer(I32) :: i2bt, j2bt
    complex(R64), allocatable :: orbTmps(:, :), orbTmpsSwap(:, :)

    nG = Grid_nPoints
    nOS = Orbs_nOrbsInState
    nO = Method_Mb_OrbBased_nOrbsSum

    if (.not. allocated(h2)) allocate (h2(nO, nO, nO, nO))
    h2(:, :, :, :) = 0.0_R64

    allocate (orbTmps, mold=orbs)
    allocate (orbTmpsSwap, mold=orbs)

    ! Loop over unordered second-particle pairs (i2 >= j2); the swapped
    ! pair (j2,i2) reuses the potential via the conjugation symmetry
    do j2 = 1, nOS
      do i2 = j2, nOS
        i2bt = Method_Mb_OrbBased_bodyTypeOfOrb(i2)
        j2bt = Method_Mb_OrbBased_bodyTypeOfOrb(j2)

        ! Skip cross-body-type pairs for second particle
        if (i2bt .ne. j2bt) cycle

        orbTmps(:, :) = 0.0_R64

        if (i2 .eq. j2) then
          call Method_Mb_OrbBased_ApplyInteractionOp(orbTmps, orbs, i2, j2, time)

          call FillH2Block(h2, orbs, orbTmps, i2, j2)
        else
          orbTmpsSwap(:, :) = 0.0_R64
          call ApplyInteractionOpPair(orbTmps, orbTmpsSwap, orbs, i2, j2, time)

          call FillH2Block(h2, orbs, orbTmps, i2, j2)
          call FillH2Block(h2, orbs, orbTmpsSwap, j2, i2)
        end if
      end do
    end do

    ! For restricted calculations, duplicate blocks for spin combinations
    ! (αα, αβ, βα, ββ all have same spatial integrals)
    if (Orbs_restrictedQ) then
      h2(nOS + 1:, :nOS, nOS + 1:, :nOS) = h2(:nOS, :nOS, :nOS, :nOS)
      h2(:nOS, nOS + 1:, :nOS, nOS + 1:) = h2(:nOS, :nOS, :nOS, :nOS)
      h2(nOS + 1:, nOS + 1:, nOS + 1:, nOS + 1:) = h2(:nOS, :nOS, :nOS, :nOS)
    end if
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillH2Block(h2, orbs, orbTmps, i2, j2)
    !---------------------------------------------------------------------------
    ! Fills the first-particle matrix elements for a fixed second-particle
    ! source pair (i2,j2):
    !   h²_{i₁,i2,j₁,j2} = ½ ⟨φ_i₁| V_(i2,j2) |φ_j₁⟩
    !---------------------------------------------------------------------------
    use M_Grid
    use M_Orbs

    complex(R64), intent(inout), contiguous :: h2(:, :, :, :)
    complex(R64), intent(in), contiguous    :: orbs(:, :)
    complex(R64), intent(in), contiguous    :: orbTmps(:, :)
    integer(I32), intent(in)                :: i2
    integer(I32), intent(in)                :: j2

    integer(I32) :: i1, j1, i1bt, j1bt, nOS

    nOS = Orbs_nOrbsInState

    !$omp parallel do default(shared) private(i1, j1, i1bt, j1bt) collapse(2)
    do j1 = 1, nOS
      do i1 = 1, nOS
        i1bt = Method_Mb_OrbBased_bodyTypeOfOrb(i1)
        j1bt = Method_Mb_OrbBased_bodyTypeOfOrb(j1)

        ! Skip cross-body-type pairs for first particle
        if (i1bt .ne. j1bt) cycle

        h2(i1, i2, j1, j2) = 0.5_R64 * Grid_InnerProduct(orbs(:, i1), orbTmps(:, j1))
      end do
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ApplyKineticOp(dOrbs, orbs, time)
    !---------------------------------------------------------------------------
    ! Applies the kinetic energy operator to all orbitals:
    !   dOrbs(:,i) += T̂ · orbs(:,i)  for i = 1, ..., nOrbsInState
    !
    ! The kinetic operator is body-type-dependent (different masses possible).
    !---------------------------------------------------------------------------
    use M_Method
    use M_Grid
    use M_Orbs
    use M_SysKinetic

    complex(R64), intent(inout), contiguous :: dOrbs(:, :)
    complex(R64), intent(in), contiguous    :: orbs(:, :)
    real(R64), intent(in)                  :: time

    complex(R64), allocatable :: orbTmp(:)

    integer(I32) :: nG, nOS
    integer(I32) :: i1, ibt

    nG = Grid_nPoints
    nOS = Orbs_nOrbsInState

    !$omp parallel default(shared) private(i1, ibt, orbTmp)
    allocate (orbTmp(nG))

    !$omp do
    do i1 = 1, nOS

      ibt = Method_Mb_OrbBased_bodyTypeOfOrb(i1)

      orbTmp(:) = 0.0_R64
      call SysKinetic_MultiplyWithKineticOp(orbTmp, orbs(:, i1), time, ibt)

      dOrbs(:, i1) = dOrbs(:, i1) + orbTmp(:)

    end do
    !$omp end do

    deallocate (orbTmp)
    !$omp end parallel

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ApplyGaugeOp(dOrbs, orbs, time)
    !---------------------------------------------------------------------------
    ! Applies the velocity-gauge coupling operator to all orbitals:
    !   dOrbs(:,i) += (A(t)·p̂/m + A²(t)/(2m)) · orbs(:,i)
    !
    ! The gauge operator is body-type-dependent through the mass, mirroring
    ! the kinetic operator.
    !---------------------------------------------------------------------------
    use M_Method
    use M_Grid
    use M_Orbs
    use M_SysGauge

    complex(R64), intent(inout), contiguous :: dOrbs(:, :)
    complex(R64), intent(in), contiguous    :: orbs(:, :)
    real(R64), intent(in)                  :: time

    complex(R64), allocatable :: orbTmp(:)

    integer(I32) :: nG, nOS
    integer(I32) :: i1, ibt

    nG = Grid_nPoints
    nOS = Orbs_nOrbsInState

    !$omp parallel default(shared) private(i1, ibt, orbTmp)
    allocate (orbTmp(nG))

    !$omp do
    do i1 = 1, nOS

      ibt = Method_Mb_OrbBased_bodyTypeOfOrb(i1)

      orbTmp(:) = 0.0_R64
      call SysGauge_MultiplyWithGaugeOp(orbTmp, orbs(:, i1), time, ibt)

      dOrbs(:, i1) = dOrbs(:, i1) + orbTmp(:)

    end do
    !$omp end do

    deallocate (orbTmp)
    !$omp end parallel

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ApplyPotentialOp(dOrbs, orbs, time)
    !---------------------------------------------------------------------------
    ! Applies the external potential operator to all orbitals:
    !   dOrbs(:,i) += V̂_ext · orbs(:,i)  for i = 1, ..., nOrbsInState
    !
    ! Optimizes by computing the potential once if body-type-independent,
    ! otherwise loops over body types and applies to corresponding orbitals.
    !---------------------------------------------------------------------------
    use M_Grid
    use M_Grid_Ylm
    use M_Orbs
    use M_Method_Mb
    use M_SysPotential

    complex(R64), intent(inout), contiguous :: dOrbs(:, :)
    complex(R64), intent(in), contiguous    :: orbs(:, :)
    real(R64), intent(in) :: time

    integer(I32) :: nG, nOS, nBt
    integer(I32) :: iOrb, bt
    complex(R64), allocatable :: orbTmp(:)
    complex(R64), allocatable :: externalPotential(:)

    nG = Grid_nPoints
    nOS = Orbs_nOrbsInState

    ! Restricted: only the shared spatial set is stored; the body-type-1
    ! potential covers all stored orbitals (spin up and down coincide)
    nBt = Method_Mb_nBodyTypes
    if (Orbs_restrictedQ) nBt = 1

    if (SysPotential_bodyTypeIndependentQ) then
      ! Body-type-independent: compute potential once, apply to all orbitals
      call SysPotential_FillExternalPotential(externalPotential, time)

      !$omp parallel default(shared) private(iOrb, orbTmp)
      allocate (orbTmp(nG))
      !$omp do
      do iOrb = 1, nOS
        call SysPotential_MultiplyWithExternalPotential(orbTmp, externalPotential, orbs(:, iOrb))
        dOrbs(:, iOrb) = dOrbs(:, iOrb) + orbTmp(:)
      end do
      !$omp end do
      deallocate (orbTmp)
      !$omp end parallel
    else
      ! Body-type-dependent: loop over types and apply to corresponding orbitals
      do bt = 1, nBt
        call SysPotential_FillExternalPotential(externalPotential, time, bt_=bt)

        !$omp parallel default(shared) private(iOrb, orbTmp)
        allocate (orbTmp(nG))
        !$omp do
        do iOrb = Method_Mb_OrbBased_nOrbsStart(bt), Method_Mb_OrbBased_nOrbsEnd(bt)
          call SysPotential_MultiplyWithExternalPotential(orbTmp, externalPotential, orbs(:, iOrb))
          dOrbs(:, iOrb) = dOrbs(:, iOrb) + orbTmp(:)
        end do
        !$omp end do
        deallocate (orbTmp)
        !$omp end parallel
      end do
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ApplyInteractionOp(dOrbs, orbs, i2, j2, time)
    !---------------------------------------------------------------------------
    ! Applies the two-body interaction operator for a given (i2,j2) source pair:
    !   dOrbs(:,i1) += ∫ W(r₁,r₂) φ*_{i2}(r₂) φ_{j2}(r₂) dr₂ · orbs(:,i1)
    !
    ! This computes the "direct" interaction contribution from the density
    ! |φ_{i2}|² (when i2=j2) or the transition density φ*_{i2}φ_{j2}.
    !---------------------------------------------------------------------------
    use M_Grid
    use M_Grid_Ylm
    use M_Method_Mb
    use M_SysInteraction
    use M_Orbs

    complex(R64), intent(inout), contiguous :: dOrbs(:, :)
    complex(R64), intent(in), contiguous    :: orbs(:, :)
    integer(I32), intent(in)               :: i2
    integer(I32), intent(in)               :: j2
    real(R64), intent(in)                  :: time

    integer(I32) :: nOS, nBt
    integer(I32) :: bt1, bt2
    complex(R64), allocatable :: src(:), interactionPotential(:)

    nOS = Orbs_nOrbsInState

    ! Restricted: only the shared spatial set is stored; the body-type-1
    ! interaction covers all stored orbitals (spin up and down coincide)
    nBt = Method_Mb_nBodyTypes
    if (Orbs_restrictedQ) nBt = 1

    ! Build source density: ρ(r₂) = φ*_{i2}(r₂) φ_{j2}(r₂)
    call SysInteraction_FillInteractionSrc(src, orbs(:, i2), orbs(:, j2))

    bt2 = Method_Mb_OrbBased_bodyTypeOfOrb(i2)

    if (SysInteraction_bodyTypeIndependentQ) then
      ! Body-type-independent: compute potential once, apply to all orbitals
      call SysInteraction_FillInteractionPotential(interactionPotential, src, time, -1, -1)
      call MultiplyPotentialOnOrbRange(dOrbs, orbs, interactionPotential, 1, nOS)
    else
      ! Body-type-dependent: loop over target body types
      do bt1 = 1, nBt
        call SysInteraction_FillInteractionPotential(interactionPotential, src, time, bt1, bt2)
        call MultiplyPotentialOnOrbRange(dOrbs, orbs, interactionPotential, &
                                         Method_Mb_OrbBased_nOrbsStart(bt1), Method_Mb_OrbBased_nOrbsEnd(bt1))
      end do
    end if

    deallocate (src, interactionPotential)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine MultiplyPotentialOnOrbRange(dOrbs, orbs, interactionPotential, iOrbStart, iOrbEnd)
    !---------------------------------------------------------------------------
    ! Accumulates dOrbs(:,i) += V · orbs(:,i) for i = iOrbStart, ..., iOrbEnd.
    !---------------------------------------------------------------------------
    use M_Grid
    use M_SysInteraction

    complex(R64), intent(inout), contiguous :: dOrbs(:, :)
    complex(R64), intent(in), contiguous    :: orbs(:, :)
    complex(R64), intent(in), contiguous    :: interactionPotential(:)
    integer(I32), intent(in)                :: iOrbStart
    integer(I32), intent(in)                :: iOrbEnd

    integer(I32) :: iOrb
    complex(R64), allocatable :: orbTmp(:)

    !$omp parallel default(shared) private(iOrb, orbTmp)
    allocate (orbTmp(Grid_nPoints))
    !$omp do
    do iOrb = iOrbStart, iOrbEnd
      call SysInteraction_MultiplyWithInteractionPotential(orbTmp, interactionPotential, orbs(:, iOrb))
      dOrbs(:, iOrb) = dOrbs(:, iOrb) + orbTmp(:)
    end do
    !$omp end do
    deallocate (orbTmp)
    !$omp end parallel

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ApplyInteractionOpPair(dOrbsA, dOrbsB, orbs, i2, j2, time)
    !---------------------------------------------------------------------------
    ! Accumulates the interaction operator for source pair (i2,j2) into dOrbsA
    ! and for the swapped pair (j2,i2) into dOrbsB, assuming i2 ≠ j2.
    !
    ! When the interaction backend guarantees the conjugation symmetry
    !   V_(j2,i2) = conjg(V_(i2,j2)),
    ! the second potential solve is replaced by a cheap conjugation, halving
    ! the number of solves over unordered orbital pairs.
    !---------------------------------------------------------------------------
    use M_Grid
    use M_Method_Mb
    use M_Orbs
    use M_SysInteraction

    complex(R64), intent(inout), contiguous :: dOrbsA(:, :)
    complex(R64), intent(inout), contiguous :: dOrbsB(:, :)
    complex(R64), intent(in), contiguous    :: orbs(:, :)
    integer(I32), intent(in)                :: i2
    integer(I32), intent(in)                :: j2
    real(R64), intent(in)                   :: time

    integer(I32) :: nOS, nBt, bt1, bt2
    complex(R64), allocatable :: src(:), interactionPotential(:), interactionPotentialConjg(:)

    if (.not. SysInteraction_conjSymmetricQ) then
      ! No usable symmetry (or no interaction configured): two full applications
      call Method_Mb_OrbBased_ApplyInteractionOp(dOrbsA, orbs, i2, j2, time)
      call Method_Mb_OrbBased_ApplyInteractionOp(dOrbsB, orbs, j2, i2, time)
      return
    end if

    nOS = Orbs_nOrbsInState

    ! Restricted: only the shared spatial set is stored; the body-type-1
    ! interaction covers all stored orbitals (spin up and down coincide)
    nBt = Method_Mb_nBodyTypes
    if (Orbs_restrictedQ) nBt = 1

    call SysInteraction_FillInteractionSrc(src, orbs(:, i2), orbs(:, j2))

    bt2 = Method_Mb_OrbBased_bodyTypeOfOrb(i2)

    if (SysInteraction_bodyTypeIndependentQ) then
      call SysInteraction_FillInteractionPotential(interactionPotential, src, time, -1, -1)

      allocate (interactionPotentialConjg(size(interactionPotential)))
      call SysInteraction_ConjugateInteractionPotential(interactionPotentialConjg, interactionPotential)

      call MultiplyPotentialOnOrbRange(dOrbsA, orbs, interactionPotential, 1, nOS)
      call MultiplyPotentialOnOrbRange(dOrbsB, orbs, interactionPotentialConjg, 1, nOS)
    else
      do bt1 = 1, nBt
        call SysInteraction_FillInteractionPotential(interactionPotential, src, time, bt1, bt2)

        if (.not. allocated(interactionPotentialConjg)) then
          allocate (interactionPotentialConjg(size(interactionPotential)))
        end if
        call SysInteraction_ConjugateInteractionPotential(interactionPotentialConjg, interactionPotential)

        call MultiplyPotentialOnOrbRange(dOrbsA, orbs, interactionPotential, &
                                         Method_Mb_OrbBased_nOrbsStart(bt1), Method_Mb_OrbBased_nOrbsEnd(bt1))
        call MultiplyPotentialOnOrbRange(dOrbsB, orbs, interactionPotentialConjg, &
                                         Method_Mb_OrbBased_nOrbsStart(bt1), Method_Mb_OrbBased_nOrbsEnd(bt1))
      end do
    end if

    deallocate (src, interactionPotential, interactionPotentialConjg)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ApplyCorrelationOp(dOrbs, orbs, invRdm1, rdm2, time, h2_)
    !---------------------------------------------------------------------------
    ! Applies the beyond-mean-field correlation operator to orbitals.
    !
    ! This implements the MCTDHX-style correlation contribution arising from
    ! the 2-RDM structure beyond what Hartree-Fock captures. The algorithm:
    !   1. Contract the spatial 2-RDM with the inverse 1-RDM invRdm1
    !      (regularized, see Method_Mb_OrbBased_FillInverseRdm1)
    !   2. Contract with interaction operator to accumulate in dOrbs
    !
    ! Optionally returns the two-body Hamiltonian matrix h2_ if requested.
    !---------------------------------------------------------------------------
    use M_Method
    use M_Grid
    use M_SysInteraction
    use M_Orbs

    complex(R64), intent(inout), contiguous, target :: dOrbs(:, :)
    complex(R64), intent(in), contiguous, target    :: orbs(:, :)
    complex(R64), intent(in), contiguous             :: invRdm1(:, :)
    complex(R64), intent(in), contiguous             :: rdm2(:, :, :, :)
    real(R64), intent(in)                  :: time
    complex(R64), intent(out), allocatable, optional :: h2_(:, :, :, :)

    complex(R64), allocatable :: orbTmps(:, :), orbTmpsSwap(:, :)
    complex(R64), allocatable :: rdm2spatial(:, :, :, :)

    integer(I32) :: nO, nOS
    integer(I32) :: k2, l2
    integer(I32) :: k2bt, l2bt

    nOS = Orbs_nOrbsInState
    nO = Method_Mb_OrbBased_nOrbsSum

    if (present(h2_)) then
      if (.not. allocated(h2_)) allocate (h2_(nO, nO, nO, nO))
      h2_(:, :, :, :) = 0.0_R64
    end if

    call FillCorrelationRdm2(rdm2spatial, invRdm1, rdm2)

    allocate (orbTmps, mold=orbs)
    allocate (orbTmpsSwap, mold=orbs)

    !------------------------------------
    ! Contract with interaction and accumulate
    !------------------------------------

    ! Loop over unordered pairs (k2 >= l2); the interaction applied for the
    ! swapped source pair reuses the potential via the conjugation symmetry
    do l2 = 1, nOS
      do k2 = l2, nOS
        l2bt = Method_Mb_OrbBased_bodyTypeOfOrb(l2)
        k2bt = Method_Mb_OrbBased_bodyTypeOfOrb(k2)

        if (l2bt .ne. k2bt) cycle

        orbTmps(:, :) = 0.0_R64

        if (k2 .eq. l2) then
          call Method_Mb_OrbBased_ApplyInteractionOp(orbTmps, orbs, l2, k2, time)

          call AccumulateCorrelation(dOrbs, orbTmps, rdm2spatial, k2, l2)
          if (present(h2_)) call FillH2CorrelationBlock(h2_, orbs, orbTmps, k2, l2)
        else
          orbTmpsSwap(:, :) = 0.0_R64
          call ApplyInteractionOpPair(orbTmps, orbTmpsSwap, orbs, l2, k2, time)

          call AccumulateCorrelation(dOrbs, orbTmps, rdm2spatial, k2, l2)
          call AccumulateCorrelation(dOrbs, orbTmpsSwap, rdm2spatial, l2, k2)
          if (present(h2_)) then
            call FillH2CorrelationBlock(h2_, orbs, orbTmps, k2, l2)
            call FillH2CorrelationBlock(h2_, orbs, orbTmpsSwap, l2, k2)
          end if
        end if

      end do
    end do

    ! For restricted calculations, duplicate the h2 blocks for the spin
    ! combinations (αα, αβ, βα, ββ share the spatial integrals); the Coeffs
    ! backend consumes h2 in the full spin-orbital basis (matches FillH2)
    if (present(h2_)) then
      if (Orbs_restrictedQ) then
        h2_(nOS + 1:, :nOS, nOS + 1:, :nOS) = h2_(:nOS, :nOS, :nOS, :nOS)
        h2_(:nOS, nOS + 1:, :nOS, nOS + 1:) = h2_(:nOS, :nOS, :nOS, :nOS)
        h2_(nOS + 1:, nOS + 1:, nOS + 1:, nOS + 1:) = h2_(:nOS, :nOS, :nOS, :nOS)
      end if
    end if

    deallocate (rdm2spatial)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ApplyCorrelationOpYlmNodes(dOrbs, orbs, invRdm1, rdm2, time, h2_)
    !---------------------------------------------------------------------------
    ! ApplyCorrelationOp for orbitals of declared m (Method_Mb_OrbBased_mOfOrb:
    ! each orbital lives in the channels of its m, e.g. m = 0 orbitals or the
    ! orbitals of an atom in a linearly polarized field) on a Ylm grid,
    ! evaluated on the Gauss-Legendre nodes in cos(theta)
    ! (Grid_Ylm_ToNodes/FromNodes).
    !
    ! The generic version multiplies every pair potential V_(l2,k2) with every
    ! orbital (one spatial product each) and then combines the results with the
    ! 2-RDM weights,
    !   dOrbs(:,j1) += sum_(k1,k2,l2) V_(l2,k2) * orbs(:,k1) rdm2spatial(k1,k2,j1,l2).
    ! On the nodes all products are pointwise: the orbitals are transformed
    ! once, each pair density/potential once, and dOrbs once at the end. The
    ! h2 elements are node quadratures of conjg(orb) V orb. The phases
    ! e^(i m phi) are bookkept: the pair (l2,k2) carries M = m_k2 - m_l2, and
    ! only the terms with m_k1 + M = m_j1 (correlation) and m_j1 + M = m_i1
    ! (h2) are kept. The dropped h2 terms vanish by the phi integral, the
    ! dropped correlation terms for an M-conserving many-body state (their
    ! 2-RDM weights vanish). The node count of the Grid_Ylm tables makes every
    ! angular integral exact, so the result equals the generic version up to
    ! round-off.
    !
    ! Bound by Method_Mb_OrbBased_Fabricate for declared m on non-ECS Ylm
    ! grids with a body-type-independent Ylm interaction (YlmNodesFastPathQ).
    !---------------------------------------------------------------------------
    use M_Grid
    use M_Grid_Ylm
    use M_Orbs

    complex(R64), intent(inout), contiguous, target :: dOrbs(:, :)
    complex(R64), intent(in), contiguous, target    :: orbs(:, :)
    complex(R64), intent(in), contiguous             :: invRdm1(:, :)
    complex(R64), intent(in), contiguous             :: rdm2(:, :, :, :)
    real(R64), intent(in)                  :: time
    complex(R64), intent(out), allocatable, optional :: h2_(:, :, :, :)

    complex(R64), allocatable :: rdm2spatial(:, :, :, :)
    complex(R64), allocatable, target :: orbNodes(:, :), dOrbNodes(:, :), potNodes(:, :)
    complex(R64), allocatable :: h2Cols(:, :), dOrb(:)
    complex(R64), contiguous, pointer :: nodes(:, :)

    integer(I32) :: nOS, nRad, nNodes, i1

    nOS = Orbs_nOrbsInState
    nRad = Grid_Ylm_nRadial
    nNodes = Grid_Ylm_nNodes

    call FillCorrelationRdm2(rdm2spatial, invRdm1, rdm2)
    call FillTargetGroupWeights(rdm2spatial)

    call FillOrbNodes(orbNodes, orbs)
    call FillPairPotentialsYlmNodes(potNodes, orbNodes, time)

    allocate (dOrbNodes, mold=orbNodes)
    dOrbNodes(:, :) = 0.0_R64

    if (present(h2_)) then
      call AccumulateTargetSumsYlmNodes(orbNodes, potNodes, h2Cols_=h2Cols, dOrbNodes_=dOrbNodes)
      call FillH2FromTargetColumns(h2_, h2Cols)
    else
      call AccumulateTargetSumsYlmNodes(orbNodes, potNodes, dOrbNodes_=dOrbNodes)
    end if

    ! Back to the Ylm representation
    !$omp parallel default(shared) private(i1, dOrb, nodes)
    allocate (dOrb(Grid_nPoints))
    !$omp do
    do i1 = 1, nOS
      nodes(1:nRad, 1:nNodes) => dOrbNodes(:, i1)
      call Grid_Ylm_FromNodes(dOrb, nodes, Grid_Ylm_lmax, Method_Mb_OrbBased_mOfOrb(i1))
      dOrbs(:, i1) = dOrbs(:, i1) + dOrb(:)
    end do
    !$omp end do
    deallocate (dOrb)
    !$omp end parallel

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillH2YlmNodes(h2, orbs, time)
    !---------------------------------------------------------------------------
    ! FillH2 for orbitals of declared m by the node quadrature of
    ! ApplyCorrelationOpYlmNodes (pair potentials and h2 columns only):
    !   h2(i1,l2,k1,k2) = 1/2 <orb_i1| V_(l2,k2) |orb_k1>,  m_i1 = m_k1 + m_k2 - m_l2.
    ! The other elements vanish by the phi integral. Bound together with
    ! ApplyCorrelationOpYlmNodes.
    !---------------------------------------------------------------------------
    complex(R64), intent(out), allocatable :: h2(:, :, :, :)
    complex(R64), intent(in), contiguous   :: orbs(:, :)
    real(R64), intent(in)                  :: time

    complex(R64), allocatable :: orbNodes(:, :), potNodes(:, :), h2Cols(:, :)

    call FillOrbNodes(orbNodes, orbs)
    call FillPairPotentialsYlmNodes(potNodes, orbNodes, time)
    call AccumulateTargetSumsYlmNodes(orbNodes, potNodes, h2Cols_=h2Cols)
    call FillH2FromTargetColumns(h2, h2Cols)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine SetupYlmNodes
    !---------------------------------------------------------------------------
    ! Builds the skeleton of the node path once (Method_Mb_OrbBased_Setup)
    ! from the declared m of the stored orbitals: the pair list, the target groups (m, body type) of the
    ! selection rule m_k + M_iP = m_j with their columns (iP, k) and orbitals
    ! j, and the quadrature weights on the node points.
    !---------------------------------------------------------------------------
    use M_Grid_Ylm
    use M_SysInteraction_Ylm
    use M_Orbs

    integer(I32), allocatable :: mOrb(:), btOrb(:), groupOfOrb(:), nColsOfGroup(:)
    integer(I32) :: nOS, nPairs, nG, iG, iP, i1, k1, k2, l2, iC, offset

    nOS = Orbs_nOrbsInState
    mOrb = Method_Mb_OrbBased_mOfOrb(1:nOS)
    btOrb = Method_Mb_OrbBased_bodyTypeOfOrb(1:nOS)

    ! Ordered pairs (k2, l2) of equal body type within the |m| cap of the
    ! potentials (pairs beyond it are left out, as in the generic operator)
    nPairs = 0
    do l2 = 1, nOS
      do k2 = 1, nOS
        if (btOrb(l2) .ne. btOrb(k2)) cycle
        if (abs(mOrb(k2) - mOrb(l2)) > SysInteraction_Ylm_mmax) cycle
        nPairs = nPairs + 1
      end do
    end do
    allocate (nodePairK2(nPairs), nodePairL2(nPairs), nodePairM(nPairs), nodePairPot(nPairs))

    iP = 0
    nodeNPots = 0
    do l2 = 1, nOS
      do k2 = l2, nOS
        if (btOrb(l2) .ne. btOrb(k2)) cycle
        if (abs(mOrb(k2) - mOrb(l2)) > SysInteraction_Ylm_mmax) cycle
        iP = iP + 1
        nodeNPots = nodeNPots + 1
        nodePairK2(iP) = k2
        nodePairL2(iP) = l2
        nodePairM(iP) = mOrb(k2) - mOrb(l2)
        nodePairPot(iP) = nodeNPots
        if (k2 .ne. l2) then
          iP = iP + 1
          nodePairK2(iP) = l2
          nodePairL2(iP) = k2
          nodePairM(iP) = -nodePairM(iP - 1)
          nodePairPot(iP) = -nodeNPots
        end if
      end do
    end do

    ! Groups of orbitals with equal (m, body type)
    allocate (groupOfOrb(nOS))
    nG = 0
    do i1 = 1, nOS
      groupOfOrb(i1) = 0
      do k1 = 1, i1 - 1
        if (mOrb(k1) .eq. mOrb(i1) .and. btOrb(k1) .eq. btOrb(i1)) then
          groupOfOrb(i1) = groupOfOrb(k1)
          exit
        end if
      end do
      if (groupOfOrb(i1) .eq. 0) then
        nG = nG + 1
        groupOfOrb(i1) = nG
      end if
    end do

    ! Columns (iP, k) of every group: the target of V_iP * orbs(:,k) has
    ! m_k + M_iP and the body type of k; count, then fill
    allocate (nColsOfGroup(nG))
    nColsOfGroup(:) = 0
    do iP = 1, nPairs
      do k1 = 1, nOS
        iG = TargetGroup(k1, iP)
        if (iG > 0) nColsOfGroup(iG) = nColsOfGroup(iG) + 1
      end do
    end do

    allocate (nodeGroups(nG))
    offset = 0
    do iG = 1, nG
      associate (g => nodeGroups(iG))
        g % orbs = pack([(i1, i1=1, nOS)], groupOfOrb .eq. iG)
        allocate (g % colPair(nColsOfGroup(iG)), g % colOrb(nColsOfGroup(iG)))
        allocate (g % weights(nColsOfGroup(iG), size(g % orbs)))
        g % colOffset = offset
        offset = offset + nColsOfGroup(iG)
      end associate
    end do

    nColsOfGroup(:) = 0
    do iP = 1, nPairs
      do k1 = 1, nOS
        iG = TargetGroup(k1, iP)
        if (iG .eq. 0) cycle
        nColsOfGroup(iG) = nColsOfGroup(iG) + 1
        iC = nColsOfGroup(iG)
        nodeGroups(iG) % colPair(iC) = iP
        nodeGroups(iG) % colOrb(iC) = k1
      end do
    end do

    ! Quadrature weights of the full measure on the node points
    nodeWeights = reshape(spread(Grid_Ylm_radialWeights(:), 2, Grid_Ylm_nNodes) * &
                          spread(Grid_Ylm_nodeWeights(:), 1, Grid_Ylm_nRadial), [Grid_Ylm_nRadial * Grid_Ylm_nNodes])

  contains

    !> Group of the target of V_iP * orbs(:,k), 0 if no orbital has its m.
    function TargetGroup(k, iPair) result(res)
      integer(I32), intent(in) :: k
      integer(I32), intent(in) :: iPair
      integer(I32) :: res

      integer(I32) :: j

      res = 0
      do j = 1, nOS
        if (mOrb(j) .eq. mOrb(k) + nodePairM(iPair) .and. btOrb(j) .eq. btOrb(k)) then
          res = groupOfOrb(j)
          return
        end if
      end do

    end function

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillTargetGroupWeights(rdm2spatial)
    !---------------------------------------------------------------------------
    ! Weights of the correlation sum of every group of nodeGroups:
    !   weights((iP,k), j) = rdm2spatial(k, k2, j, l2),  iP = (k2, l2).
    !---------------------------------------------------------------------------
    complex(R64), intent(in), contiguous :: rdm2spatial(:, :, :, :)

    integer(I32) :: iG, iC, iJ, iP

    do iG = 1, size(nodeGroups)
      associate (g => nodeGroups(iG))
        do iJ = 1, size(g % orbs)
          do iC = 1, size(g % colPair)
            iP = g % colPair(iC)
            g % weights(iC, iJ) = rdm2spatial(g % colOrb(iC), nodePairK2(iP), g % orbs(iJ), nodePairL2(iP))
          end do
        end do
      end associate
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillOrbNodes(orbNodes, orbs)
    !---------------------------------------------------------------------------
    ! Orbitals on the angular nodes (channels of their declared m), flattened
    ! to (node point, orbital).
    !---------------------------------------------------------------------------
    use M_Grid_Ylm
    use M_Orbs

    complex(R64), intent(out), allocatable, target :: orbNodes(:, :)
    complex(R64), intent(in), contiguous           :: orbs(:, :)

    complex(R64), contiguous, pointer :: nodes(:, :)
    integer(I32) :: nRad, nNodes, i1

    nRad = Grid_Ylm_nRadial
    nNodes = Grid_Ylm_nNodes
    allocate (orbNodes(nRad * nNodes, Orbs_nOrbsInState))

    !$omp parallel do default(shared) private(i1, nodes)
    do i1 = 1, Orbs_nOrbsInState
      nodes(1:nRad, 1:nNodes) => orbNodes(:, i1)
      call Grid_Ylm_ToNodes(nodes, orbs(:, i1), Grid_Ylm_lmax, Method_Mb_OrbBased_mOfOrb(i1))
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillPairPotentialsYlmNodes(potNodes, orbNodes, time)
    !---------------------------------------------------------------------------
    ! Pair potentials V_(l2,k2) of the sources conjg(phi_l2) phi_k2 of the
    ! unordered pairs on the nodes (the swapped pair has
    ! V_(k2,l2) = conjg(V_(l2,k2))), in parallel over the pairs: on the
    ! channel matrix of the single m = m_k2 - m_l2 of the pair, with radial
    ! metric weights, by SysInteraction_Ylm_FillInteractionPotentialM.
    !---------------------------------------------------------------------------
    use M_Grid_Ylm
    use M_SysInteraction_Ylm

    complex(R64), intent(out), allocatable, target :: potNodes(:, :)
    complex(R64), intent(in), contiguous, target   :: orbNodes(:, :)
    real(R64), intent(in)                          :: time

    complex(R64), allocatable, target :: srcNodes(:)
    complex(R64), allocatable :: srcM(:, :), potM(:, :)
    complex(R64), contiguous, pointer :: nodes(:, :)
    integer(I32) :: nRad, nNodes, nL, iP, iL, mP

    nRad = Grid_Ylm_nRadial
    nNodes = Grid_Ylm_nNodes
    allocate (potNodes(nRad * nNodes, nodeNPots))

    !$omp parallel default(shared) private(iP, iL, mP, nL, srcNodes, srcM, potM, nodes)
    allocate (srcNodes(nRad * nNodes))
    !$omp do schedule(dynamic)
    do iP = 1, size(nodePairPot)
      if (nodePairPot(iP) < 0) cycle
      mP = nodePairM(iP)
      nL = SysInteraction_Ylm_lmax - abs(mP) + 1

      srcNodes(:) = conjg(orbNodes(:, nodePairL2(iP))) * orbNodes(:, nodePairK2(iP))
      nodes(1:nRad, 1:nNodes) => srcNodes
      allocate (srcM(nRad, nL), potM(nRad, nL))
      call Grid_Ylm_NodesToChannels(srcM, nodes, mP)
      do iL = 1, nL
        srcM(:, iL) = srcM(:, iL) * Grid_Ylm_radialMetricWeights(:)
      end do

      call SysInteraction_Ylm_FillInteractionPotentialM(potM, srcM, mP, time, bt1_=-1, bt2_=-1)

      nodes(1:nRad, 1:nNodes) => potNodes(:, nodePairPot(iP))
      call Grid_Ylm_ChannelsToNodes(nodes, potM, mP)
      deallocate (srcM, potM)
    end do
    !$omp end do
    deallocate (srcNodes)
    !$omp end parallel

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine AccumulateTargetSumsYlmNodes(orbNodes, potNodes, h2Cols_, dOrbNodes_)
    !---------------------------------------------------------------------------
    ! Pointwise sums of the node path over blocks of nodes, per target group:
    ! the products X = V_iP * orbs(:,k) of the columns (iP, k) of the group (the
    ! only elementwise work) enter two ZGEMMs,
    !   dOrbNodes_(:,j) += X weights            (correlation sum over (iP, k)),
    !   h2Cols_(i, colOffset + (iP,k)) = sum w conjg(orbs(:,i)) X(:,(iP,k)),
    ! i and j in the group, which touch only the terms allowed by the selection
    ! rule. Only the potentials of unordered pairs are stored; the swapped pair
    ! uses the conjugate.
    !---------------------------------------------------------------------------
    use M_Utils_BlasLib, only: ZGEMM
    use M_Orbs

    complex(R64), intent(in), contiguous              :: orbNodes(:, :)
    complex(R64), intent(in), contiguous              :: potNodes(:, :)
    complex(R64), intent(out), allocatable, optional  :: h2Cols_(:, :)
    complex(R64), intent(inout), contiguous, optional :: dOrbNodes_(:, :)

    !> Nodes per block of the pointwise sums (block work arrays stay in cache).
    integer(I32), parameter :: blockSize = 128

    complex(R64), allocatable :: x(:, :), orbW(:, :), dOrbBlock(:, :), h2Private(:, :)
    integer(I32) :: nOS, nPts, nCols, iB, nB, iG, iC, iJ, nC, nJ, iP, k1, j1
    logical :: h2Q, dOrbQ

    nOS = Orbs_nOrbsInState
    nPts = size(orbNodes, 1)
    nCols = sum([(size(nodeGroups(iG) % colPair), iG=1, size(nodeGroups))])
    h2Q = present(h2Cols_)
    dOrbQ = present(dOrbNodes_)

    if (h2Q) then
      allocate (h2Cols_(nOS, nCols))
      h2Cols_(:, :) = 0.0_R64
    end if

    !$omp parallel default(shared) private(iB, nB, iG, iC, iJ, nC, nJ, iP, k1, j1, x, orbW, dOrbBlock, h2Private)
    allocate (x(blockSize, maxval([(size(nodeGroups(iG) % colPair), iG=1, size(nodeGroups))])))
    allocate (orbW(blockSize, nOS), dOrbBlock(blockSize, nOS))
    if (h2Q) then
      allocate (h2Private(nOS, nCols))
      h2Private(:, :) = 0.0_R64
    end if

    !$omp do schedule(static)
    do iB = 1, nPts, blockSize
      nB = min(blockSize, nPts - iB + 1)

      do iG = 1, size(nodeGroups)
        nC = size(nodeGroups(iG) % colPair)
        nJ = size(nodeGroups(iG) % orbs)
        if (nC .eq. 0) cycle

        do iC = 1, nC
          iP = nodeGroups(iG) % colPair(iC)
          k1 = nodeGroups(iG) % colOrb(iC)
          if (nodePairPot(iP) > 0) then
            x(1:nB, iC) = potNodes(iB:iB + nB - 1, nodePairPot(iP)) * orbNodes(iB:iB + nB - 1, k1)
          else
            x(1:nB, iC) = conjg(potNodes(iB:iB + nB - 1, -nodePairPot(iP))) * orbNodes(iB:iB + nB - 1, k1)
          end if
        end do

        if (dOrbQ) then
          call ZGEMM('N', 'N', nB, nJ, nC, (1.0_R64, 0.0_R64), x, blockSize, nodeGroups(iG) % weights, nC, &
                     (0.0_R64, 0.0_R64), dOrbBlock, blockSize)
          do iJ = 1, nJ
            j1 = nodeGroups(iG) % orbs(iJ)
            dOrbNodes_(iB:iB + nB - 1, j1) = dOrbNodes_(iB:iB + nB - 1, j1) + dOrbBlock(1:nB, iJ)
          end do
        end if

        if (h2Q) then
          do iJ = 1, nJ
            orbW(1:nB, iJ) = orbNodes(iB:iB + nB - 1, nodeGroups(iG) % orbs(iJ)) * nodeWeights(iB:iB + nB - 1)
          end do
          call ZGEMM('C', 'N', nJ, nC, nB, (1.0_R64, 0.0_R64), orbW, blockSize, x, blockSize, &
                     (1.0_R64, 0.0_R64), h2Private(1, nodeGroups(iG) % colOffset + 1), nOS)
        end if
      end do
    end do
    !$omp end do

    if (h2Q) then
      !$omp critical
      h2Cols_(:, :) = h2Cols_(:, :) + h2Private(:, :)
      !$omp end critical
    end if
    !$omp end parallel

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillH2FromTargetColumns(h2, h2Cols)
    !---------------------------------------------------------------------------
    !   h2(i1,l2,k,k2) = 1/2 <orb_i1| V_(l2,k2) |orb_k>  for the column (iP,k),
    ! zero elsewhere; restricted: duplicated for the spin combinations (see
    ! ApplyCorrelationOp).
    !---------------------------------------------------------------------------
    use M_Orbs

    complex(R64), intent(out), allocatable :: h2(:, :, :, :)
    complex(R64), intent(in), contiguous   :: h2Cols(:, :)

    integer(I32) :: nO, nOS, iG, iC, iJ, iP

    nOS = Orbs_nOrbsInState
    nO = Method_Mb_OrbBased_nOrbsSum

    allocate (h2(nO, nO, nO, nO))
    h2(:, :, :, :) = 0.0_R64

    do iG = 1, size(nodeGroups)
      associate (g => nodeGroups(iG))
        do iC = 1, size(g % colPair)
          iP = g % colPair(iC)
          do iJ = 1, size(g % orbs)
            h2(g % orbs(iJ), nodePairL2(iP), g % colOrb(iC), nodePairK2(iP)) = 0.5_R64 * h2Cols(iJ, g % colOffset + iC)
          end do
        end do
      end associate
    end do

    if (Orbs_restrictedQ) then
      h2(nOS + 1:, :nOS, nOS + 1:, :nOS) = h2(:nOS, :nOS, :nOS, :nOS)
      h2(:nOS, nOS + 1:, :nOS, nOS + 1:) = h2(:nOS, :nOS, :nOS, :nOS)
      h2(nOS + 1:, nOS + 1:, nOS + 1:, nOS + 1:) = h2(:nOS, :nOS, :nOS, :nOS)
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillCorrelationRdm2(rdm2spatial, invRdm1, rdm2)
    !---------------------------------------------------------------------------
    ! Spatial 2-RDM contracted with the inverse 1-RDM M = invRdm1, the weights
    ! of the correlation operator:
    !   rdm2spatial(i1,i2,j1,j2) = sum_k rho2(i1,i2,k,j2) M(k,j1),
    ! with the spin-summed 2-RDM in the restricted case.
    !---------------------------------------------------------------------------
    use M_Orbs

    complex(R64), intent(out), allocatable :: rdm2spatial(:, :, :, :)
    complex(R64), intent(in), contiguous   :: invRdm1(:, :)
    complex(R64), intent(in), contiguous   :: rdm2(:, :, :, :)

    complex(R64), allocatable :: correlOp(:, :, :, :)

    integer(I32) :: nO, nOS
    integer(I32) :: j1, k1

    nOS = Orbs_nOrbsInState
    nO = Method_Mb_OrbBased_nOrbsSum

    !------------------------------------
    ! Build spatial 2-RDM (spin-summed for restricted case)
    !------------------------------------

    allocate (rdm2spatial(nOS, nOS, nOS, nOS))

    if (nO .eq. nOS) then
      ! Unrestricted or single-species: direct copy
      rdm2spatial(:, :, :, :) = rdm2(1:nOS, 1:nOS, 1:nOS, 1:nOS)

    else if (nO .eq. 2 * nOS) then
      ! Restricted: sum over spin blocks (αα + αβ + βα + ββ)
      rdm2spatial(:, :, :, :) = rdm2(1:nOS, 1:nOS, 1:nOS, 1:nOS) + &
                                rdm2(nOS + 1:, :nOS, nOS + 1:, :nOS) + &
                                rdm2(:nOS, nOS + 1:, :nOS, nOS + 1:) + &
                                rdm2(nOS + 1:, nOS + 1:, nOS + 1:, nOS + 1:)
    else
      error stop "invalid nO to nOS relation in FillCorrelationRdm2"
    end if

    allocate (correlOp, mold=rdm2spatial)

    ! Apply: rdm2spatial_{i₁,i₂,j₁,l₂} = Σ_k correlOp_{i₁,i₂,k,l₂} · M_{k,j₁}
    correlOp(:, :, :, :) = rdm2spatial(:, :, :, :)
    rdm2spatial(:, :, :, :) = 0.0_R64
    do j1 = 1, nOS
      do k1 = 1, nOS
        rdm2spatial(:, :, j1, :) = rdm2spatial(:, :, j1, :) + correlOp(:, :, k1, :) * invRdm1(k1, j1)
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine AccumulateCorrelation(dOrbs, orbTmps, rdm2spatial, k2, l2)
    !---------------------------------------------------------------------------
    ! Accumulates the correlation contribution of one (k2,l2) source pair:
    !   dOrbs(:,j₁) += Σ_k₁ orbTmps(:,k₁) · ρ²_{k₁,k2,j₁,l2}
    ! where orbTmps holds the interaction applied for source pair (l2,k2).
    !---------------------------------------------------------------------------
    use M_Orbs

    complex(R64), intent(inout), contiguous :: dOrbs(:, :)
    complex(R64), intent(in), contiguous    :: orbTmps(:, :)
    complex(R64), intent(in), contiguous    :: rdm2spatial(:, :, :, :)
    integer(I32), intent(in)                :: k2
    integer(I32), intent(in)                :: l2

    integer(I32) :: j1, k1, j1bt, k1bt, nOS

    nOS = Orbs_nOrbsInState

    do j1 = 1, nOS
      j1bt = Method_Mb_OrbBased_bodyTypeOfOrb(j1)

      do k1 = 1, nOS
        k1bt = Method_Mb_OrbBased_bodyTypeOfOrb(k1)

        if (j1bt .ne. k1bt) cycle

        dOrbs(:, j1) = dOrbs(:, j1) + orbTmps(:, k1) * rdm2spatial(k1, k2, j1, l2)
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillH2CorrelationBlock(h2, orbs, orbTmps, k2, l2)
    !---------------------------------------------------------------------------
    ! Fills the h2 matrix elements of one (k2,l2) source pair:
    !   h²_{i₁,l2,j₁,k2} = ½ ⟨φ_i₁| V_(l2,k2) |φ_j₁⟩
    !---------------------------------------------------------------------------
    use M_Grid
    use M_Orbs

    complex(R64), intent(inout), contiguous :: h2(:, :, :, :)
    complex(R64), intent(in), contiguous    :: orbs(:, :)
    complex(R64), intent(in), contiguous    :: orbTmps(:, :)
    integer(I32), intent(in)                :: k2
    integer(I32), intent(in)                :: l2

    integer(I32) :: i1, j1, i1bt, j1bt, nOS

    nOS = Orbs_nOrbsInState

    !$omp parallel do default(shared) private(i1, j1, i1bt, j1bt) collapse(2)
    do j1 = 1, nOS
      do i1 = 1, nOS
        i1bt = Method_Mb_OrbBased_bodyTypeOfOrb(i1)
        j1bt = Method_Mb_OrbBased_bodyTypeOfOrb(j1)

        if (i1bt .ne. j1bt) cycle

        h2(i1, l2, j1, k2) = 0.5_R64 * Grid_InnerProduct(orbs(:, i1), orbTmps(:, j1))
      end do
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ApplyHartreeFockOp(dOrbs, orbs, rdm1, time)
    !---------------------------------------------------------------------------
    ! Applies the mean-field (Hartree-Fock) operator to orbitals.
    !
    ! Computes both direct (Hartree) and exchange (Fock) contributions:
    !   dOrbs_j += Σ_{k,l} ρ¹_{kl} [W_{jklm} δ_{jm} - W_{kjlm} δ_{km}] · orbs_m
    !
    ! The first term is the direct (Coulomb/Hartree) contribution,
    ! the second is the exchange contribution (same-spin only).
    !---------------------------------------------------------------------------
    use M_Method
    use M_Grid
    use M_SysInteraction
    use M_Orbs

    complex(R64), intent(inout), contiguous, target :: dOrbs(:, :)
    complex(R64), intent(in), contiguous, target    :: orbs(:, :)
    complex(R64), intent(in), contiguous             :: rdm1(:, :)
    real(R64), intent(in)                  :: time

    complex(R64), allocatable :: orbTmps(:, :), orbTmpsSwap(:, :)
    integer(I32) :: nOS
    integer(I32) :: k2, l2
    integer(I32) :: l2bt, k2bt

    allocate (orbTmps, mold=orbs)
    allocate (orbTmpsSwap, mold=orbs)

    nOS = Orbs_nOrbsInState

    ! Loop over unordered pairs (k2 >= l2); the interaction applied for the
    ! swapped source pair reuses the potential via the conjugation symmetry
    do l2 = 1, nOS
      do k2 = l2, nOS
        l2bt = Method_Mb_OrbBased_bodyTypeOfOrb(l2)
        k2bt = Method_Mb_OrbBased_bodyTypeOfOrb(k2)

        ! Skip cross-body-type pairs (vanishing density matrix elements)
        if (l2bt .ne. k2bt) cycle

        orbTmps(:, :) = 0.0_R64

        if (k2 .eq. l2) then
          ! Apply interaction for source pair (l2,k2)
          call Method_Mb_OrbBased_ApplyInteractionOp(orbTmps, orbs, l2, k2, time)

          call AccumulateHartreeFock(dOrbs, orbTmps, rdm1, k2, l2)
        else
          orbTmpsSwap(:, :) = 0.0_R64
          call ApplyInteractionOpPair(orbTmps, orbTmpsSwap, orbs, l2, k2, time)

          call AccumulateHartreeFock(dOrbs, orbTmps, rdm1, k2, l2)
          call AccumulateHartreeFock(dOrbs, orbTmpsSwap, rdm1, l2, k2)
        end if
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine AccumulateHartreeFock(dOrbs, orbTmps, rdm1, k2, l2)
    !---------------------------------------------------------------------------
    ! Accumulates the mean-field contribution of one (k2,l2) source pair,
    ! where orbTmps holds the interaction applied for source pair (l2,k2):
    !   direct:   dOrbs(:,j₁) += orbTmps(:,j₁) · ρ¹_{k2,l2}
    !   exchange: dOrbs(:,k2) -= orbTmps(:,j₁) · ρ¹_{j₁,l2}  (same body type)
    !---------------------------------------------------------------------------
    use M_Orbs

    complex(R64), intent(inout), contiguous :: dOrbs(:, :)
    complex(R64), intent(in), contiguous    :: orbTmps(:, :)
    complex(R64), intent(in), contiguous    :: rdm1(:, :)
    integer(I32), intent(in)                :: k2
    integer(I32), intent(in)                :: l2

    integer(I32) :: j1, j1bt, k2bt, nOS
    complex(R64) :: rdm1Direct

    nOS = Orbs_nOrbsInState
    k2bt = Method_Mb_OrbBased_bodyTypeOfOrb(k2)

    ! Direct (Hartree) density weight; restricted: the stored source pair
    ! (k2,l2) carries both spins, so sum the spin-up and spin-down blocks
    rdm1Direct = rdm1(k2, l2)
    if (Orbs_restrictedQ) rdm1Direct = rdm1Direct + rdm1(nOS + k2, nOS + l2)

    do j1 = 1, nOS
      j1bt = Method_Mb_OrbBased_bodyTypeOfOrb(j1)

      ! Direct (Hartree) term: all j1 get contribution from ρ_{k2,l2}
      dOrbs(:, j1) = dOrbs(:, j1) + orbTmps(:, j1) * rdm1Direct

      ! Exchange (Fock) term: only same-body-type (same-spin) pairs
      if (k2bt .ne. j1bt) cycle
      dOrbs(:, k2) = dOrbs(:, k2) - orbTmps(:, j1) * rdm1(j1, l2)

    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ApplySingleBodyOpLin(dOrbs, orbs, time)
    !---------------------------------------------------------------------------
    ! Linear split-step operator: the full one-body Hamiltonian
    ! (Method_Mb_OrbBased_ApplyLinOp without exact field-free propagation).
    !---------------------------------------------------------------------------
    complex(R64), intent(inout), contiguous :: dOrbs(:, :)
    complex(R64), intent(in), contiguous    :: orbs(:, :)
    real(R64), intent(in)                   :: time

    call Method_Mb_OrbBased_ApplySingleBodyOp(dOrbs, orbs, time)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ApplyGaugeOpLin(dOrbs, orbs, time)
    !---------------------------------------------------------------------------
    ! Linear split-step operator with exact field-free propagation: only the
    ! gauge coupling (T + V_ext are propagated by Method_PropagateExactPart).
    !---------------------------------------------------------------------------
    complex(R64), intent(inout), contiguous :: dOrbs(:, :)
    complex(R64), intent(in), contiguous    :: orbs(:, :)
    real(R64), intent(in)                   :: time

    call Method_Mb_OrbBased_ApplyGaugeOp(dOrbs, orbs, time)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ApplySingleBodyOp(dOrbs, orbs, time, h1_)
    !---------------------------------------------------------------------------
    ! Applies the one-body Hamiltonian (kinetic + external potential + gauge
    ! coupling) to orbitals:
    !   dOrbs(:,i) += (T̂ + V̂_ext + Ĥ_A) · orbs(:,i)
    !
    ! Optionally returns the one-body Hamiltonian matrix h1_ if requested.
    ! For restricted calculations, duplicates the computed block.
    !---------------------------------------------------------------------------
    use M_Method
    use M_Orbs
    use M_Grid
    use M_SysKinetic
    use M_SysPotential

    complex(R64), intent(inout), contiguous, target :: dOrbs(:, :)
    complex(R64), intent(in), contiguous, target    :: orbs(:, :)
    real(R64), intent(in)               :: time
    complex(R64), intent(out), allocatable, optional :: h1_(:, :)

    integer(I32) :: i1, j1, i1bt, j1bt, nOS, nO

    nOS = Orbs_nOrbsInState
    nO = Method_Mb_OrbBased_nOrbsSum

    ! Apply operators (results accumulate in dOrbs)
    call Method_Mb_OrbBased_ApplyKineticOp(dOrbs, orbs, time)
    call Method_Mb_OrbBased_ApplyPotentialOp(dOrbs, orbs, time)
    call Method_Mb_OrbBased_ApplyGaugeOp(dOrbs, orbs, time)

    ! Optionally compute matrix elements from the accumulated result
    if (present(h1_)) then

      if (.not. allocated(h1_)) allocate (h1_(nO, nO))
      h1_(:, :) = 0.0_R64

      !$omp parallel do default(shared) private(i1, j1, i1bt, j1bt) collapse(2)
      do j1 = 1, nOS
        do i1 = 1, nOS

          i1bt = Method_Mb_OrbBased_bodyTypeOfOrb(i1)
          j1bt = Method_Mb_OrbBased_bodyTypeOfOrb(j1)

          if (i1bt .ne. j1bt) cycle

          h1_(i1, j1) = Grid_InnerProduct(orbs(:, i1), dOrbs(:, j1))

        end do
      end do
      !$omp end parallel do

      ! Duplicate for second spin block in restricted case
      if (Orbs_restrictedQ) then
        h1_(nOS + 1:, nOS + 1:) = h1_(:nOS, :nOS)
      end if

    end if

  end subroutine

end submodule
