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

  !> Target group (m, body type) of the definite-m node path: its columns
  !> (pair colPair, orbital colOrb), its orbitals, the offset of its columns
  !> in the h2 columns of all groups and the 2-RDM weights (columns x orbitals).
  type :: T_TargetGroup
    integer(I32), allocatable :: colPair(:)
    integer(I32), allocatable :: colOrb(:)
    integer(I32), allocatable :: orbs(:)
    integer(I32) :: colOffset = 0
    complex(R64), allocatable :: weights(:, :)
  end type

  !> Whether the definite-m node path is bound (Method_Mb_OrbBased_Fabricate).
  logical :: ylmNodesQ = .false.

  !> Skeleton of the node path (SetupYlmNodes): ordered pairs (k2, l2), every
  !> unordered pair followed by its swap, with M = m_k2 - m_l2 and the index of
  !> the stored potential (negative for the swap: conjugate); the target groups
  !> with their total and largest column count; the quadrature weights.
  integer(I32), allocatable :: nodePairK2(:), nodePairL2(:), nodePairM(:), nodePairPot(:)
  integer(I32) :: nodeNPots = 0
  type(T_TargetGroup), allocatable :: nodeGroups(:)
  integer(I32) :: nodeNCols = 0, nodeMaxCols = 0
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
    ylmNodesQ = .false.
    if (Method_Mb_OrbBased_definiteMQ) ylmNodesQ = YlmNodesFastPathQ()
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
    ! then only propagates the gauge coupling (the target bound above); the
    ! methods with a split-step linear part bind Method_PropagateExactPart
    Method_Mb_OrbBased_exactFieldFreeQ = Json_Get("method.mb.orbBased.exactFieldFreeQ", .false.)
    if (Method_Mb_OrbBased_exactFieldFreeQ) then
      Method_Mb_OrbBased_ApplyLinOp => Method_Mb_OrbBased_ApplyGaugeOp
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
    use M_Method

    if (ylmNodesQ) call SetupYlmNodes
    if (associated(Method_PropagateExactPart)) call Method_Mb_OrbBased_SetupFieldFree

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Mb_OrbBased_FillSpatialRdm1(rdm1Spatial, rdm1)
    use M_Orbs

    complex(R64), intent(out), allocatable :: rdm1Spatial(:, :)
    complex(R64), intent(in), contiguous   :: rdm1(:, :)

    integer(I32) :: nOS

    nOS = Orbs_nOrbsInState

    ! Restricted: sum over the spin blocks
    if (Orbs_restrictedQ) then
      rdm1Spatial = rdm1(:nOS, :nOS) + rdm1(nOS + 1:, nOS + 1:)
    else
      rdm1Spatial = rdm1(:nOS, :nOS)
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
    ! Whether the definite-m node quadrature applies: Hermitian Ylm grid and a
    ! body-type-independent Ylm interaction within the node tables.
    !---------------------------------------------------------------------------
    use M_Utils_Json
    use M_Grid_Ylm
    use M_SysInteraction
    use M_SysInteraction_Ylm

    logical :: res

    res = .false.
    if (Grid_Ylm_cProductQ) return
    if (.not. Json_GetExistence("sysInteraction.ylm")) return
    if (.not. SysInteraction_bodyTypeIndependentQ) return
    if (SysInteraction_Ylm_lmax > Grid_Ylm_lmaxNodes) return
    if (SysInteraction_Ylm_mmax > Grid_Ylm_mmaxNodes) return

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
    ! One-body Hamiltonian matrix in the orbital basis (ApplySingleBodyOp):
    !   h¹_{ij} = ⟨φ_i| T̂ + V̂_ext + Ĥ_A |φ_j⟩
    ! where Ĥ_A is the velocity-gauge coupling (zero unless sysGauge is set).
    !---------------------------------------------------------------------------
    complex(R64), intent(out), allocatable :: h1(:, :)
    complex(R64), intent(in), contiguous   :: orbs(:, :)
    real(R64), intent(in)                  :: time

    complex(R64), allocatable :: orbTmps(:, :)

    allocate (orbTmps, mold=orbs)
    orbTmps(:, :) = 0.0_R64
    call Method_Mb_OrbBased_ApplySingleBodyOp(orbTmps, orbs, time, h1_=h1)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillH2(h2, orbs, time)
    !---------------------------------------------------------------------------
    ! Builds the two-body Hamiltonian tensor in the orbital basis:
    !   h²_{i₁i₂j₁j₂} = ½ ∫∫ φ*_{i₁}(r₁) φ*_{i₂}(r₂) W(r₁,r₂) φ_{j₁}(r₁) φ_{j₂}(r₂) dr₁dr₂
    !
    ! The factor ½ accounts for double-counting in the Hamiltonian. Elements
    ! with orbitals of different body type within a pair are zero.
    !---------------------------------------------------------------------------
    complex(R64), intent(out), allocatable :: h2(:, :, :, :)
    complex(R64), intent(in), contiguous   :: orbs(:, :)
    real(R64), intent(in)                  :: time

    integer(I32) :: nO

    nO = Method_Mb_OrbBased_nOrbsSum
    allocate (h2(nO, nO, nO, nO))
    h2(:, :, :, :) = 0.0_R64

    call AccumulatePairTerms(orbs, time, h2_=h2)
    call DuplicateRestrictedSpinBlocks(h2)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillH2Block(h2, orbs, orbTmps, i2, j2)
    !---------------------------------------------------------------------------
    ! First-particle elements of one second-particle pair (i2,j2), orbTmps =
    ! V_(i2,j2) orbs:  h²_{i₁,i2,j₁,j2} = ½ ⟨φ_i₁| V_(i2,j2) |φ_j₁⟩
    !---------------------------------------------------------------------------
    use M_Orbs

    complex(R64), intent(inout), contiguous :: h2(:, :, :, :)
    complex(R64), intent(in), contiguous    :: orbs(:, :)
    complex(R64), intent(in), contiguous    :: orbTmps(:, :)
    integer(I32), intent(in)                :: i2
    integer(I32), intent(in)                :: j2

    complex(R64), allocatable :: mat(:, :)
    integer(I32) :: nOS

    nOS = Orbs_nOrbsInState
    call FillOrbitalMatrix(mat, orbs, orbTmps)
    h2(:nOS, i2, :nOS, j2) = 0.5_R64 * mat(:, :)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillOrbitalMatrix(mat, orbs, orbTmps)
    !---------------------------------------------------------------------------
    ! mat(i,j) = ⟨φ_i|orbTmps(:,j)⟩ for orbitals of equal body type, zero
    ! otherwise (nOrbsInState x nOrbsInState).
    !---------------------------------------------------------------------------
    use M_Grid
    use M_Orbs

    complex(R64), intent(out), allocatable :: mat(:, :)
    complex(R64), intent(in), contiguous   :: orbs(:, :)
    complex(R64), intent(in), contiguous   :: orbTmps(:, :)

    integer(I32) :: i1, j1, nOS

    nOS = Orbs_nOrbsInState
    allocate (mat(nOS, nOS))
    mat(:, :) = 0.0_R64

    !$omp parallel do default(shared) private(i1, j1) collapse(2)
    do j1 = 1, nOS
      do i1 = 1, nOS
        if (Method_Mb_OrbBased_bodyTypeOfOrb(i1) .ne. Method_Mb_OrbBased_bodyTypeOfOrb(j1)) cycle
        mat(i1, j1) = Grid_InnerProduct(orbs(:, i1), orbTmps(:, j1))
      end do
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine DuplicateRestrictedSpinBlocks(h2)
    !---------------------------------------------------------------------------
    ! Restricted orbitals: copies the spatial block of h2 to the spin
    ! combinations αβ, βα and ββ (they share the spatial integrals; the Coeffs
    ! backend consumes h2 in the full spin-orbital basis).
    !---------------------------------------------------------------------------
    use M_Orbs

    complex(R64), intent(inout), contiguous :: h2(:, :, :, :)

    integer(I32) :: nOS

    if (.not. Orbs_restrictedQ) return

    nOS = Orbs_nOrbsInState
    h2(nOS + 1:, :nOS, nOS + 1:, :nOS) = h2(:nOS, :nOS, :nOS, :nOS)
    h2(:nOS, nOS + 1:, :nOS, nOS + 1:) = h2(:nOS, :nOS, :nOS, :nOS)
    h2(nOS + 1:, nOS + 1:, nOS + 1:, nOS + 1:) = h2(:nOS, :nOS, :nOS, :nOS)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ApplyKineticOp(dOrbs, orbs, time)
    !---------------------------------------------------------------------------
    ! dOrbs(:,i) += T̂ · orbs(:,i), with the mass of the body type of orbital i.
    !---------------------------------------------------------------------------
    use M_SysKinetic

    complex(R64), intent(inout), contiguous :: dOrbs(:, :)
    complex(R64), intent(in), contiguous    :: orbs(:, :)
    real(R64), intent(in)                   :: time

    call AccumulatePerOrbitalOp(dOrbs, orbs, time, SysKinetic_MultiplyWithKineticOp)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ApplyGaugeOp(dOrbs, orbs, time)
    !---------------------------------------------------------------------------
    ! dOrbs(:,i) += (A(t)·p̂/m + A²(t)/(2m)) · orbs(:,i), with the mass of the
    ! body type of orbital i.
    !---------------------------------------------------------------------------
    use M_SysGauge

    complex(R64), intent(inout), contiguous :: dOrbs(:, :)
    complex(R64), intent(in), contiguous    :: orbs(:, :)
    real(R64), intent(in)                   :: time

    call AccumulatePerOrbitalOp(dOrbs, orbs, time, SysGauge_MultiplyWithGaugeOp)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine AccumulatePerOrbitalOp(dOrbs, orbs, time, MultiplyWithOp)
    !---------------------------------------------------------------------------
    ! dOrbs(:,i) += O · orbs(:,i) for every orbital, with O applied for the
    ! body type of orbital i (kinetic and gauge operators share the interface).
    !---------------------------------------------------------------------------
    use M_Grid
    use M_Orbs
    use M_SysKinetic

    complex(R64), intent(inout), contiguous :: dOrbs(:, :)
    complex(R64), intent(in), contiguous    :: orbs(:, :)
    real(R64), intent(in)                   :: time
    procedure(I_SysKinetic_MultiplyWithKineticOp) :: MultiplyWithOp

    complex(R64), allocatable :: orbTmp(:)
    integer(I32) :: i1

    !$omp parallel default(shared) private(i1, orbTmp)
    allocate (orbTmp(Grid_nPoints))
    !$omp do
    do i1 = 1, Orbs_nOrbsInState
      call MultiplyWithOp(orbTmp, orbs(:, i1), time, bt_=Method_Mb_OrbBased_bodyTypeOfOrb(i1))
      dOrbs(:, i1) = dOrbs(:, i1) + orbTmp(:)
    end do
    !$omp end do
    deallocate (orbTmp)
    !$omp end parallel

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ApplyPotentialOp(dOrbs, orbs, time)
    !---------------------------------------------------------------------------
    ! dOrbs(:,i) += V̂_ext · orbs(:,i); the potential is computed once if
    ! body-type-independent, otherwise per body type.
    !---------------------------------------------------------------------------
    use M_Orbs
    use M_Method_Mb
    use M_SysPotential

    complex(R64), intent(inout), contiguous :: dOrbs(:, :)
    complex(R64), intent(in), contiguous    :: orbs(:, :)
    real(R64), intent(in)                   :: time

    integer(I32) :: nBt, bt
    complex(R64), allocatable :: externalPotential(:)

    if (SysPotential_bodyTypeIndependentQ) then
      call SysPotential_FillExternalPotential(externalPotential, time)
      call MultiplyPotentialOnOrbRange(dOrbs, orbs, externalPotential, 1, Orbs_nOrbsInState, &
                                       SysPotential_MultiplyWithExternalPotential)
    else
      ! Restricted: only the shared spatial set is stored; the body-type-1
      ! potential covers all stored orbitals (spin up and down coincide)
      nBt = Method_Mb_nBodyTypes
      if (Orbs_restrictedQ) nBt = 1
      do bt = 1, nBt
        call SysPotential_FillExternalPotential(externalPotential, time, bt_=bt)
        call MultiplyPotentialOnOrbRange(dOrbs, orbs, externalPotential, Method_Mb_OrbBased_nOrbsStart(bt), &
                                         Method_Mb_OrbBased_nOrbsEnd(bt), SysPotential_MultiplyWithExternalPotential)
      end do
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ApplyInteractionOp(dOrbs, orbs, i2, j2, time)
    !---------------------------------------------------------------------------
    ! Two-body interaction operator of the source pair (i2,j2):
    !   dOrbs(:,i1) += ∫ W(r₁,r₂) φ*_{i2}(r₂) φ_{j2}(r₂) dr₂ · orbs(:,i1)
    !---------------------------------------------------------------------------
    complex(R64), intent(inout), contiguous :: dOrbs(:, :)
    complex(R64), intent(in), contiguous    :: orbs(:, :)
    integer(I32), intent(in)                :: i2
    integer(I32), intent(in)                :: j2
    real(R64), intent(in)                   :: time

    call ApplyInteraction(dOrbs, orbs, i2, j2, time)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ApplyInteraction(dOrbs, orbs, i2, j2, time, dOrbsSwap_)
    !---------------------------------------------------------------------------
    ! dOrbs += V_(i2,j2) orbs (ApplyInteractionOp) and, if present,
    ! dOrbsSwap_ += V_(j2,i2) orbs from the conjugate potential
    ! V_(j2,i2) = conjg(V_(i2,j2)) (requires SysInteraction_conjSymmetricQ),
    ! which halves the potential solves over unordered pairs. The potential is
    ! computed once if body-type-independent, otherwise per target body type.
    !---------------------------------------------------------------------------
    use M_Method_Mb
    use M_SysInteraction
    use M_Orbs

    complex(R64), intent(inout), contiguous           :: dOrbs(:, :)
    complex(R64), intent(in), contiguous              :: orbs(:, :)
    integer(I32), intent(in)                          :: i2
    integer(I32), intent(in)                          :: j2
    real(R64), intent(in)                             :: time
    complex(R64), intent(inout), contiguous, optional :: dOrbsSwap_(:, :)

    integer(I32) :: nBt, bt
    complex(R64), allocatable :: src(:)

    ! Source density ρ(r₂) = φ*_{i2}(r₂) φ_{j2}(r₂)
    call SysInteraction_FillInteractionSrc(src, orbs(:, i2), orbs(:, j2))

    if (SysInteraction_bodyTypeIndependentQ) then
      call ApplyPotential(-1, -1, 1, Orbs_nOrbsInState)
    else
      ! Restricted: only the shared spatial set is stored; the body-type-1
      ! interaction covers all stored orbitals (spin up and down coincide)
      nBt = Method_Mb_nBodyTypes
      if (Orbs_restrictedQ) nBt = 1
      do bt = 1, nBt
        call ApplyPotential(bt, Method_Mb_OrbBased_bodyTypeOfOrb(i2), &
                            Method_Mb_OrbBased_nOrbsStart(bt), Method_Mb_OrbBased_nOrbsEnd(bt))
      end do
    end if

  contains

    !> Potential for the target body type btTarget, applied to the orbitals iStart..iEnd.
    subroutine ApplyPotential(btTarget, btSource, iStart, iEnd)
      integer(I32), intent(in) :: btTarget
      integer(I32), intent(in) :: btSource
      integer(I32), intent(in) :: iStart
      integer(I32), intent(in) :: iEnd

      complex(R64), allocatable :: pot(:), potConjg(:)

      call SysInteraction_FillInteractionPotential(pot, src, time, btTarget, btSource)
      call MultiplyPotentialOnOrbRange(dOrbs, orbs, pot, iStart, iEnd, SysInteraction_MultiplyWithInteractionPotential)
      if (.not. present(dOrbsSwap_)) return

      allocate (potConjg, mold=pot)
      call SysInteraction_ConjugateInteractionPotential(potConjg, pot)
      call MultiplyPotentialOnOrbRange(dOrbsSwap_, orbs, potConjg, iStart, iEnd, &
                                       SysInteraction_MultiplyWithInteractionPotential)

    end subroutine

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine MultiplyPotentialOnOrbRange(dOrbs, orbs, pot, iOrbStart, iOrbEnd, MultiplyWithPotential)
    !---------------------------------------------------------------------------
    ! Accumulates dOrbs(:,i) += V · orbs(:,i) for i = iOrbStart, ..., iOrbEnd
    ! (external or interaction potential, the same interface).
    !---------------------------------------------------------------------------
    use M_Grid
    use M_SysPotential

    complex(R64), intent(inout), contiguous :: dOrbs(:, :)
    complex(R64), intent(in), contiguous    :: orbs(:, :)
    complex(R64), intent(in), contiguous    :: pot(:)
    integer(I32), intent(in)                :: iOrbStart
    integer(I32), intent(in)                :: iOrbEnd
    procedure(I_SysPotential_MultiplyWithExternalPotential) :: MultiplyWithPotential

    integer(I32) :: iOrb
    complex(R64), allocatable :: orbTmp(:)

    !$omp parallel default(shared) private(iOrb, orbTmp)
    allocate (orbTmp(Grid_nPoints))
    !$omp do
    do iOrb = iOrbStart, iOrbEnd
      call MultiplyWithPotential(orbTmp, pot, orbs(:, iOrb))
      dOrbs(:, iOrb) = dOrbs(:, iOrb) + orbTmp(:)
    end do
    !$omp end do
    deallocate (orbTmp)
    !$omp end parallel

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ApplyCorrelationOp(dOrbs, orbs, invRdm1, rdm2, time, h2_)
    !---------------------------------------------------------------------------
    ! Beyond-mean-field correlation operator of the orbital equations:
    !   dOrbs(:,j₁) += Σ_{k₁,k2,l2} V_(l2,k2) orbs(:,k₁) rdm2spatial(k₁,k2,j₁,l2)
    ! with the spatial 2-RDM contracted with the inverse 1-RDM
    ! (FillCorrelationRdm2). Optionally returns the two-body Hamiltonian h2_
    ! from the same pair potentials.
    !---------------------------------------------------------------------------
    complex(R64), intent(inout), contiguous, target :: dOrbs(:, :)
    complex(R64), intent(in), contiguous, target    :: orbs(:, :)
    complex(R64), intent(in), contiguous             :: invRdm1(:, :)
    complex(R64), intent(in), contiguous             :: rdm2(:, :, :, :)
    real(R64), intent(in)                  :: time
    complex(R64), intent(out), allocatable, optional :: h2_(:, :, :, :)

    complex(R64), allocatable :: rdm2spatial(:, :, :, :)
    integer(I32) :: nO

    call FillCorrelationRdm2(rdm2spatial, invRdm1, rdm2)

    if (present(h2_)) then
      nO = Method_Mb_OrbBased_nOrbsSum
      allocate (h2_(nO, nO, nO, nO))
      h2_(:, :, :, :) = 0.0_R64
    end if

    call AccumulatePairTerms(orbs, time, dOrbs_=dOrbs, h2_=h2_, rdm2spatial_=rdm2spatial)

    if (present(h2_)) call DuplicateRestrictedSpinBlocks(h2_)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine AccumulatePairTerms(orbs, time, dOrbs_, h2_, rdm2spatial_, rdm1_)
    !---------------------------------------------------------------------------
    ! Loop over the source pairs (l2,k2) of equal body type shared by h2, the
    ! correlation and the mean-field operator: with orbTmps = V_(l2,k2) orbs,
    ! fills the h2_ block of the pair (FillH2Block) and accumulates into dOrbs_
    ! (required with rdm2spatial_ or rdm1_) the correlation term with the
    ! weights rdm2spatial_ (AccumulateCorrelation) and/or the mean-field term
    ! with rdm1_ (AccumulateHartreeFock). Unordered pairs (k2 >= l2); the
    ! swapped pair reuses the potential via the conjugation symmetry if the
    ! backend has it.
    !---------------------------------------------------------------------------
    use M_SysInteraction
    use M_Orbs

    complex(R64), intent(in), contiguous              :: orbs(:, :)
    real(R64), intent(in)                             :: time
    complex(R64), intent(inout), contiguous, optional :: dOrbs_(:, :)
    complex(R64), intent(inout), contiguous, optional :: h2_(:, :, :, :)
    complex(R64), intent(in), contiguous, optional    :: rdm2spatial_(:, :, :, :)
    complex(R64), intent(in), contiguous, optional    :: rdm1_(:, :)

    complex(R64), allocatable :: orbTmps(:, :), orbTmpsSwap(:, :)
    integer(I32) :: k2, l2

    allocate (orbTmps, mold=orbs)
    allocate (orbTmpsSwap, mold=orbs)

    do l2 = 1, Orbs_nOrbsInState
      do k2 = l2, Orbs_nOrbsInState
        if (Method_Mb_OrbBased_bodyTypeOfOrb(l2) .ne. Method_Mb_OrbBased_bodyTypeOfOrb(k2)) cycle

        orbTmps(:, :) = 0.0_R64

        if (k2 .eq. l2) then
          call Method_Mb_OrbBased_ApplyInteractionOp(orbTmps, orbs, l2, k2, time)
          call AccumulatePair(orbTmps, k2, l2)
        else
          orbTmpsSwap(:, :) = 0.0_R64
          if (SysInteraction_conjSymmetricQ) then
            call ApplyInteraction(orbTmps, orbs, l2, k2, time, dOrbsSwap_=orbTmpsSwap)
          else
            call Method_Mb_OrbBased_ApplyInteractionOp(orbTmps, orbs, l2, k2, time)
            call Method_Mb_OrbBased_ApplyInteractionOp(orbTmpsSwap, orbs, k2, l2, time)
          end if
          call AccumulatePair(orbTmps, k2, l2)
          call AccumulatePair(orbTmpsSwap, l2, k2)
        end if
      end do
    end do

  contains

    !> Terms of the source pair (lP,kP), vOrbs = V_(lP,kP) orbs.
    subroutine AccumulatePair(vOrbs, kP, lP)
      complex(R64), intent(in), contiguous :: vOrbs(:, :)
      integer(I32), intent(in)             :: kP
      integer(I32), intent(in)             :: lP

      if (present(h2_)) call FillH2Block(h2_, orbs, vOrbs, lP, kP)
      if (present(rdm2spatial_)) call AccumulateCorrelation(dOrbs_, vOrbs, rdm2spatial_, kP, lP)
      if (present(rdm1_)) call AccumulateHartreeFock(dOrbs_, vOrbs, rdm1_, kP, lP)

    end subroutine

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ApplyCorrelationOpYlmNodes(dOrbs, orbs, invRdm1, rdm2, time, h2_)
    !---------------------------------------------------------------------------
    ! ApplyCorrelationOp for orbitals of declared m on the Gauss-Legendre
    ! nodes in cos(theta), where all products are pointwise: the orbitals are
    ! transformed once, each pair potential once and dOrbs once at the end.
    ! The pair (l2,k2) carries M = m_k2 - m_l2; only the terms with
    ! m_k1 + M = m_j1 (correlation) and m_j1 + M = m_i1 (h2) are kept, the
    ! others vanish. Equal to ApplyCorrelationOp up to round-off (A_Method.md).
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
      call Grid_Ylm_FromNodes(dOrb, nodes, Method_Mb_OrbBased_mOfOrb(i1))
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
    ! Skeleton of the node path from the declared m of the stored orbitals: the
    ! pair list, the target groups (the blocks Orbs_blockOfOrb) of the
    ! selection rule m_k + M_iP = m_j with their columns (iP, k), and the
    ! quadrature weights on the node points.
    !---------------------------------------------------------------------------
    use M_Grid_Ylm
    use M_SysInteraction_Ylm
    use M_Orbs

    integer(I32), allocatable :: mOrb(:), btOrb(:), targetGroup(:, :), pairOfCol(:, :), orbOfCol(:, :)
    integer(I32) :: nOS, nPairs, iG, iP, j1, k1, k2, l2

    nOS = Orbs_nOrbsInState
    mOrb = Method_Mb_OrbBased_mOfOrb(1:nOS)
    btOrb = Method_Mb_OrbBased_bodyTypeOfOrb(1:nOS)

    ! Ordered pairs (k2, l2) of equal body type within the |m| cap of the
    ! potentials (pairs beyond it are left out, as in the generic operator)
    nPairs = count(spread(btOrb, 1, nOS) .eq. spread(btOrb, 2, nOS) .and. &
                   abs(spread(mOrb, 1, nOS) - spread(mOrb, 2, nOS)) <= SysInteraction_Ylm_mmax)
    allocate (nodePairK2(nPairs), nodePairL2(nPairs), nodePairM(nPairs), nodePairPot(nPairs))

    iP = 0
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

    ! Target group of the column (iP, k): the block of the orbitals with
    ! m_k + M_iP and the body type of k (0 if there is none)
    allocate (targetGroup(nOS, nPairs))
    do iP = 1, nPairs
      do k1 = 1, nOS
        j1 = findloc(mOrb .eq. mOrb(k1) + nodePairM(iP) .and. btOrb .eq. btOrb(k1), .true., dim=1)
        targetGroup(k1, iP) = 0
        if (j1 > 0) targetGroup(k1, iP) = Orbs_blockOfOrb(j1)
      end do
    end do
    orbOfCol = spread([(k1, k1=1, nOS)], 2, nPairs)
    pairOfCol = spread([(iP, iP=1, nPairs)], 1, nOS)

    allocate (nodeGroups(maxval(Orbs_blockOfOrb)))
    do iG = 1, size(nodeGroups)
      associate (g => nodeGroups(iG))
        g % orbs = pack([(k1, k1=1, nOS)], Orbs_blockOfOrb .eq. iG)
        g % colPair = pack(pairOfCol, targetGroup .eq. iG)
        g % colOrb = pack(orbOfCol, targetGroup .eq. iG)
        allocate (g % weights(size(g % colPair), size(g % orbs)))
        g % colOffset = nodeNCols
        nodeNCols = nodeNCols + size(g % colPair)
        nodeMaxCols = max(nodeMaxCols, size(g % colPair))
      end associate
    end do

    ! Quadrature weights of the full measure on the node points
    nodeWeights = reshape(spread(Grid_Ylm_radialWeights(:), 2, Grid_Ylm_nNodes) * &
                          spread(Grid_Ylm_nodeWeights(:), 1, Grid_Ylm_nRadial), [Grid_Ylm_nRadial * Grid_Ylm_nNodes])

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
      call Grid_Ylm_ToNodes(nodes, orbs(:, i1), Method_Mb_OrbBased_mOfOrb(i1))
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillPairPotentialsYlmNodes(potNodes, orbNodes, time)
    !---------------------------------------------------------------------------
    ! Pair potentials V_(l2,k2) of the sources conjg(phi_l2) phi_k2 of the
    ! unordered pairs on the nodes, in parallel over the pairs, on the channel
    ! matrix of the single m = m_k2 - m_l2 of the pair.
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

      call SysInteraction_Ylm_FillInteractionPotentialM(potM, srcM, mP, time)

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
    ! the products X = V_iP * orbs(:,k) of the columns (iP, k) of the group
    ! enter two ZGEMMs (i, j in the group),
    !   dOrbNodes_(:,j) += X weights,
    !   h2Cols_(i, colOffset + (iP,k)) = sum w conjg(orbs(:,i)) X(:,(iP,k)).
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
    integer(I32) :: nOS, nPts, iB, nB, iG, iC, iJ, nC, nJ, iP, k1, j1
    logical :: h2Q, dOrbQ

    nOS = Orbs_nOrbsInState
    nPts = size(orbNodes, 1)
    h2Q = present(h2Cols_)
    dOrbQ = present(dOrbNodes_)

    if (h2Q) then
      allocate (h2Cols_(nOS, nodeNCols))
      h2Cols_(:, :) = 0.0_R64
    end if

    !$omp parallel default(shared) private(iB, nB, iG, iC, iJ, nC, nJ, iP, k1, j1, x, orbW, dOrbBlock, h2Private)
    allocate (x(blockSize, nodeMaxCols))
    allocate (orbW(blockSize, nOS), dOrbBlock(blockSize, nOS))
    if (h2Q) then
      allocate (h2Private(nOS, nodeNCols))
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
    ! zero elsewhere.
    !---------------------------------------------------------------------------
    complex(R64), intent(out), allocatable :: h2(:, :, :, :)
    complex(R64), intent(in), contiguous   :: h2Cols(:, :)

    integer(I32) :: nO, iG, iC, iJ, iP

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

    call DuplicateRestrictedSpinBlocks(h2)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillCorrelationRdm2(rdm2spatial, invRdm1, rdm2)
    !---------------------------------------------------------------------------
    ! Spatial 2-RDM contracted with the inverse 1-RDM M = invRdm1, the weights
    ! of the correlation operator:
    !   rdm2spatial(i1,i2,j1,j2) = sum_k rho2(i1,i2,k,j2) M(k,j1),
    ! with the spin-summed 2-RDM (αα + αβ + βα + ββ) in the restricted case.
    !---------------------------------------------------------------------------
    use M_Orbs

    complex(R64), intent(out), allocatable :: rdm2spatial(:, :, :, :)
    complex(R64), intent(in), contiguous   :: invRdm1(:, :)
    complex(R64), intent(in), contiguous   :: rdm2(:, :, :, :)

    complex(R64), allocatable :: rho2(:, :, :, :)
    integer(I32) :: nOS, j1, k1

    nOS = Orbs_nOrbsInState

    if (Orbs_restrictedQ) then
      rho2 = rdm2(:nOS, :nOS, :nOS, :nOS) + rdm2(nOS + 1:, :nOS, nOS + 1:, :nOS) + &
             rdm2(:nOS, nOS + 1:, :nOS, nOS + 1:) + rdm2(nOS + 1:, nOS + 1:, nOS + 1:, nOS + 1:)
    else
      rho2 = rdm2(:nOS, :nOS, :nOS, :nOS)
    end if

    allocate (rdm2spatial, mold=rho2)
    rdm2spatial(:, :, :, :) = 0.0_R64
    do j1 = 1, nOS
      do k1 = 1, nOS
        rdm2spatial(:, :, j1, :) = rdm2spatial(:, :, j1, :) + rho2(:, :, k1, :) * invRdm1(k1, j1)
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
    complex(R64), intent(inout), contiguous, target :: dOrbs(:, :)
    complex(R64), intent(in), contiguous, target    :: orbs(:, :)
    complex(R64), intent(in), contiguous             :: rdm1(:, :)
    real(R64), intent(in)                  :: time

    call AccumulatePairTerms(orbs, time, dOrbs_=dOrbs, rdm1_=rdm1)

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
  subroutine ApplySingleBodyOp(dOrbs, orbs, time, h1_)
    !---------------------------------------------------------------------------
    ! Applies the one-body Hamiltonian (kinetic + external potential + gauge
    ! coupling) to orbitals:
    !   dOrbs(:,i) += (T̂ + V̂_ext + Ĥ_A) · orbs(:,i)
    !
    ! Optionally returns h1_(i,j) = ⟨φ_i|dOrbs(:,j)⟩ from the accumulated
    ! result (equal body type; duplicated for the second spin block if
    ! restricted).
    !---------------------------------------------------------------------------
    use M_Orbs

    complex(R64), intent(inout), contiguous, target :: dOrbs(:, :)
    complex(R64), intent(in), contiguous, target    :: orbs(:, :)
    real(R64), intent(in)               :: time
    complex(R64), intent(out), allocatable, optional :: h1_(:, :)

    complex(R64), allocatable :: h1Spatial(:, :)
    integer(I32) :: nOS, nO

    call Method_Mb_OrbBased_ApplyKineticOp(dOrbs, orbs, time)
    call Method_Mb_OrbBased_ApplyPotentialOp(dOrbs, orbs, time)
    call Method_Mb_OrbBased_ApplyGaugeOp(dOrbs, orbs, time)

    if (.not. present(h1_)) return

    nOS = Orbs_nOrbsInState
    nO = Method_Mb_OrbBased_nOrbsSum
    call FillOrbitalMatrix(h1Spatial, orbs, dOrbs)

    allocate (h1_(nO, nO))
    h1_(:, :) = 0.0_R64
    h1_(:nOS, :nOS) = h1Spatial(:, :)
    if (Orbs_restrictedQ) h1_(nOS + 1:, nOS + 1:) = h1Spatial(:, :)

  end subroutine

end submodule
