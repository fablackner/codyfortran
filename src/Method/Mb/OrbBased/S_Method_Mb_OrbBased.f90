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
  ! Delegates to specific methods: TDCI, TDHX, MCTDHX, TD-2RDM, TD-2RDM-StaticOrbs, TDKS
  !
  ! The definite-m node path of the correlation operator and of h2 (Ylm grids)
  ! is the child submodule S_Method_Mb_OrbBased_YlmNodes.
  !-----------------------------------------------------------------------------

  implicit none

  !=============================================================================
  ! local data
  !=============================================================================

  !> Whether the definite-m node path is bound (Method_Mb_OrbBased_Fabricate).
  logical :: ylmNodesQ = .false.

  !> Whether the gauge coupling is bound (a sysGauge block exists;
  !> Method_Mb_OrbBased_Fabricate): it is then a term of h1
  !> (Method_Mb_OrbBased_H1TimeIndependentQ).
  logical :: gaugeQ = .false.

  !> Static part of h1 of the static orbitals (Method_Mb_OrbBased_RefillStaticH1:
  !> h1 at its first call), reset by Method_Mb_OrbBased_Setup.
  complex(R64), allocatable :: h0Static(:, :)

  !=============================================================================
  ! procedures of the child submodule S_Method_Mb_OrbBased_YlmNodes
  !=============================================================================

  interface
    !> Whether the definite-m node quadrature applies: Hermitian Ylm grid and
    !> a body-type-independent Ylm interaction within the node tables.
    module function YlmNodesFastPathQ() result(res)
      logical :: res
    end function

    !> Skeleton of the node path from the declared m of the stored orbitals.
    module subroutine SetupYlmNodes
    end subroutine

    !> Method_Mb_OrbBased_ApplyCorrelationOp for orbitals of declared m.
    module subroutine ApplyCorrelationOpYlmNodes(dOrbs, orbs, invRdm1, rdm2, time, h2_)
      complex(R64), intent(inout), contiguous, target :: dOrbs(:, :)
      complex(R64), intent(in), contiguous, target    :: orbs(:, :)
      complex(R64), intent(in), contiguous             :: invRdm1(:, :)
      complex(R64), intent(in), contiguous             :: rdm2(:, :, :, :)
      real(R64), intent(in)                  :: time
      complex(R64), intent(out), allocatable, optional :: h2_(:, :, :, :)
    end subroutine

    !> Method_Mb_OrbBased_FillH2 for orbitals of declared m.
    module subroutine FillH2YlmNodes(h2, orbs, time)
      complex(R64), intent(out), allocatable :: h2(:, :, :, :)
      complex(R64), intent(in), contiguous   :: orbs(:, :)
      real(R64), intent(in)                  :: time
    end subroutine
  end interface

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Mb_OrbBased_Fabricate
    !---------------------------------------------------------------------------
    ! Initializes orbital-based infrastructure and dispatches to concrete method.
    !
    ! JSON keys read:
    !   - method.mb.orbBased.nOrbs:                 Orbitals per body type
    !   - method.mb.orbBased.regularizationParameter: Numerical stability (default 1e-10)
    !   - method.mb.orbBased.mOfOrbs:               Declared m of every spin
    !                                               orbital (Ylm grids, optional)
    !   - method.mb.orbBased.adiabaticSwitching.switchTime: interaction switched
    !                                               on over -T < t < 0
    !                                               (S_Method_Mb_OrbBased_Switching)
    !
    ! With a "fieldFree" block (M_FieldFree, fabricated before Method) the
    ! linear split-step operator carries only the gauge coupling.
    !
    ! Binds default operator-application procedures, then branches to the
    ! selected method (TDCI, TDHX, MCTDHX, TD-2RDM, TD-2RDM-StaticOrbs or TDKS).
    !---------------------------------------------------------------------------
    use M_Utils_Json
    use M_Utils_Say
    use M_Utils_NoOpProcedures
    use M_Grid_Ylm
    use M_FieldFree
    use M_Method
    use M_Method_Mb
    use M_Method_Mb_OrbBased_Tdci
    use M_Method_Mb_OrbBased_Tdhx
    use M_Method_Mb_OrbBased_Mctdhx
    use M_Method_Mb_OrbBased_Td2rdm
    use M_Method_Mb_OrbBased_Td2rdmStaticOrbs
    use M_Method_Mb_OrbBased_Tdks

    !> Endings of the method.mb.orbBased.exactFieldFree* keys that became the
    !> fieldFree block.
    character(len=11), parameter :: FORMER_KEYS(7) = &
                                    [character(len=11) :: "Q", "GpuQ", "Cut", "Tol", "KrylovDim", "MaxRestarts", "MeanFieldQ"]

    integer(I32) :: i, ibt, index, counter

    call Say_Fabricate("method.mb.orbBased")

    !------------------------------------
    ! parse orbital configuration
    !------------------------------------

    Method_Mb_OrbBased_nOrbs = Json_Get("method.mb.orbBased.nOrbs", [(0, i=1, Method_Mb_nBodyTypes)])

    Method_Mb_OrbBased_nOrbsSum = sum(Method_Mb_OrbBased_nOrbs)

    !------------------------------------
    ! build orbital index mappings
    !------------------------------------

    ! A repeated fabricate replaces the index mappings of the earlier one
    if (allocated(Method_Mb_OrbBased_bodyTypeOfOrb)) deallocate (Method_Mb_OrbBased_bodyTypeOfOrb)
    if (allocated(Method_Mb_OrbBased_nOrbsStart)) deallocate (Method_Mb_OrbBased_nOrbsStart)
    if (allocated(Method_Mb_OrbBased_nOrbsEnd)) deallocate (Method_Mb_OrbBased_nOrbsEnd)
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
    Method_Mb_OrbBased_OrthonormalizeKeepingState => OrthonormalizeKeepingStateUnsupported

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

    gaugeQ = Json_GetExistence("sysGauge")
    if (gaugeQ) then
      Method_Mb_OrbBased_ApplyGaugeOp => ApplyGaugeOp
    else
      Method_Mb_OrbBased_ApplyGaugeOp => NoOpProcedures_ApplyPotentialOp
    end if

    if (Json_GetExistence("sysInteraction")) then
      Method_Mb_OrbBased_ApplyInteractionOp => ApplyInteractionOp
    else
      Method_Mb_OrbBased_ApplyInteractionOp => NoOpProcedures_ApplyInteractionOp
    end if

    ! Exact field-free part (M_FieldFree, a "fieldFree" block; a bound
    ! FieldFree_Propagate): the linear split-step integrator then only
    ! propagates the gauge coupling (the target bound above) and the methods
    ! with a split-step linear part bind Method_PropagateExactPart to
    ! FieldFree_Propagate. A fabricate of FieldFree that did not read this
    ! input would leave the exact part of another one (or none) bound.
    if (Json_GetExistence("fieldFree") .neqv. associated(FieldFree_Propagate)) &
      error stop "method.mb.orbBased: call FieldFree_Fabricate before Method_Fabricate (it reads the fieldFree block)"
    if (associated(FieldFree_Propagate)) then
      Method_Mb_OrbBased_ApplyLinOp => ApplyDrivingOp
    else
      Method_Mb_OrbBased_ApplyLinOp => ApplySingleBodyOpLin
    end if

    ! The exact part was configured here before it became the fieldFree block:
    ! an input with the former keys would silently run without it
    do i = 1, size(FORMER_KEYS)
      if (Json_GetExistence("method.mb.orbBased.exactFieldFree"//trim(FORMER_KEYS(i)))) then
        write (*, '(3A)') " method.mb.orbBased.exactFieldFree", trim(FORMER_KEYS(i)), " is now in the fieldFree block: &
                          &fieldFree.eigen (cut, gpuQ) on Ylm grids, fieldFree.cayley (tolerance, krylovDim, &
                          &maxRestarts, meanFieldQ) on grid.fem; see src/FieldFree/A_FieldFree.md"
        error stop "method.mb.orbBased: the exactFieldFree* keys moved to the fieldFree block"
      end if
    end do

    ! Set by the static-orbital methods only
    Method_Mb_OrbBased_staticOrbsQ = .false.
    ! Set by the methods that implement the option (checked after the dispatch)
    Method_Mb_OrbBased_implementsMeanFieldReferenceQ = .false.
    Method_Mb_OrbBased_implementsSpectralCutQ = .false.
    Method_Mb_OrbBased_implementsSwitchingQ = .false.

    ! Adiabatic switching of the interaction (S_Method_Mb_OrbBased_Switching)
    Method_Mb_OrbBased_switchTime = 0.0_R64
    if (Json_GetExistence("method.mb.orbBased.adiabaticSwitching")) then
      Method_Mb_OrbBased_switchTime = Json_Get("method.mb.orbBased.adiabaticSwitching.switchTime", 0.0_R64)
      if (.not. (Method_Mb_OrbBased_switchTime > 0.0_R64)) &
        error stop "method.mb.orbBased.adiabaticSwitching.switchTime must be positive"
      if (.not. Json_GetExistence("sysInteraction")) &
        error stop "method.mb.orbBased.adiabaticSwitching requires sysInteraction (the interaction it switches on)"
      if (any(Method_Mb_bodyStatistics .ne. 'f')) &
        error stop "method.mb.orbBased.adiabaticSwitching requires fermions (the reference is a Hartree-Fock mean field)"
      if (FieldFree_meanFieldQ) &
        error stop "method.mb.orbBased.adiabaticSwitching is not available with fieldFree.cayley.meanFieldQ"
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

    else if (Json_GetExistence("method.mb.orbBased.tdks")) then
      call Method_Mb_OrbBased_Tdks_Fabricate

    else
      error stop "method.mb.orbBased is missing one of: tdci, tdhx, mctdhx, td2rdm, td2rdmStaticOrbs, tdks"
    end if

    !------------------------------------
    ! options that only some methods implement
    !------------------------------------

    ! The method declares them at its fabricate: the exact part by binding
    ! Method_PropagateExactPart, the others by their implements*Q flags
    if (associated(FieldFree_Propagate) .and. .not. associated(Method_PropagateExactPart)) &
      error stop "fieldFree (the exact field-free part) is not implemented for the selected method"
    if (FieldFree_meanFieldQ .and. .not. Method_Mb_OrbBased_implementsMeanFieldReferenceQ) &
      error stop "fieldFree.cayley.meanFieldQ is not implemented for the selected method"
    if (FieldFree_cutQ .and. .not. Method_Mb_OrbBased_implementsSpectralCutQ) &
      error stop "fieldFree.eigen.cut is not implemented for the selected method"
    if (Method_Mb_OrbBased_switchTime > 0.0_R64 .and. .not. Method_Mb_OrbBased_implementsSwitchingQ) &
      error stop "method.mb.orbBased.adiabaticSwitching is not implemented for the selected method"

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine OrthonormalizeKeepingStateUnsupported(state)
    use M_Utils_UnusedVariables
    complex(R64), intent(inout), contiguous, target :: state(:)

    if (.false.) call UnusedVariables_Mark(state)
    error stop "Method_Mb_OrbBased_OrthonormalizeKeepingState (absorber.keepOrthonormalQ) requires method.mb.orbBased.mctdhx"

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Mb_OrbBased_Setup

    if (allocated(h0Static)) deallocate (h0Static)
    if (ylmNodesQ) call SetupYlmNodes

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
  module subroutine Method_Mb_OrbBased_FillStateWeights(weights)
    use M_Grid
    use M_Orbs
    use M_Method

    real(R64), intent(out), allocatable :: weights(:)

    integer(I32) :: nHead, iOrb

    allocate (weights(size(Method_state)))
    weights(:) = 1.0_R64

    ! Static orbitals are not part of the state
    if (Method_Mb_OrbBased_staticOrbsQ) return

    nHead = size(Method_state) - Grid_nPoints * Orbs_nOrbsInState
    do iOrb = 1, Orbs_nOrbsInState
      weights(nHead + (iOrb - 1) * Grid_nPoints + 1:nHead + iOrb * Grid_nPoints) = Grid_weights(:)
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Mb_OrbBased_ApplyCiHamiltonian(dCoeffs, coeffs, time)
    use M_Utils_UnusedVariables
    use M_Coeffs

    complex(R64), intent(out), contiguous, target :: dCoeffs(:)
    complex(R64), intent(in), contiguous, target  :: coeffs(:)
    real(R64), intent(in) :: time

    if (.false.) call UnusedVariables_Mark(time)

    dCoeffs(:) = 0.0_R64
    call Coeffs_ApplyH1FillRdm1(coeffs, dCoeffs_=dCoeffs, h1_=Method_Mb_OrbBased_h1)
    call Coeffs_ApplyH2FillRdm2(coeffs, dCoeffs_=dCoeffs, h2_=Method_Mb_OrbBased_h2)

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
  module subroutine Method_Mb_OrbBased_FillH0(h0, orbs, time)
    complex(R64), intent(out), allocatable :: h0(:, :)
    complex(R64), intent(in), contiguous   :: orbs(:, :)
    real(R64), intent(in)                  :: time

    complex(R64), allocatable :: hOrbs(:, :)

    allocate (hOrbs, mold=orbs)
    hOrbs(:, :) = 0.0_R64
    call Method_Mb_OrbBased_ApplyKineticOp(hOrbs, orbs, time)
    call Method_Mb_OrbBased_ApplyPotentialOp(hOrbs, orbs, time)
    call Method_Mb_OrbBased_FillOneBodyMatrix(h0, orbs, hOrbs)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Mb_OrbBased_RefillStaticH1(h1, orbs, time)
    use M_SysKinetic
    use M_SysPotential

    complex(R64), intent(inout), allocatable :: h1(:, :)
    complex(R64), intent(in), contiguous     :: orbs(:, :)
    real(R64), intent(in)                    :: time

    complex(R64), allocatable :: gOrbs(:, :), h1Gauge(:, :)

    if (Method_Mb_OrbBased_H1TimeIndependentQ()) return
    if (.not. (SysKinetic_timeIndependentQ .and. SysPotential_timeIndependentQ)) then
      call Method_Mb_OrbBased_FillH1(h1, orbs, time)
    else
      ! the static part: h1 as the setup (FillH0) or a program left it
      if (.not. allocated(h0Static)) h0Static = h1
      h1 = h0Static
      if (gaugeQ) then
        allocate (gOrbs, mold=orbs)
        gOrbs(:, :) = 0.0_R64
        call Method_Mb_OrbBased_ApplyGaugeOp(gOrbs, orbs, time)
        call Method_Mb_OrbBased_FillOneBodyMatrix(h1Gauge, orbs, gOrbs)
        h1 = h1 + h1Gauge
      end if
    end if
    call Method_Mb_OrbBased_AddSwitchingReferenceMatrix(h1, time)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module function Method_Mb_OrbBased_H1TimeIndependentQ() result(res)
    !---------------------------------------------------------------------------
    ! h1 = T + V_ext (+ Ĥ_A if the gauge coupling is bound, gaugeQ) of
    ! ApplySingleBodyOp is time-independent if all of its terms are. The flags
    ! are read at call time, so the result does not depend on the order of
    ! the Sys* fabricates.
    !---------------------------------------------------------------------------
    use M_SysKinetic
    use M_SysPotential
    use M_SysGauge

    logical :: res

    res = SysPotential_timeIndependentQ .and. SysKinetic_timeIndependentQ .and. &
          (.not. gaugeQ .or. SysGauge_timeIndependentQ) .and. Method_Mb_OrbBased_switchTime <= 0.0_R64

  end function

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
  module subroutine Method_Mb_OrbBased_FillOneBodyMatrix(h1, orbs, hOrbs)
    !---------------------------------------------------------------------------
    ! h1(i,j) = ⟨φ_i|hOrbs(:,j)⟩ (FillOrbitalMatrix: equal body type) in the
    ! nOrbsSum layout; restricted: duplicated for the second spin block.
    !---------------------------------------------------------------------------
    use M_Orbs

    complex(R64), intent(out), allocatable :: h1(:, :)
    complex(R64), intent(in), contiguous   :: orbs(:, :)
    complex(R64), intent(in), contiguous   :: hOrbs(:, :)

    complex(R64), allocatable :: h1Spatial(:, :)
    integer(I32) :: nOS, nO

    nOS = Orbs_nOrbsInState
    nO = Method_Mb_OrbBased_nOrbsSum
    call FillOrbitalMatrix(h1Spatial, orbs, hOrbs)

    allocate (h1(nO, nO))
    h1(:, :) = 0.0_R64
    h1(:nOS, :nOS) = h1Spatial(:, :)
    if (Orbs_restrictedQ) h1(nOS + 1:, nOS + 1:) = h1Spatial(:, :)

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
    real(R64) :: lambda
    complex(R64), allocatable :: src(:)

    ! Source density ρ(r₂) = φ*_{i2}(r₂) φ_{j2}(r₂), times the switching
    ! factor of the interaction (S_Method_Mb_OrbBased_Switching)
    call SysInteraction_FillInteractionSrc(src, orbs(:, i2), orbs(:, j2))
    lambda = Method_Mb_OrbBased_SwitchFactor(time)
    if (lambda < 1.0_R64) src(:) = lambda * src(:)

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

    allocate (rho2(nOS, nOS, nOS, nOS))
    if (Orbs_restrictedQ) then
      rho2(:, :, :, :) = rdm2(:nOS, :nOS, :nOS, :nOS) + rdm2(nOS + 1:, :nOS, nOS + 1:, :nOS) + &
                         rdm2(:nOS, nOS + 1:, :nOS, nOS + 1:) + rdm2(nOS + 1:, nOS + 1:, nOS + 1:, nOS + 1:)
    else
      rho2(:, :, :, :) = rdm2(:nOS, :nOS, :nOS, :nOS)
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
    ! restricted; Method_Mb_OrbBased_FillOneBodyMatrix).
    !---------------------------------------------------------------------------
    complex(R64), intent(inout), contiguous, target :: dOrbs(:, :)
    complex(R64), intent(in), contiguous, target    :: orbs(:, :)
    real(R64), intent(in)               :: time
    complex(R64), intent(out), allocatable, optional :: h1_(:, :)

    call Method_Mb_OrbBased_ApplyKineticOp(dOrbs, orbs, time)
    call Method_Mb_OrbBased_ApplyPotentialOp(dOrbs, orbs, time)
    call ApplyDrivingOp(dOrbs, orbs, time)

    if (present(h1_)) call Method_Mb_OrbBased_FillOneBodyMatrix(h1_, orbs, dOrbs)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ApplyDrivingOp(dOrbs, orbs, time)
    !---------------------------------------------------------------------------
    ! The time-dependent one-body terms outside h0 = T + V_ext: the gauge
    ! coupling and the switching reference (1 - lambda) U0 (moving orbitals;
    ! S_Method_Mb_OrbBased_Switching). The linear split-step operator with
    ! an exact field-free part.
    !---------------------------------------------------------------------------
    complex(R64), intent(inout), contiguous :: dOrbs(:, :)
    complex(R64), intent(in), contiguous    :: orbs(:, :)
    real(R64), intent(in)                   :: time

    call Method_Mb_OrbBased_ApplyGaugeOp(dOrbs, orbs, time)
    call Method_Mb_OrbBased_AddSwitchingReference(dOrbs, orbs, time)

  end subroutine

end submodule
