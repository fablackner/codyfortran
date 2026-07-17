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

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Mb_OrbBased_Fabricate
    !---------------------------------------------------------------------------
    ! Initializes orbital-based infrastructure and dispatches to concrete method.
    !
    ! JSON keys read:
    !   - method.mb.orbBased.nOrbs:                 Orbitals per body type
    !   - method.mb.orbBased.regularizationParameter: Numerical stability (default 1e-10)
    !
    ! Binds default operator-application procedures, then branches to the
    ! selected method (TDCI, TDHX, MCTDHX, TD-2RDM, or TD-2RDM-StaticOrbs).
    !---------------------------------------------------------------------------
    use M_Utils_Json
    use M_Utils_Say
    use M_Utils_NoOpProcedures
    use M_Method
    use M_Method_Mb
    use M_Method_Mb_OrbBased_Tdci
    use M_Method_Mb_OrbBased_Tdhx
    use M_Method_Mb_OrbBased_Mctdhx

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
    ! bind default procedure pointers
    !------------------------------------

    Method_GetEnergy => GetEnergy
    Method_Mb_OrbBased_FillH1 => FillH1
    Method_Mb_OrbBased_FillH2 => FillH2
    Method_Mb_OrbBased_ApplyKineticOp => ApplyKineticOp
    Method_Mb_OrbBased_ApplyCorrelationOp => ApplyCorrelationOp
    Method_Mb_OrbBased_ApplyHartreeFockOp => ApplyHartreeFockOp
    Method_Mb_OrbBased_ApplySingleBodyOp => ApplySingleBodyOp

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

    !------------------------------------
    ! branch: method selection
    !------------------------------------

    if (Json_GetExistence("method.mb.orbBased.tdci")) then
      call Method_Mb_OrbBased_Tdci_Fabricate

    else if (Json_GetExistence("method.mb.orbBased.tdhx")) then
      call Method_Mb_OrbBased_Tdhx_Fabricate

    else if (Json_GetExistence("method.mb.orbBased.mctdhx")) then
      call Method_Mb_OrbBased_Mctdhx_Fabricate

    else
      error stop "method.mb.orbBased is missing one of: tdci, tdhx, mctdhx"
    end if

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
  subroutine ApplyCorrelationOp(dOrbs, orbs, rdm1, rdm2, time, h2_)
    !---------------------------------------------------------------------------
    ! Applies the beyond-mean-field correlation operator to orbitals.
    !
    ! This implements the MCTDHX-style correlation contribution arising from
    ! the 2-RDM structure beyond what Hartree-Fock captures. The algorithm:
    !   1. Diagonalize 1-RDM → natural orbitals and occupations
    !   2. Transform 2-RDM to natural orbital basis
    !   3. Apply regularized inverse occupations to avoid small-denominator issues
    !   4. Contract with interaction operator to accumulate in dOrbs
    !
    ! Optionally returns the two-body Hamiltonian matrix h2_ if requested.
    !---------------------------------------------------------------------------
    use M_Utils_RdmDiagonalize
    use M_Method
    use M_Grid
    use M_SysInteraction
    use M_Orbs

    complex(R64), intent(inout), contiguous, target :: dOrbs(:, :)
    complex(R64), intent(in), contiguous, target    :: orbs(:, :)
    complex(R64), intent(in), contiguous             :: rdm1(:, :)
    complex(R64), intent(in), contiguous             :: rdm2(:, :, :, :)
    real(R64), intent(in)                  :: time
    complex(R64), intent(out), allocatable, optional :: h2_(:, :, :, :)

    real(R64), allocatable :: natocc(:)
    complex(R64), allocatable :: natorb(:, :)
    complex(R64), allocatable :: correlOp(:, :, :, :)

    complex(R64), allocatable :: orbTmps(:, :), orbTmpsSwap(:, :)

    complex(R64), allocatable :: rdm1spatial(:, :)
    complex(R64), allocatable :: rdm2spatial(:, :, :, :)

    integer(I32) :: nO, nOS
    integer(I32) :: j1, k1, k2, l2
    integer(I32) :: k2bt, l2bt
    real(R64) :: reg

    nOS = Orbs_nOrbsInState
    nO = Method_Mb_OrbBased_nOrbsSum

    if (present(h2_)) then
      if (.not. allocated(h2_)) allocate (h2_(nO, nO, nO, nO))
      h2_(:, :, :, :) = 0.0_R64
    end if

    !------------------------------------
    ! Build spatial RDMs (spin-summed for restricted case)
    !------------------------------------

    allocate (rdm1spatial(nOS, nOS))
    allocate (rdm2spatial(nOS, nOS, nOS, nOS))

    if (nO .eq. nOS) then
      ! Unrestricted or single-species: direct copy
      rdm1spatial(:, :) = rdm1(1:nOS, 1:nOS)
      rdm2spatial(:, :, :, :) = rdm2(1:nOS, 1:nOS, 1:nOS, 1:nOS)

    else if (nO .eq. 2 * nOS) then
      ! Restricted: sum over spin blocks (αα + αβ + βα + ββ)
      rdm1spatial(:, :) = rdm1(1:nOS, 1:nOS) + rdm1(nOS + 1:, nOS + 1:)
      rdm2spatial(:, :, :, :) = rdm2(1:nOS, 1:nOS, 1:nOS, 1:nOS) + &
                                rdm2(nOS + 1:, :nOS, nOS + 1:, :nOS) + &
                                rdm2(:nOS, nOS + 1:, :nOS, nOS + 1:) + &
                                rdm2(nOS + 1:, nOS + 1:, nOS + 1:, nOS + 1:)
    else
      error stop "invalid nO to nOS relation in ApplyCorrelationOp"
    end if

    !------------------------------------
    ! Diagonalize 1-RDM → natural orbitals
    !------------------------------------

    call RdmDiagonalize_Rdm1(natocc, natorb, rdm1spatial)

    reg = Method_Mb_OrbBased_regularizationParameter

    !------------------------------------
    ! Transform 2-RDM and apply regularized inverse
    !------------------------------------

    allocate (correlOp, mold=rdm2spatial)
    allocate (orbTmps, mold=orbs)
    allocate (orbTmpsSwap, mold=orbs)

    ! The transform to natural orbitals, regularized inversion, and
    ! back-transform combine into a single matrix acting on index 3:
    !   M_{k₁,j₁} = Σ_j U_{k₁,j} · n_j/(n_j² + ε) · U†_{j,j₁}
    block
      complex(R64) :: regInvMat(nOS, nOS)

      do j1 = 1, nOS
        do k1 = 1, nOS
          regInvMat(k1, j1) = sum(natorb(k1, :) * (natocc(:) / (natocc(:)**2 + reg)) * conjg(natorb(j1, :)))
        end do
      end do

      ! Apply: rdm2spatial_{i₁,i₂,j₁,l₂} = Σ_k correlOp_{i₁,i₂,k,l₂} · M_{k,j₁}
      correlOp(:, :, :, :) = rdm2spatial(:, :, :, :)
      rdm2spatial(:, :, :, :) = 0.0_R64
      do j1 = 1, nOS
        do k1 = 1, nOS
          rdm2spatial(:, :, j1, :) = rdm2spatial(:, :, j1, :) + correlOp(:, :, k1, :) * regInvMat(k1, j1)
        end do
      end do
    end block

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

    deallocate (rdm1spatial)
    deallocate (rdm2spatial)

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
