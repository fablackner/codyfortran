! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Time-dependent Kohn-Sham method: spin blocks, Kohn-Sham operator and energy.
!>
!> Spin blocks (SetupSpinBlocks): ranges of stored orbitals sharing one
!> density column and one exchange operator. Restricted storage: one block
!> (all stored orbitals), density weights R↑ + R↓, exchange weights R↑ and
!> multiplicity 2 (both spins); unrestricted: one block per body type with
!> its 1-RDM block R_bt for both. The Kohn-Sham operator of a block is
!> v_KS,s (SysXc, one call for all spins) plus -a K_s with
!> K_s φ_i = Σ_l V[conj(φ_l) φ_i] g_l, g = φ R (pair potentials of
!> SysInteraction over unordered pairs, the swapped one by conjugation if
!> SysInteraction_conjSymmetricQ; pairs of two unoccupied orbitals skipped).
submodule(M_Method_Mb_OrbBased_Tdks) S_Method_Mb_OrbBased_Tdks

  implicit none

  !=============================================================================
  ! local data
  !=============================================================================

  !> Static CI coefficients (the occupations of the determinant)
  complex(R64), allocatable, target :: staticComponents(:)

  !> Spin blocks: stored-orbital ranges, density column and the 1-RDM blocks
  !> rdm1(r0:, r0:), r0 = rdmStart(k, b), k ≤ nRdmBlocks(b), whose sum weights
  !> the density (restricted: R↑ + R↓); the first one weights the exchange
  integer(I32) :: nBlocks = 0
  integer(I32), allocatable :: blockStart(:), blockEnd(:), blockSpin(:)
  integer(I32), allocatable :: rdmStart(:, :), nRdmBlocks(:)

  !> Energy multiplicity of the exchange of a block (2: restricted storage)
  real(R64) :: exchangeMultiplicity = 1.0_R64

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Mb_OrbBased_Tdks_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Method
    use M_Method_Mb
    use M_Method_Mb_OrbBased
    use M_FieldFree

    call Say_Fabricate("method.mb.orbBased.tdks")

    if (.not. Json_GetExistence("sysXc")) &
      error stop "method.mb.orbBased.tdks requires sysXc (the exchange-correlation functional)"
    if (any(Method_Mb_bodyStatistics .ne. 'f')) error stop "method.mb.orbBased.tdks requires fermions"
    if (Method_Mb_nBodyTypes > 2) error stop "method.mb.orbBased.tdks supports one or two body types (spins)"

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    Method_kohnShamQ = .true.
    Method_Mb_OrbBased_implementsMeanFieldReferenceQ = .true.
    Method_Mb_OrbBased_implementsSpectralCutQ = .true.
    Method_Setup => Setup
    Method_GetEnergy => GetEnergy
    Method_Mb_OrbBased_FillRdm1 => FillRdm1
    Method_Mb_OrbBased_FillRdm2 => FillRdm2
    Method_TimeDerivative => TimeDerivative
    Method_Mb_OrbBased_TimeDerivativeOrbsLin => TimeDerivativeOrbsLin
    Method_Mb_OrbBased_TimeDerivativeCoeffsPlusOrbsNonLin => TimeDerivativeOrbsNonLin
    ! The exact part propagates h0 on the orbitals; phi functions only for
    ! the eigen-decomposed exact part (fieldFree.eigen)
    if (associated(FieldFree_Propagate)) then
      Method_PropagateExactPart => PropagateExactPart
      if (associated(FieldFree_ApplyPhi)) Method_ApplyExactPartPhi => ApplyExactPartPhi
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Setup
    use M_Utils_Say
    use M_Grid
    use M_Orbs
    use M_OrbsInit
    use M_Coeffs
    use M_CoeffsInit
    use M_SysXc
    use M_Method
    use M_Method_Mb_OrbBased

    integer(I32) :: nG, nOS

    call Say_Setup("method.mb.orbBased.tdks")
    call Method_Mb_OrbBased_Setup

    if (SysXc_nPoints .le. 0 .or. .not. allocated(SysXc_weights)) &
      error stop "method.mb.orbBased.tdks: call SysXc_Fabricate and SysXc_Setup before Method_Setup"

    nG = Grid_nPoints
    nOS = Orbs_nOrbsInState

    allocate (staticComponents(Coeffs_nCoeffs))
    allocate (Method_state(nG * nOS))

    Coeffs_coeffs(1:Coeffs_nCoeffs) => staticComponents(1:Coeffs_nCoeffs)
    call CoeffsInit_Initialize(Coeffs_coeffs)

    Orbs_orbs(1:nG, 1:nOS) => Method_state(1:nG * nOS)
    call OrbsInit_Initialize(Orbs_orbs)

    call Coeffs_ApplyH1FillRdm1(Coeffs_coeffs, rdm1_=Method_Mb_OrbBased_rdm1)
    call Coeffs_ApplyH2FillRdm2(Coeffs_coeffs, rdm2_=Method_Mb_OrbBased_rdm2)

    call SetupSpinBlocks(Method_Mb_OrbBased_rdm1)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Spin blocks and spin channels of the stored orbitals (submodule header);
  !> restricted storage requires equal spin blocks of the 1-RDM.
  subroutine SetupSpinBlocks(rdm1)
    use M_Orbs
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    complex(R64), intent(in), contiguous :: rdm1(:, :)

    integer(I32) :: nOS, bt

    nOS = Orbs_nOrbsInState
    if (allocated(blockStart)) deallocate (blockStart, blockEnd, blockSpin, rdmStart, nRdmBlocks)

    if (Orbs_restrictedQ) then
      if (Method_Mb_nBodyTypes .ne. 2) error stop "method.mb.orbBased.tdks: restricted orbitals need two body types (spins)"
      if (maxval(abs(rdm1(:nOS, :nOS) - rdm1(nOS + 1:, nOS + 1:))) > 1.0e-10_R64) &
        error stop "method.mb.orbBased.tdks: restricted orbitals need a closed-shell state (equal spin blocks of the 1-RDM)"
      Method_Mb_OrbBased_Tdks_nSpin = 1
      nBlocks = 1
      exchangeMultiplicity = 2.0_R64
      allocate (blockStart(1), blockEnd(1), blockSpin(1), rdmStart(2, 1), nRdmBlocks(1))
      blockStart(1) = 1
      blockEnd(1) = nOS
      blockSpin(1) = 1
      nRdmBlocks(1) = 2
      rdmStart(:, 1) = [1, nOS + 1]
    else
      Method_Mb_OrbBased_Tdks_nSpin = 2
      nBlocks = Method_Mb_nBodyTypes
      exchangeMultiplicity = 1.0_R64
      allocate (blockStart(nBlocks), blockEnd(nBlocks), blockSpin(nBlocks), rdmStart(2, nBlocks), nRdmBlocks(nBlocks))
      do bt = 1, nBlocks
        blockStart(bt) = Method_Mb_OrbBased_nOrbsStart(bt)
        blockEnd(bt) = Method_Mb_OrbBased_nOrbsEnd(bt)
        blockSpin(bt) = bt
        nRdmBlocks(bt) = 1
        rdmStart(:, bt) = blockStart(bt)
      end do
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine TimeDerivative(dState, state, time)
    use M_Utils_Constants
    use M_Method_Mb_OrbBased
    use M_FieldFree

    complex(R64), intent(out), contiguous, target :: dState(:)
    complex(R64), intent(in), contiguous, target  :: state(:)
    real(R64), intent(in)             :: time

    complex(R64), contiguous, pointer :: orbs(:, :), dOrbs(:, :)

    dState(:) = 0.0_R64
    call GetOrbs(orbs, state)
    call GetOrbs(dOrbs, dState)

    call Method_Mb_OrbBased_ApplySingleBodyOp(dOrbs, orbs, time)
    call Method_Mb_OrbBased_Tdks_ApplyKohnShamOp(dOrbs, orbs, Method_Mb_OrbBased_rdm1, time)
    call FieldFree_ProjectBelowCut(dOrbs)

    dState(:) = -IU * dState(:) ! out is the time-derivative!

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Linear split-step part: -i L orbs with L = Method_Mb_OrbBased_ApplyLinOp
  !> (h, or only the gauge coupling with an exact field-free part).
  subroutine TimeDerivativeOrbsLin(dState, state, time)
    use M_Utils_Constants
    use M_Method_Mb_OrbBased
    use M_FieldFree

    complex(R64), intent(out), contiguous, target :: dState(:)
    complex(R64), intent(in), contiguous, target  :: state(:)
    real(R64), intent(in)             :: time

    complex(R64), contiguous, pointer :: orbs(:, :), dOrbs(:, :)

    dState(:) = 0.0_R64
    call GetOrbs(orbs, state)
    call GetOrbs(dOrbs, dState)

    call Method_Mb_OrbBased_ApplyLinOp(dOrbs, orbs, time)
    call FieldFree_ProjectBelowCut(dOrbs)

    dState(:) = -IU * dState(:) ! out is the time-derivative!

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Non-linear split-step part: -i (v_KS - a K) orbs (the mean field of the
  !> instantaneous densities); with the frozen mean-field reference
  !> (FieldFree_meanFieldQ), whose U = v_KS[rho0] - a K0 the exact part
  !> propagates with h0, -i (v_KS - a K - U) orbs, which vanishes at the
  !> reference state.
  subroutine TimeDerivativeOrbsNonLin(dState, state, time)
    use M_Utils_Constants
    use M_Method_Mb_OrbBased
    use M_FieldFree

    complex(R64), intent(out), contiguous, target :: dState(:)
    complex(R64), intent(in), contiguous, target  :: state(:)
    real(R64), intent(in)             :: time

    complex(R64), contiguous, pointer :: orbs(:, :), dOrbs(:, :)
    complex(R64), allocatable :: uOrbs(:, :)

    dState(:) = 0.0_R64
    call GetOrbs(orbs, state)
    call GetOrbs(dOrbs, dState)

    call Method_Mb_OrbBased_Tdks_ApplyKohnShamOp(dOrbs, orbs, Method_Mb_OrbBased_rdm1, time)

    if (FieldFree_meanFieldQ) then
      call EnsureFieldFreeMeanField(state, time)
      allocate (uOrbs, mold=orbs)
      uOrbs(:, :) = (0.0_R64, 0.0_R64)
      call FieldFree_ApplyMeanField(uOrbs, orbs)
      dOrbs(:, :) = dOrbs(:, :) - uOrbs(:, :)
    end if
    call FieldFree_ProjectBelowCut(dOrbs)

    dState(:) = -IU * dState(:) ! out is the time-derivative!

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Builds the frozen Kohn-Sham reference of the exact part from state (its
  !> orbitals and the 1-RDM) unless it exists since the last FieldFree_Setup:
  !> at the first exact-part call (or the first nonlinear call if a
  !> propagator starts with that part), unless a program built it before
  !> (FieldFree_BuildMeanField, e.g. from the ground state
  !> before a kick).
  subroutine EnsureFieldFreeMeanField(state, time)
    use M_Method_Mb_OrbBased
    use M_FieldFree

    complex(R64), intent(in), contiguous, target :: state(:)
    real(R64), intent(in) :: time

    complex(R64), contiguous, pointer :: orbs(:, :)

    if (FieldFree_meanFieldBuiltQ) return

    call GetOrbs(orbs, state)
    call FieldFree_BuildMeanField(orbs, Method_Mb_OrbBased_rdm1, time)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Exact field-free propagation of the orbitals (FieldFree_Propagate).
  subroutine PropagateExactPart(state, t0, t1, imagTimeQ_)
    use M_Method_Mb_OrbBased
    use M_FieldFree

    complex(R64), intent(inout), contiguous, target :: state(:)
    real(R64), intent(in) :: t0
    real(R64), intent(in) :: t1
    logical, intent(in), optional :: imagTimeQ_

    complex(R64), contiguous, pointer :: orbs(:, :)
    logical :: imagTimeQ

    imagTimeQ = .false.
    if (present(imagTimeQ_)) imagTimeQ = imagTimeQ_

    if (FieldFree_meanFieldQ .and. .not. imagTimeQ) call EnsureFieldFreeMeanField(state, t0)

    call GetOrbs(orbs, state)
    call FieldFree_Propagate(orbs, t1 - t0, imagTimeQ)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Phi functions of the exact part (h0 on the orbitals, the whole state).
  subroutine ApplyExactPartPhi(res, vecs, dt)
    use M_Grid
    use M_Orbs
    use M_Method_Mb_OrbBased
    use M_FieldFree

    complex(R64), intent(out), contiguous, target :: res(:)
    complex(R64), intent(in), contiguous, target :: vecs(:, :)
    real(R64), intent(in) :: dt

    complex(R64), contiguous, pointer :: orbs(:, :), orbsM(:, :)
    complex(R64), allocatable :: orbsList(:, :, :)
    integer(I32) :: m

    call GetOrbs(orbs, res)
    allocate (orbsList(Grid_nPoints, Orbs_nOrbsInState, 0:size(vecs, 2) - 1))
    do m = 0, size(vecs, 2) - 1
      call GetOrbs(orbsM, vecs(:, m + 1))
      orbsList(:, :, m) = orbsM(:, :)
    end do

    call FieldFree_ApplyPhi(orbs, orbsList, dt)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Mb_OrbBased_Tdks_ApplyKohnShamOp(dOrbs, orbs, rdm1, time)
    complex(R64), intent(inout), contiguous :: dOrbs(:, :)
    complex(R64), intent(in), contiguous :: orbs(:, :)
    complex(R64), intent(in), contiguous :: rdm1(:, :)
    real(R64), intent(in) :: time

    call EvaluateKohnSham(orbs, rdm1, time, dOrbs_=dOrbs)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Kohn-Sham energy of Method_state: Σ_b Σ_i ⟨g_i|h φ_i⟩ + E_H + E_xc - a E_x.
  function GetEnergy(time) result(res)
    use M_Grid
    use M_Orbs
    use M_SysXc
    use M_Method
    use M_Method_Mb_OrbBased

    real(R64) :: res
    real(R64), intent(in) :: time

    complex(R64), allocatable :: hOrbs(:, :), gD(:, :)
    real(R64) :: oneBody
    integer(I32) :: b, i

    call FillRdm1(Method_Mb_OrbBased_rdm1, Method_state)

    allocate (hOrbs, mold=Orbs_orbs)
    hOrbs(:, :) = 0.0_R64
    call Method_Mb_OrbBased_ApplySingleBodyOp(hOrbs, Orbs_orbs, time)

    oneBody = 0.0_R64
    do b = 1, nBlocks
      call FillBlockWeights(gD, Orbs_orbs, Method_Mb_OrbBased_rdm1, b, nRdmBlocks(b))
      do i = blockStart(b), blockEnd(b)
        oneBody = oneBody + real(Grid_InnerProduct(gD(:, i - blockStart(b) + 1), hOrbs(:, i)), kind=R64)
      end do
    end do

    call EvaluateKohnSham(Orbs_orbs, Method_Mb_OrbBased_rdm1, time)

    Method_Mb_OrbBased_Tdks_oneBodyEnergy = oneBody
    res = oneBody + Method_Mb_OrbBased_Tdks_hartreeEnergy + Method_Mb_OrbBased_Tdks_xcEnergy - &
          SysXc_exactExchangeFraction * Method_Mb_OrbBased_Tdks_exactExchangeEnergy

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Densities of the spin blocks, Kohn-Sham potentials, energies and, if
  !> dOrbs_ is present, dOrbs_ += (v_KS,s - a K_s) orbs.
  subroutine EvaluateKohnSham(orbs, rdm1, time, dOrbs_)
    use M_Grid
    use M_Orbs
    use M_SysXc

    complex(R64), intent(in), contiguous :: orbs(:, :)
    complex(R64), intent(in), contiguous :: rdm1(:, :)
    real(R64), intent(in) :: time
    complex(R64), intent(inout), contiguous, optional :: dOrbs_(:, :)

    real(R64), allocatable :: kohnShamPot(:, :)
    complex(R64), allocatable :: g(:, :), kOrbs(:, :), orbTmp(:)
    real(R64) :: hartreeEnergy, xcEnergy, exchangeEnergy, a
    integer(I32) :: b, i, iS, iE

    call FillSpinKohnShamPotentials(kohnShamPot, hartreeEnergy, xcEnergy, orbs, rdm1, time)

    if (present(dOrbs_)) then
      allocate (orbTmp(Grid_nPoints))
      do b = 1, nBlocks
        do i = blockStart(b), blockEnd(b)
          call SysXc_MultiplyWithPotential(orbTmp, kohnShamPot(:, blockSpin(b)), orbs(:, i))
          dOrbs_(:, i) = dOrbs_(:, i) + orbTmp(:)
        end do
      end do
    end if

    !------------------------------------
    ! exact exchange of a hybrid
    !------------------------------------

    exchangeEnergy = 0.0_R64
    a = SysXc_exactExchangeFraction
    if (abs(a) > 0.0_R64) then
      do b = 1, nBlocks
        iS = blockStart(b)
        iE = blockEnd(b)
        call FillBlockWeights(g, orbs, rdm1, b, 1)
        allocate (kOrbs(Grid_nPoints, iE - iS + 1))
        call ApplyExchange(kOrbs, orbs(:, iS:iE), g, time)
        do i = 1, iE - iS + 1
          exchangeEnergy = exchangeEnergy + 0.5_R64 * exchangeMultiplicity * &
                           real(Grid_InnerProduct(g(:, i), kOrbs(:, i)), kind=R64)
        end do
        if (present(dOrbs_)) dOrbs_(:, iS:iE) = dOrbs_(:, iS:iE) - a * kOrbs(:, :)
        deallocate (kOrbs)
      end do
    end if

    Method_Mb_OrbBased_Tdks_hartreeEnergy = hartreeEnergy
    Method_Mb_OrbBased_Tdks_xcEnergy = xcEnergy
    Method_Mb_OrbBased_Tdks_exactExchangeEnergy = exchangeEnergy

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Mb_OrbBased_Tdks_FillKohnShamPotential(kohnShamPot, spinOfBodyType, orbs, rdm1, time)
    use M_Method_Mb

    real(R64), intent(out), allocatable :: kohnShamPot(:, :)
    integer(I32), intent(out), allocatable :: spinOfBodyType(:)
    complex(R64), intent(in), contiguous :: orbs(:, :)
    complex(R64), intent(in), contiguous :: rdm1(:, :)
    real(R64), intent(in) :: time

    real(R64) :: hartreeEnergy, xcEnergy
    integer(I32) :: bt

    call FillSpinKohnShamPotentials(kohnShamPot, hartreeEnergy, xcEnergy, orbs, rdm1, time)

    ! Column of every body type: restricted storage uses the total density
    ! (one column) for both spins
    allocate (spinOfBodyType(Method_Mb_nBodyTypes))
    do bt = 1, Method_Mb_nBodyTypes
      spinOfBodyType(bt) = min(bt, Method_Mb_OrbBased_Tdks_nSpin)
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Spin densities of the blocks and their local Kohn-Sham potentials
  !> v_KS,s (SysXc_nPoints, nSpin) with the Hartree and xc energies.
  subroutine FillSpinKohnShamPotentials(kohnShamPot, hartreeEnergy, xcEnergy, orbs, rdm1, time)
    use M_SysXc

    real(R64), intent(out), allocatable :: kohnShamPot(:, :)
    real(R64), intent(out) :: hartreeEnergy
    real(R64), intent(out) :: xcEnergy
    complex(R64), intent(in), contiguous :: orbs(:, :)
    complex(R64), intent(in), contiguous :: rdm1(:, :)
    real(R64), intent(in) :: time

    real(R64), allocatable :: density(:, :)
    complex(R64), allocatable :: g(:, :)
    integer(I32) :: b

    if (nBlocks .eq. 0) error stop "method.mb.orbBased.tdks: Method_Setup was not called"

    allocate (density(SysXc_nPoints, Method_Mb_OrbBased_Tdks_nSpin), &
              kohnShamPot(SysXc_nPoints, Method_Mb_OrbBased_Tdks_nSpin))
    density(:, :) = 0.0_R64
    do b = 1, nBlocks
      call FillBlockWeights(g, orbs, rdm1, b, nRdmBlocks(b))
      call SysXc_AddDensity(density(:, blockSpin(b)), orbs(:, blockStart(b):blockEnd(b)), g)
    end do
    call SysXc_FillKohnShamPotential(kohnShamPot, hartreeEnergy, xcEnergy, density, time)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> g = φ_b R of block b, R the sum of its first nRdm 1-RDM blocks (the
  !> density: nRdmBlocks(b), restricted R↑ + R↓; the exchange: 1, R_s).
  subroutine FillBlockWeights(g, orbs, rdm1, b, nRdm)
    complex(R64), intent(out), allocatable :: g(:, :)
    complex(R64), intent(in), contiguous :: orbs(:, :)
    complex(R64), intent(in), contiguous :: rdm1(:, :)
    integer(I32), intent(in) :: b
    integer(I32), intent(in) :: nRdm

    complex(R64), allocatable :: weights(:, :)
    integer(I32) :: n, k, r0

    n = blockEnd(b) - blockStart(b) + 1
    allocate (weights(n, n))
    weights(:, :) = 0.0_R64
    do k = 1, nRdm
      r0 = rdmStart(k, b)
      weights(:, :) = weights(:, :) + rdm1(r0:r0 + n - 1, r0:r0 + n - 1)
    end do
    g = matmul(orbs(:, blockStart(b):blockEnd(b)), weights)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> kOrbs(:, i) = K φ_i = Σ_l V[conj(φ_l) φ_i] g_l within one spin block:
  !> one potential solve per unordered pair (the swapped pair by
  !> conjugation if conjSymmetricQ); pairs of two unweighted orbitals
  !> (g = 0) contribute nothing and are skipped.
  subroutine ApplyExchange(kOrbs, orbs, g, time)
    use M_Grid
    use M_SysInteraction

    complex(R64), intent(out), contiguous :: kOrbs(:, :)
    complex(R64), intent(in), contiguous :: orbs(:, :)
    complex(R64), intent(in), contiguous :: g(:, :)
    real(R64), intent(in) :: time

    complex(R64), allocatable :: src(:), pot(:), potSwap(:), orbTmp(:)
    logical, allocatable :: weightedQ(:)
    integer(I32) :: n, i, l

    n = size(orbs, 2)
    kOrbs(:, :) = 0.0_R64
    allocate (orbTmp(Grid_nPoints), weightedQ(n))
    do l = 1, n
      weightedQ(l) = any(abs(g(:, l)) > 0.0_R64)
    end do

    do i = 1, n
      do l = 1, i
        if (.not. (weightedQ(l) .or. weightedQ(i))) cycle

        ! V_li = V[conj(φ_l) φ_i]: K φ_i += V_li g_l
        call SysInteraction_FillInteractionSrc(src, orbs(:, l), orbs(:, i))
        call SysInteraction_FillInteractionPotential(pot, src, time, bt1_=-1, bt2_=-1)
        if (weightedQ(l)) then
          call SysInteraction_MultiplyWithInteractionPotential(orbTmp, pot, g(:, l))
          kOrbs(:, i) = kOrbs(:, i) + orbTmp(:)
        end if
        if (l .eq. i .or. .not. weightedQ(i)) cycle

        ! V_il = V[conj(φ_i) φ_l]: K φ_l += V_il g_i
        if (SysInteraction_conjSymmetricQ) then
          allocate (potSwap, mold=pot)
          call SysInteraction_ConjugateInteractionPotential(potSwap, pot)
        else
          call SysInteraction_FillInteractionSrc(src, orbs(:, i), orbs(:, l))
          call SysInteraction_FillInteractionPotential(potSwap, src, time, bt1_=-1, bt2_=-1)
        end if
        call SysInteraction_MultiplyWithInteractionPotential(orbTmp, potSwap, g(:, i))
        kOrbs(:, l) = kOrbs(:, l) + orbTmp(:)
        deallocate (potSwap)
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillRdm1(rdm1, state)
    use M_Coeffs
    use M_Utils_UnusedVariables

    complex(R64), intent(out), allocatable :: rdm1(:, :)
    complex(R64), intent(in), contiguous, target :: state(:)

    complex(R64), contiguous, pointer :: coeffs(:)

    if (.false.) call UnusedVariables_Mark(state)
    coeffs(1:Coeffs_nCoeffs) => staticComponents(1:Coeffs_nCoeffs)
    call Coeffs_ApplyH1FillRdm1(coeffs, rdm1_=rdm1)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillRdm2(rdm2, state)
    use M_Coeffs
    use M_Utils_UnusedVariables

    complex(R64), intent(out), allocatable :: rdm2(:, :, :, :)
    complex(R64), intent(in), contiguous, target :: state(:)

    complex(R64), contiguous, pointer :: coeffs(:)

    if (.false.) call UnusedVariables_Mark(state)
    coeffs(1:Coeffs_nCoeffs) => staticComponents(1:Coeffs_nCoeffs)
    call Coeffs_ApplyH2FillRdm2(coeffs, rdm2_=rdm2)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine GetOrbs(orbs, state)
    use M_Grid
    use M_Orbs

    complex(R64), pointer, intent(out) :: orbs(:, :)
    complex(R64), intent(in), contiguous, target :: state(:)

    orbs(1:Grid_nPoints, 1:Orbs_nOrbsInState) => state(1:)

  end subroutine

end submodule
