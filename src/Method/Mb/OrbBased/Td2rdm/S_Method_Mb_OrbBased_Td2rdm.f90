! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Method_Mb_OrbBased_Td2rdm) S_Method_Mb_OrbBased_Td2rdm
  !-----------------------------------------------------------------------------
  ! Active/weak natural-orbital split (nWeakOrbs > 0; derivation, energy
  ! conservation and usage in src/Method/A_Method.md). Orbital EOM
  !   dphi_i/dt = -i Q sum_u G_u M(u,i),  G_u = sum_j h phi_j D1(j,u) + Gamma_u,
  !   M = sum_{u active} |u> n_u/(n_u^2 + eps) <u|.
  ! Split step: the linear part rotates with H_L = (1-P_w) h (1-P_w) + P_w h P_w,
  ! the non-linear part adds the active-weak blocks of h1 to the 2-RDM.
  ! Gamma M is Method_Mb_OrbBased_ApplyCorrelationOp with invRdm1 = M.
  !-----------------------------------------------------------------------------

  implicit none

  !=============================================================================
  ! local data
  !=============================================================================

  !> Weak natural orbitals on the grid from the last refresh, stored in the
  !> leading columns of their body type (layout of Orbs_orbs, so that the
  !> one-body operator applies); only columns with weakColQ hold one.
  complex(R64), allocatable :: weakOrbs(:, :)
  !> Whether column j of weakOrbs holds a weak natural orbital.
  logical, allocatable :: weakColQ(:)
  !> h weakOrbs and <weakOrbs_w|h|weakOrbs_v> at hWeakOrbsTime (the linear
  !> part is evaluated repeatedly at one time, e.g. for every Krylov vector).
  complex(R64), allocatable :: hWeakOrbs(:, :)
  complex(R64), allocatable :: hWeak(:, :)
  real(R64) :: hWeakOrbsTime
  logical :: hWeakOrbsValidQ = .false.

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Mb_OrbBased_Td2rdm_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Orbs
    use M_TwoRdm
    use M_Method
    use M_Method_Mb_OrbBased

    call Say_Fabricate("method.mb.orbbased.td2rdm")

    !------------------------------------
    ! read parameters
    !------------------------------------

    Method_Mb_OrbBased_Td2rdm_nWeakOrbs = Json_Get("method.mb.orbBased.td2rdm.nWeakOrbs", 0)
    if (Method_Mb_OrbBased_Td2rdm_nWeakOrbs < 0 .or. Method_Mb_OrbBased_Td2rdm_nWeakOrbs >= minval(Method_Mb_OrbBased_nOrbs)) &
      error stop "method.mb.orbBased.td2rdm.nWeakOrbs must be >= 0 and < nOrbs"

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    Method_Setup => Setup
    Method_Mb_OrbBased_FillRdm1 => FillRdm1
    Method_Mb_OrbBased_FillRdm2 => FillRdm2
    Method_TimeDerivative => TimeDerivative
    Method_Mb_OrbBased_TimeDerivativeOrbsLin => TimeDerivativeOrbsLin
    Method_Mb_OrbBased_TimeDerivativeCoeffsPlusOrbsNonLin => TimeDerivativeTwoRdmPlusOrbsNonLin

    ! With weak orbitals the linear part is H_L, not h0 plus the gauge coupling
    if (Method_Mb_OrbBased_exactFieldFreeQ) then
      if (Method_Mb_OrbBased_Td2rdm_nWeakOrbs > 0) &
        error stop "method.mb.orbBased.exactFieldFreeQ is not available with td2rdm.nWeakOrbs > 0"
      Method_PropagateExactPart => PropagateExactPart
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Mb_OrbBased_Td2rdm_RefreshWeakOrbs
    use M_TwoRdm
    use M_Orbs

    complex(R64), allocatable :: rdm1(:, :), invRdm1(:, :), activeProj(:, :), weakNatOrbs(:, :)

    if (Method_Mb_OrbBased_Td2rdm_nWeakOrbs .eq. 0) return

    call TwoRdm_FillRdm1FromTwoRdmState(rdm1, TwoRdm_twoRdmState)
    call FillWeakSplit(invRdm1, activeProj, weakNatOrbs, rdm1)
    call StoreWeakOrbs(Orbs_orbs, weakNatOrbs)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Setup
    use M_Utils_Say
    use M_Grid
    use M_Orbs
    use M_OrbsInit
    use M_TwoRdm
    use M_TwoRdmInit
    use M_Method
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    integer(I32) :: nVals1, nVals2, nG, nOS, nE

    nG = Grid_nPoints
    nOS = Orbs_nOrbsInState
    nE = TwoRdm_nElements

    nVals1 = nE
    nVals2 = nG * nOS

    call Say_Setup("method.mb.orbbased.td2rdm")
    call Method_Mb_OrbBased_Setup

    allocate (Method_state(1:(nVals1 + nVals2)))

    call GetTwoRdmState(TwoRdm_twoRdmState, Method_state)
    call TwoRdmInit_Initialize(TwoRdm_twoRdmState)

    call GetOrbs(Orbs_orbs, Method_state)
    call OrbsInit_Initialize(Orbs_orbs)

    ! weak natural orbitals of the initial state for the first linear step
    call Method_Mb_OrbBased_Td2rdm_RefreshWeakOrbs

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine TimeDerivative(dState, state, time)
    use M_Utils_Constants
    use M_Reconstruction
    use M_TwoRdm
    use M_Orbs
    use M_Method
    use M_Method_Mb_OrbBased

    complex(R64), intent(out), contiguous, target :: dState(:)
    complex(R64), intent(in), contiguous, target  :: state(:)
    real(R64), intent(in)             :: time

    complex(R64), contiguous, pointer :: twoRdmState(:)
    complex(R64), contiguous, pointer :: dTwoRdmState(:)

    complex(R64), contiguous, pointer :: orbs(:, :)
    complex(R64), contiguous, pointer :: dOrbs(:, :)

    complex(R64), allocatable :: h1(:, :)
    complex(R64), allocatable :: h2(:, :, :, :)
    complex(R64), allocatable :: rdm1(:, :)
    complex(R64), allocatable :: rdm2(:, :, :, :)
    complex(R64), allocatable :: invRdm1(:, :), activeProj(:, :), weakNatOrbs(:, :), hOrbs(:, :)

    dState(:) = 0.0_R64

    call GetTwoRdmState(twoRdmState, state)
    call GetTwoRdmState(dTwoRdmState, dState)

    call GetOrbs(orbs, state)
    call GetOrbs(dOrbs, dState)

    call TwoRdm_FillRdm1FromTwoRdmState(rdm1, twoRdmState)
    call TwoRdm_FillRdm2FromTwoRdmState(rdm2, twoRdmState)

    if (Method_Mb_OrbBased_Td2rdm_nWeakOrbs .eq. 0) then
      call Method_Mb_OrbBased_ApplySingleBodyOp(dOrbs, orbs, time, h1_=h1)
      call Method_Mb_OrbBased_FillInverseRdm1(invRdm1, rdm1)
    else
      ! h phi D1 M + Gamma M
      call FillWeakSplit(invRdm1, activeProj, weakNatOrbs, rdm1)
      allocate (hOrbs, mold=orbs)
      hOrbs(:, :) = 0.0_R64
      call Method_Mb_OrbBased_ApplySingleBodyOp(hOrbs, orbs, time, h1_=h1)
      dOrbs(:, :) = matmul(hOrbs, activeProj)
    end if
    call Method_Mb_OrbBased_ApplyCorrelationOp(dOrbs, orbs, invRdm1, rdm2, time, h2_=h2)

    call TwoRdm_ApplyCommutatorH1(dTwoRdmState, twoRdmState, h1)
    call TwoRdm_ApplyCommutatorH2(dTwoRdmState, twoRdmState, h2)
    call Reconstruction_AddCollisionOp(dTwoRdmState, twoRdmState, h2)

    call Orbs_ProjectOnSubspace(dOrbs, orbs)

    dState(:) = -IU * dState(:) ! out is the time-derivative!

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Linear split-step part: the one-body operator on the orbitals only (the
  !> 2-RDM is constant), dOrbs = -i L orbs without projection,
  !> L = Method_Mb_OrbBased_ApplyLinOp (h, or only the gauge coupling if
  !> exactFieldFreeQ). With weak orbitals h is replaced by H_L (header).
  subroutine TimeDerivativeOrbsLin(dState, state, time)
    use M_Utils_Constants
    use M_Grid
    use M_TwoRdm
    use M_Orbs
    use M_Method_Mb_OrbBased

    complex(R64), intent(out), contiguous, target :: dState(:)
    complex(R64), intent(in), contiguous, target  :: state(:)
    real(R64), intent(in)             :: time

    complex(R64), contiguous, pointer :: orbs(:, :)
    complex(R64), contiguous, pointer :: dOrbs(:, :)

    dState(:) = 0.0_R64

    call GetOrbs(orbs, state)
    call GetOrbs(dOrbs, dState)

    call Method_Mb_OrbBased_ApplyLinOp(dOrbs, orbs, time)
    if (Method_Mb_OrbBased_Td2rdm_nWeakOrbs > 0) call AddWeakProjectorTerms(dOrbs, orbs, time)

    dState(:) = -IU * dState(:) ! out is the time-derivative!

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Exact field-free propagation of the orbitals (the 2-RDM is constant in
  !> the linear split step).
  subroutine PropagateExactPart(state, t0, t1, imagTimeQ_)
    use M_Method_Mb_OrbBased

    complex(R64), intent(inout), contiguous, target :: state(:)
    real(R64), intent(in) :: t0
    real(R64), intent(in) :: t1
    logical, intent(in), optional :: imagTimeQ_

    complex(R64), contiguous, pointer :: orbs(:, :)
    logical :: imagTimeQ

    imagTimeQ = .false.
    if (present(imagTimeQ_)) imagTimeQ = imagTimeQ_

    call GetOrbs(orbs, state)
    call Method_Mb_OrbBased_PropagateFieldFree(orbs, t1 - t0, imagTimeQ)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Non-linear split-step part: the two-body interaction only. The 2-RDM
  !> evolves with [h2, D2] plus the collision operator, the orbitals with the
  !> projected correlation term. The one-body commutator is absent because the
  !> linear part rotates the orbitals with the full (unprojected) h. With weak
  !> orbitals the linear part rotates with H_L instead, the 2-RDM gets the
  !> active-weak blocks of h1, and the weak orbitals are refreshed (header).
  subroutine TimeDerivativeTwoRdmPlusOrbsNonLin(dState, state, time)
    use M_Utils_Constants
    use M_Grid
    use M_Reconstruction
    use M_TwoRdm
    use M_Orbs
    use M_Method_Mb_OrbBased

    complex(R64), intent(out), contiguous, target :: dState(:)
    complex(R64), intent(in), contiguous, target  :: state(:)
    real(R64), intent(in)             :: time

    complex(R64), contiguous, pointer :: twoRdmState(:)
    complex(R64), contiguous, pointer :: dTwoRdmState(:)

    complex(R64), contiguous, pointer :: orbs(:, :)
    complex(R64), contiguous, pointer :: dOrbs(:, :)

    complex(R64), allocatable :: h2(:, :, :, :)
    complex(R64), allocatable :: rdm1(:, :)
    complex(R64), allocatable :: rdm2(:, :, :, :)
    complex(R64), allocatable :: invRdm1(:, :), activeProj(:, :), weakNatOrbs(:, :)

    dState(:) = 0.0_R64

    call GetTwoRdmState(twoRdmState, state)
    call GetTwoRdmState(dTwoRdmState, dState)

    call GetOrbs(orbs, state)
    call GetOrbs(dOrbs, dState)

    call TwoRdm_FillRdm1FromTwoRdmState(rdm1, twoRdmState)
    call TwoRdm_FillRdm2FromTwoRdmState(rdm2, twoRdmState)

    if (Method_Mb_OrbBased_Td2rdm_nWeakOrbs .eq. 0) then
      call Method_Mb_OrbBased_FillInverseRdm1(invRdm1, rdm1)
    else
      call FillWeakSplit(invRdm1, activeProj, weakNatOrbs, rdm1)
    end if
    call Method_Mb_OrbBased_ApplyCorrelationOp(dOrbs, orbs, invRdm1, rdm2, time, h2_=h2)

    call TwoRdm_ApplyCommutatorH2(dTwoRdmState, twoRdmState, h2)
    call Reconstruction_AddCollisionOp(dTwoRdmState, twoRdmState, h2)

    if (Method_Mb_OrbBased_Td2rdm_nWeakOrbs > 0) then
      call AddActiveWeakCommutator(dTwoRdmState, twoRdmState, orbs, weakNatOrbs, time)
      ! the last RK4 stage approximates the state at the end of the step; this
      ! keeps splitStep order2 second order (programs that modify the state call
      ! Method_Mb_OrbBased_Td2rdm_RefreshWeakOrbs)
      call StoreWeakOrbs(orbs, weakNatOrbs)
    end if

    call Orbs_ProjectOnSubspace(dOrbs, orbs)

    dState(:) = -IU * dState(:) ! out is the time-derivative!

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Active/weak split of the spatial 1-RDM (spin-summed for restricted
  !> orbitals), per body type: the truncated inverse M of the header
  !> (invRdm1), the active projector D1 M (activeProj), and the weak natural orbitals
  !> in coefficient space (weakNatOrbs, nOrbsInState x nOrbsInState), stored in
  !> the leading columns of their body type (layout of weakOrbs), zero else.
  subroutine FillWeakSplit(invRdm1, activeProj, weakNatOrbs, rdm1)
    use M_Utils_RdmDiagonalize
    use M_Orbs
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    complex(R64), intent(out), allocatable :: invRdm1(:, :)
    complex(R64), intent(out), allocatable :: activeProj(:, :)
    complex(R64), intent(out), allocatable :: weakNatOrbs(:, :)
    complex(R64), intent(in), contiguous   :: rdm1(:, :)

    real(R64), allocatable :: natOcc(:), weights(:)
    complex(R64), allocatable :: natOrbs(:, :), rdm1Spatial(:, :)
    integer(I32) :: nOS, nB, nW, ibt, iStart, iEnd, j1, k1, w
    real(R64) :: reg

    nOS = Orbs_nOrbsInState
    reg = Method_Mb_OrbBased_regularizationParameter
    nW = Method_Mb_OrbBased_Td2rdm_nWeakOrbs

    call Method_Mb_OrbBased_FillSpatialRdm1(rdm1Spatial, rdm1)

    allocate (invRdm1(nOS, nOS), activeProj(nOS, nOS), weakNatOrbs(nOS, nOS))
    invRdm1(:, :) = 0.0_R64
    activeProj(:, :) = 0.0_R64
    weakNatOrbs(:, :) = 0.0_R64

    ! the 1-RDM is block diagonal in the body types
    do ibt = 1, Method_Mb_nBodyTypes
      iStart = Method_Mb_OrbBased_nOrbsStart(ibt)
      if (iStart > nOS) exit
      iEnd = Method_Mb_OrbBased_nOrbsEnd(ibt)
      nB = iEnd - iStart + 1

      ! descending occupations: the last nW natural orbitals are weak; per
      ! block of declared m, so that the weak orbitals keep a definite m
      call RdmDiagonalize_Rdm1(natOcc, natOrbs, rdm1Spatial(iStart:iEnd, iStart:iEnd), &
                               blockOf_=Orbs_blockOfOrb(iStart:iEnd))
      weights = natOcc(:) / (natOcc(:)**2 + reg)
      weights(nB - nW + 1:nB) = 0.0_R64

      do j1 = 1, nB
        do k1 = 1, nB
          invRdm1(iStart + k1 - 1, iStart + j1 - 1) = sum(natOrbs(k1, :) * weights(:) * conjg(natOrbs(j1, :)))
          activeProj(iStart + k1 - 1, iStart + j1 - 1) = sum(natOrbs(k1, :) * natOcc(:) * weights(:) * conjg(natOrbs(j1, :)))
        end do
      end do

      do w = 1, nW
        weakNatOrbs(iStart:iEnd, iStart + w - 1) = natOrbs(:, nB - nW + w)
      end do

      deallocate (natOcc, natOrbs, weights)
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Stores the weak natural orbitals on the grid, weakOrbs = orbs weakNatOrbs,
  !> orthonormalized per body type.
  subroutine StoreWeakOrbs(orbs, weakNatOrbs)
    use M_Grid
    use M_Orbs
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    complex(R64), intent(in), contiguous :: orbs(:, :)
    complex(R64), intent(in), contiguous :: weakNatOrbs(:, :)

    integer(I32) :: nOS, ibt, iStart, nW

    nOS = Orbs_nOrbsInState
    nW = Method_Mb_OrbBased_Td2rdm_nWeakOrbs

    weakOrbs = matmul(orbs, weakNatOrbs)
    if (allocated(weakColQ)) deallocate (weakColQ)
    allocate (weakColQ(nOS))
    weakColQ(:) = .false.

    do ibt = 1, Method_Mb_nBodyTypes
      iStart = Method_Mb_OrbBased_nOrbsStart(ibt)
      if (iStart > nOS) exit

      weakColQ(iStart:iStart + nW - 1) = .true.
      call Grid_Orthonormalize(weakOrbs(:, iStart:iStart + nW - 1))
    end do

    hWeakOrbsValidQ = .false.

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> dOrbs += (H_L - h) x = (- P_w h - h P_w + 2 P_w h P_w) x for every column
  !> x of orbs (any vector, e.g. a Krylov vector), P_w the projector on the
  !> weak orbitals of the same body type. On entry dOrbs = h orbs.
  subroutine AddWeakProjectorTerms(dOrbs, orbs, time)
    use M_Grid
    use M_Orbs
    use M_Method_Mb_OrbBased

    complex(R64), intent(inout), contiguous :: dOrbs(:, :)
    complex(R64), intent(in), contiguous    :: orbs(:, :)
    real(R64), intent(in)                   :: time

    complex(R64), allocatable :: c(:), g(:), hc(:)
    logical, allocatable :: weakOfJQ(:)
    integer(I32) :: nOS, j, w, v

    nOS = Orbs_nOrbsInState

    ! Cached H weak orbitals are valid only for exactly the same time
    if (.not. hWeakOrbsValidQ .or. abs(time - hWeakOrbsTime) > 0.0_R64) then
      if (.not. allocated(hWeakOrbs)) allocate (hWeakOrbs, mold=weakOrbs)
      if (.not. allocated(hWeak)) allocate (hWeak(nOS, nOS))
      hWeakOrbs(:, :) = 0.0_R64
      call Method_Mb_OrbBased_ApplySingleBodyOp(hWeakOrbs, weakOrbs, time)
      hWeak(:, :) = 0.0_R64
      do v = 1, nOS
        do w = 1, nOS
          if (.not. (weakColQ(w) .and. weakColQ(v))) cycle
          if (Method_Mb_OrbBased_bodyTypeOfOrb(w) .ne. Method_Mb_OrbBased_bodyTypeOfOrb(v)) cycle
          hWeak(w, v) = Grid_InnerProduct(weakOrbs(:, w), hWeakOrbs(:, v))
        end do
      end do
      hWeakOrbsTime = time
      hWeakOrbsValidQ = .true.
    end if

    allocate (c(nOS), g(nOS), hc(nOS))

    do j = 1, nOS
      ! c_w = <chi_w|x_j>, g_w = <chi_w|h x_j> for the weak w of the body type of j
      weakOfJQ = weakColQ(:) .and. Method_Mb_OrbBased_bodyTypeOfOrb(1:nOS) .eq. Method_Mb_OrbBased_bodyTypeOfOrb(j)
      c(:) = 0.0_R64
      g(:) = 0.0_R64
      do w = 1, nOS
        if (.not. weakOfJQ(w)) cycle
        c(w) = Grid_InnerProduct(weakOrbs(:, w), orbs(:, j))
        g(w) = Grid_InnerProduct(weakOrbs(:, w), dOrbs(:, j))
      end do
      hc(:) = matmul(hWeak, c)

      do w = 1, nOS
        if (.not. weakOfJQ(w)) cycle
        dOrbs(:, j) = dOrbs(:, j) - weakOrbs(:, w) * g(w) - hWeakOrbs(:, w) * c(w) &
                      + 2.0_R64 * weakOrbs(:, w) * hc(w)
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> dTwoRdmState += [h1aw, D2] with the active-weak blocks of h1,
  !>   h1aw = Pi h1 (1 - Pi) + (1 - Pi) h1 Pi,  Pi = weakNatOrbs weakNatOrbs^dagger,
  !> the part of h1 that the linear step with H_L does not rotate.
  subroutine AddActiveWeakCommutator(dTwoRdmState, twoRdmState, orbs, weakNatOrbs, time)
    use M_TwoRdm
    use M_Orbs
    use M_Method_Mb_OrbBased

    complex(R64), intent(inout), contiguous :: dTwoRdmState(:)
    complex(R64), intent(in), contiguous    :: twoRdmState(:)
    complex(R64), intent(in), contiguous    :: orbs(:, :)
    complex(R64), intent(in), contiguous    :: weakNatOrbs(:, :)
    real(R64), intent(in)                   :: time

    complex(R64), allocatable :: h1(:, :), h1aw(:, :), weakProj(:, :), activeProj(:, :)
    integer(I32) :: i, nO, nOS

    nOS = Orbs_nOrbsInState
    nO = Method_Mb_OrbBased_nOrbsSum

    call Method_Mb_OrbBased_FillH1(h1, orbs, time)

    weakProj = matmul(weakNatOrbs, transpose(conjg(weakNatOrbs)))
    activeProj = -weakProj
    do i = 1, nOS
      activeProj(i, i) = activeProj(i, i) + 1.0_R64
    end do

    allocate (h1aw, mold=h1)
    h1aw(:, :) = 0.0_R64
    h1aw(1:nOS, 1:nOS) = matmul(weakProj, matmul(h1(1:nOS, 1:nOS), activeProj)) + &
                         matmul(activeProj, matmul(h1(1:nOS, 1:nOS), weakProj))
    if (nO .eq. 2 * nOS) h1aw(nOS + 1:, nOS + 1:) = h1aw(1:nOS, 1:nOS)

    call TwoRdm_ApplyCommutatorH1(dTwoRdmState, twoRdmState, h1aw)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine FillRdm1(rdm1, state)
    use M_TwoRdm

    complex(R64), intent(out), allocatable :: rdm1(:, :)
    complex(R64), intent(in), contiguous, target :: state(:)

    complex(R64), contiguous, pointer :: twoRdmState(:)

    call GetTwoRdmState(twoRdmState, state)

    call TwoRdm_FillRdm1FromTwoRdmState(rdm1, twoRdmState)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillRdm2(rdm2, state)
    use M_TwoRdm

    complex(R64), intent(out), allocatable :: rdm2(:, :, :, :)
    complex(R64), intent(in), contiguous, target :: state(:)

    complex(R64), contiguous, pointer :: twoRdmState(:)

    call GetTwoRdmState(twoRdmState, state)

    call TwoRdm_FillRdm2FromTwoRdmState(rdm2, twoRdmState)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine GetTwoRdmState(twoRdmState, state)
    use M_TwoRdm

    complex(R64), pointer, intent(out) :: twoRdmState(:)
    complex(R64), intent(in), contiguous, target :: state(:)

    twoRdmState(1:TwoRdm_nElements) => state(1:TwoRdm_nElements)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine GetOrbs(orbs, state)
    use M_Grid
    use M_Orbs
    use M_TwoRdm

    complex(R64), pointer, intent(out) :: orbs(:, :)
    complex(R64), intent(in), contiguous, target :: state(:)

    integer(I32) :: nG, nOS, nE

    nG = Grid_nPoints
    nOS = Orbs_nOrbsInState
    nE = TwoRdm_nElements

    orbs(1:nG, 1:nOS) => state(nE + 1:)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end submodule
