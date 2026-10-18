! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Method_Mb_OrbBased_Td2rdm) S_Method_Mb_OrbBased_Td2rdm

  implicit none

  !=============================================================================
  ! local data
  !=============================================================================

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
    ! set values and procedure pointers
    !------------------------------------

    Method_Setup => Setup
    Method_Mb_OrbBased_FillRdm1 => FillRdm1
    Method_Mb_OrbBased_FillRdm2 => FillRdm2
    Method_TimeDerivative => TimeDerivative
    Method_Mb_OrbBased_TimeDerivativeOrbsLin => TimeDerivativeOrbsLin
    Method_Mb_OrbBased_TimeDerivativeCoeffsPlusOrbsNonLin => TimeDerivativeTwoRdmPlusOrbsNonLin

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

    integer(I32) :: nVals1, nVals2, nG, nOS, nE

    nG = Grid_nPoints
    nOS = Orbs_nOrbsInState
    nE = TwoRdm_nElements

    nVals1 = nE
    nVals2 = nG * nOS

    call Say_Setup("method.mb.orbbased.td2rdm")

    allocate (Method_state(1:(nVals1 + nVals2)))

    TwoRdm_twoRdmState(1:nE) => Method_state(1:nVals1)
    call TwoRdmInit_Initialize(TwoRdm_twoRdmState)

    Orbs_orbs(1:nG, 1:nOS) => Method_state(nVals1 + 1:(nVals1 + nVals2))
    call OrbsInit_Initialize(Orbs_orbs)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine TimeDerivative(dState, state, time)
    use M_Utils_Constants
    use M_Grid
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

    integer(I32) :: nG, nOS, nE

    dState(:) = 0.0_R64

    nG = Grid_nPoints
    nOS = Orbs_nOrbsInState
    nE = TwoRdm_nElements

    twoRdmState(1:nE) => state(1:)
    orbs(1:nG, 1:nOS) => state(nE + 1:)

    dTwoRdmState(1:nE) => dState(1:)
    dOrbs(1:nG, 1:nOS) => dState(nE + 1:)

    call TwoRdm_FillRdm1FromTwoRdmState(rdm1, twoRdmState)
    call TwoRdm_FillRdm2FromTwoRdmState(rdm2, twoRdmState)

    call Method_Mb_OrbBased_ApplySingleBodyOp(dOrbs, orbs, time, h1_=h1)
    call Method_Mb_OrbBased_ApplyCorrelationOp(dOrbs, orbs, rdm1, rdm2, time, h2_=h2)

    call TwoRdm_ApplyCommutatorH1(dTwoRdmState, twoRdmState, h1)
    call TwoRdm_ApplyCommutatorH2(dTwoRdmState, twoRdmState, h2)
    call Reconstruction_AddCollisionOp(dTwoRdmState, twoRdmState, h2)

    call Orbs_ProjectOnSubspace(dOrbs, orbs)

    dState(:) = -IU * dState(:) ! out is the time-derivative!

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Linear split-step part: the one-body operator on the orbitals only (the
  !> 2-RDM is constant), dOrbs = -i h orbs without projection.
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

    call Method_Mb_OrbBased_ApplySingleBodyOp(dOrbs, orbs, time)

    dState(:) = -IU * dState(:) ! out is the time-derivative!

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Non-linear split-step part: the two-body interaction only. The 2-RDM
  !> evolves with [h2, D2] plus the collision operator, the orbitals with the
  !> projected correlation term. The one-body commutator is absent because the
  !> linear part rotates the orbitals with the full (unprojected) h.
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

    dState(:) = 0.0_R64

    call GetTwoRdmState(twoRdmState, state)
    call GetTwoRdmState(dTwoRdmState, dState)

    call GetOrbs(orbs, state)
    call GetOrbs(dOrbs, dState)

    call TwoRdm_FillRdm1FromTwoRdmState(rdm1, twoRdmState)
    call TwoRdm_FillRdm2FromTwoRdmState(rdm2, twoRdmState)

    call Method_Mb_OrbBased_ApplyCorrelationOp(dOrbs, orbs, rdm1, rdm2, time, h2_=h2)

    call TwoRdm_ApplyCommutatorH2(dTwoRdmState, twoRdmState, h2)
    call Reconstruction_AddCollisionOp(dTwoRdmState, twoRdmState, h2)

    call Orbs_ProjectOnSubspace(dOrbs, orbs)

    dState(:) = -IU * dState(:) ! out is the time-derivative!

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
