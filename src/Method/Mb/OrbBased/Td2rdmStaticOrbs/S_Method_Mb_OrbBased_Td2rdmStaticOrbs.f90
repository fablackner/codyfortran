! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Method_Mb_OrbBased_Td2rdmStaticOrbs) S_Method_Mb_OrbBased_Td2rdmStaticOrbs

  implicit none

  !=============================================================================
  ! local data
  !=============================================================================

  complex(R64), allocatable, target                  :: staticComponents(:)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Mb_OrbBased_Td2rdmStaticOrbs_Fabricate
    use M_Utils_Say
    use M_Method
    use M_Method_Mb_OrbBased

    call Say_Fabricate("method.mb.orbBased.td2rdmStaticOrbs")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    Method_Setup => Setup
    Method_Mb_OrbBased_staticOrbsQ = .true.
    Method_Mb_OrbBased_implementsSwitchingQ = .true.
    Method_Mb_OrbBased_FillRdm1 => FillRdm1
    Method_Mb_OrbBased_FillRdm2 => FillRdm2
    Method_TimeDerivative => TimeDerivative
    Method_GetEnergy => GetEnergy

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Setup
    use M_Utils_Say
    use M_Grid
    use M_Orbs
    use M_TwoRdm
    use M_TwoRdmInit
    use M_Orbs
    use M_OrbsInit
    use M_Method
    use M_Method_Mb_OrbBased

    integer(I32) :: nVals1, nVals2, nG, nOS, nE

    nG = Grid_nPoints
    nOS = Orbs_nOrbsInState
    nE = TwoRdm_nElements

    nVals1 = nE
    nVals2 = nG * nOS

    call Say_Setup("method.mb.orbBased.td2rdmStaticOrbs")
    call Method_Mb_OrbBased_Setup

    allocate (Method_state(1:nVals1))
    allocate (staticComponents(1:nVals2))

    TwoRdm_twoRdmState(1:nE) => Method_state(1:nVals1)
    call TwoRdmInit_Initialize(TwoRdm_twoRdmState)

    Orbs_orbs(1:nG, 1:nOS) => staticComponents(1:nVals2)
    call OrbsInit_Initialize(Orbs_orbs)

    ! the static Hamiltonian; the time-dependent terms follow in
    ! Method_Mb_OrbBased_RefillStaticH1
    call Method_Mb_OrbBased_FillH0(Method_Mb_OrbBased_h1, Orbs_orbs, 0.0_R64)
    call Method_Mb_OrbBased_FillH2(Method_Mb_OrbBased_h2, Orbs_orbs, 0.0_R64)
    call Method_Mb_OrbBased_SetSwitchingReference

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine TimeDerivative(dTwoRdmState, twoRdmState, time)
    use M_Utils_Constants
    use M_Reconstruction
    use M_TwoRdm
    use M_Orbs
    use M_Method
    use M_Method_Mb_OrbBased

    complex(R64), intent(out), contiguous, target :: dTwoRdmState(:)
    complex(R64), intent(in), contiguous, target  :: twoRdmState(:)
    real(R64), intent(in)             :: time

    complex(R64), allocatable :: h2Switched(:, :, :, :)
    real(R64) :: lambda

    ! h1 at this time if any of its terms depends on time (as TDCI); h2 is
    ! static and enters as lambda h2 during the adiabatic switching. The
    ! reconstruction gets the Hamiltonian it closes (a closure may depend on
    ! it), so the collision operator carries lambda through its h2.
    call Method_Mb_OrbBased_RefillStaticH1(Method_Mb_OrbBased_h1, Orbs_orbs, time)
    lambda = Method_Mb_OrbBased_SwitchFactor(time)

    dTwoRdmState(:) = 0.0_R64

    if (lambda < 1.0_R64) then
      h2Switched = lambda * Method_Mb_OrbBased_h2
      call AddTwoBodyTerms(h2Switched)
    else
      call AddTwoBodyTerms(Method_Mb_OrbBased_h2)
    end if
    call TwoRdm_ApplyCommutatorH1(dTwoRdmState, twoRdmState, Method_Mb_OrbBased_h1)

    dTwoRdmState(:) = -IU * dTwoRdmState(:) ! out is the time-derivative!

  contains

    subroutine AddTwoBodyTerms(h2)
      complex(R64), intent(in), contiguous :: h2(:, :, :, :)

      call TwoRdm_ApplyCommutatorH2(dTwoRdmState, twoRdmState, h2)
      call Reconstruction_SetHamiltonian(Method_Mb_OrbBased_h1, h2)
      call Reconstruction_AddCollisionOp(dTwoRdmState, twoRdmState, h2)

    end subroutine

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Energy of the Hamiltonian at time with the stored matrices (the static
  !> Hamiltonian is not rebuilt from the grid: a program may have replaced
  !> it): h1 at time, lambda h2; fills Method_Mb_OrbBased_rdm1 and _rdm2 as
  !> the generic energy does.
  function GetEnergy(time) result(res)
    use M_Utils_RdmObservables
    use M_Orbs
    use M_Method
    use M_Method_Mb_OrbBased

    real(R64), intent(in) :: time
    real(R64) :: res

    call Method_Mb_OrbBased_RefillStaticH1(Method_Mb_OrbBased_h1, Orbs_orbs, time)
    call FillRdm1(Method_Mb_OrbBased_rdm1, Method_state)
    call FillRdm2(Method_Mb_OrbBased_rdm2, Method_state)
    res = RdmObservables_Energy(Method_Mb_OrbBased_rdm1, Method_Mb_OrbBased_h1, &
                                Method_Mb_OrbBased_SwitchFactor(time) * Method_Mb_OrbBased_rdm2, Method_Mb_OrbBased_h2)

  end function

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
    use M_Utils_UnusedVariables

    complex(R64), pointer, intent(out) :: orbs(:, :)
    complex(R64), intent(in), contiguous, target :: state(:)

    integer(I32) :: nG, nOS

    if (.false.) call UnusedVariables_Mark(state)

    if (.not. allocated(staticComponents)) then
      nullify (orbs)
      return
    end if

    nG = Grid_nPoints
    nOS = Orbs_nOrbsInState

    orbs(1:nG, 1:nOS) => staticComponents(1:nG * nOS)

  end subroutine

end submodule
