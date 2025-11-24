! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Method_Mb_OrbBased_Mctdhx) S_Method_Mb_OrbBased_Mctdhx

  implicit none

  !=============================================================================
  ! local data
  !=============================================================================

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Mb_OrbBased_Mctdhx_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Orbs
    use M_Coeffs
    use M_Method
    use M_Method_Mb_OrbBased

    call Say_Fabricate("method.manbody.orbbased.mctdhx")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    Method_Setup => Setup
    Method_Mb_OrbBased_FillRdm1 => FillRdm1
    Method_Mb_OrbBased_FillRdm2 => FillRdm2
    Method_TimeDerivative => TimeDerivative
    Method_Mb_OrbBased_TimeDerivativeOrbsLin => TimeDerivativeOrbsLin
    Method_Mb_OrbBased_TimeDerivativeCoeffsPlusOrbsNonLin => TimeDerivativeCoeffsPlusOrbsNonLin

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Setup
    use M_Utils_Say
    use M_Grid
    use M_Orbs
    use M_OrbsInit
    use M_Coeffs
    use M_CoeffsInit
    use M_Method

    integer(I32) :: nVals1, nVals2, nG, nOS, nC

    call Say_Setup("method.manbody.orbbased.mctdhx")

    nG = Grid_nPoints
    nOS = Orbs_nOrbsInState
    nC = Coeffs_nCoeffs

    nVals1 = nC
    nVals2 = nG * nOS

    allocate (Method_state(1:(nVals1 + nVals2)))

    call GetCoeffs(Coeffs_coeffs, Method_state)
    call CoeffsInit_Initialize(Coeffs_coeffs)

    call GetOrbs(Orbs_orbs, Method_state)
    call OrbsInit_Initialize(Orbs_orbs)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillRdm1(rdm1, state)
    use M_Coeffs

    complex(R64), intent(out), allocatable :: rdm1(:, :)
    complex(R64), intent(in), contiguous, target :: state(:)

    complex(R64), contiguous, pointer :: coeffs(:)

    call GetCoeffs(coeffs, state)

    call Coeffs_ApplyH1FillRdm1(coeffs, rdm1_=rdm1)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillRdm2(rdm2, state)
    use M_Coeffs

    complex(R64), intent(out), allocatable :: rdm2(:, :, :, :)
    complex(R64), intent(in), contiguous, target :: state(:)

    complex(R64), contiguous, pointer :: coeffs(:)

    call GetCoeffs(coeffs, state)

    call Coeffs_ApplyH2FillRdm2(coeffs, rdm2_=rdm2)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine TimeDerivative(dState, state, time)
    use M_Utils_Timer
    use M_Grid
    use M_Coeffs
    use M_Orbs
    use M_Method
    use M_Method_Mb_OrbBased
    use M_Utils_Constants

    complex(R64), intent(out), contiguous, target :: dState(:)
    complex(R64), intent(in), contiguous, target  :: state(:)
    real(R64), intent(in)             :: time

    complex(R64), contiguous, pointer :: coeffs(:)
    complex(R64), contiguous, pointer :: dCoeffs(:)

    complex(R64), contiguous, pointer :: orbs(:, :)
    complex(R64), contiguous, pointer :: dOrbs(:, :)

    complex(R64), allocatable :: h1(:, :)
    complex(R64), allocatable :: h2(:, :, :, :)
    complex(R64), allocatable :: rdm1(:, :)
    complex(R64), allocatable :: rdm2(:, :, :, :)

    !call Timer_Timer("rest")

    dState(:) = 0.0_R64

    call GetCoeffs(coeffs, state)
    call GetCoeffs(dCoeffs, dState)

    call GetOrbs(orbs, state)
    call GetOrbs(dOrbs, dState)

    call Coeffs_ApplyH1FillRdm1(coeffs, rdm1_=rdm1)
    call Coeffs_ApplyH2FillRdm2(coeffs, rdm2_=rdm2)

    call Method_Mb_OrbBased_ApplySingleBodyOp(dOrbs, orbs, time, h1_=h1)
    call Method_Mb_OrbBased_ApplyCorrelationOp(dOrbs, orbs, rdm1, rdm2, time, h2_=h2)

    call Coeffs_ApplyH1FillRdm1(coeffs, dCoeffs_=dCoeffs, h1_=h1)
    call Coeffs_ApplyH2FillRdm2(coeffs, dCoeffs_=dCoeffs, h2_=h2)

    call Coeffs_ProjectOnSubspace(dCoeffs, coeffs)
    call Orbs_ProjectOnSubspace(dOrbs, orbs)

    dState(:) = -IU * dState(:) ! out is the time-derivative!

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine TimeDerivativeOrbsLin(dState, state, time)
    use M_Utils_Timer
    use M_Grid
    use M_Orbs
    use M_Coeffs
    use M_Method
    use M_Method_Mb_OrbBased
    use M_Utils_Constants

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
  subroutine TimeDerivativeCoeffsPlusOrbsNonLin(dState, state, time)
    use M_Utils_Timer
    use M_Grid
    use M_Coeffs
    use M_Orbs
    use M_Method
    use M_Method_Mb_OrbBased
    use M_Utils_Constants

    complex(R64), intent(out), contiguous, target :: dState(:)
    complex(R64), intent(in), contiguous, target  :: state(:)
    real(R64), intent(in)             :: time

    complex(R64), contiguous, pointer :: coeffs(:)
    complex(R64), contiguous, pointer :: dCoeffs(:)

    complex(R64), contiguous, pointer :: orbs(:, :)
    complex(R64), contiguous, pointer :: dOrbs(:, :)

    complex(R64), allocatable :: h1(:, :)
    complex(R64), allocatable :: h2(:, :, :, :)
    complex(R64), allocatable :: rdm1(:, :)
    complex(R64), allocatable :: rdm2(:, :, :, :)

    dState(:) = 0.0_R64

    call GetCoeffs(coeffs, state)
    call GetCoeffs(dCoeffs, dState)

    call GetOrbs(orbs, state)
    call GetOrbs(dOrbs, dState)

    call Coeffs_ApplyH1FillRdm1(coeffs, rdm1_=rdm1)
    call Coeffs_ApplyH2FillRdm2(coeffs, rdm2_=rdm2)

    call Method_Mb_OrbBased_FillH1(h1, orbs, time)
    call Method_Mb_OrbBased_ApplyCorrelationOp(dOrbs, orbs, rdm1, rdm2, time, h2_=h2)

    call Coeffs_ApplyH1FillRdm1(coeffs, dCoeffs_=dCoeffs, h1_=h1)
    call Coeffs_ApplyH2FillRdm2(coeffs, dCoeffs_=dCoeffs, h2_=h2)

    call Coeffs_ProjectOnSubspace(dCoeffs, coeffs)
    call Orbs_ProjectOnSubspace(dOrbs, orbs)

    dState(:) = -IU * dState(:) ! out is the time-derivative!

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine GetCoeffs(coeffs, state)
    use M_Coeffs
    complex(R64), pointer, intent(out) :: coeffs(:)
    complex(R64), intent(in), contiguous, target :: state(:)

    coeffs(1:Coeffs_nCoeffs) => state(1:Coeffs_nCoeffs)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine GetOrbs(orbs, state)
    use M_Coeffs
    use M_Grid
    use M_Orbs
    use M_Method_Mb_OrbBased
    complex(R64), pointer, intent(out) :: orbs(:, :)
    complex(R64), intent(in), contiguous, target :: state(:)

    orbs(1:Grid_nPoints, 1:Orbs_nOrbsInState) => state(Coeffs_nCoeffs + 1:)

  end subroutine

end submodule
