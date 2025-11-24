! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Method_Mb_OrbBased_Tdhx) S_Method_Mb_OrbBased_Tdhx

  implicit none

  !=============================================================================
  ! local data
  !=============================================================================

  complex(R64), allocatable, target                  :: staticComponents(:)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Mb_OrbBased_Tdhx_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Orbs
    use M_Coeffs
    use M_Method
    use M_Method_Mb_OrbBased

    call Say_Fabricate("method.mb.orbbased.tdhx")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    Method_Setup => Setup
    Method_Mb_OrbBased_FillRdm1 => FillRdm1
    Method_Mb_OrbBased_FillRdm2 => FillRdm2
    Method_TimeDerivative => TimeDerivative

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
    use M_Method_Mb_OrbBased

    integer(I32) :: nVals1, nVals2, nG, nOS

    call Say_Setup("method.mb.orbbased.tdhx")

    nG = Grid_nPoints
    nOS = Orbs_nOrbsInState

    nVals1 = Coeffs_nCoeffs
    nVals2 = nG * nOS

    allocate (staticComponents(1:nVals1))
    allocate (Method_state(1:nVals2))

    Coeffs_coeffs(1:nVals1) => staticComponents(1:nVals1)
    call CoeffsInit_Initialize(Coeffs_coeffs)

    Orbs_orbs(1:nG, 1:nOS) => Method_state(1:nVals2)
    call OrbsInit_Initialize(Orbs_orbs)

    call Coeffs_ApplyH1FillRdm1(Coeffs_coeffs, rdm1_=Method_Mb_OrbBased_rdm1)
    call Coeffs_ApplyH2FillRdm2(Coeffs_coeffs, rdm2_=Method_Mb_OrbBased_rdm2)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine TimeDerivative(dState, state, time)
    use M_Utils_Constants
    use M_Utils_Timer
    use M_Grid
    use M_Orbs
    use M_Method
    use M_Method_Mb_OrbBased

    complex(R64), intent(out), contiguous, target :: dState(:)
    complex(R64), intent(in), contiguous, target  :: state(:)
    real(R64), intent(in)             :: time

    complex(R64), contiguous, pointer :: orbs(:, :)
    complex(R64), contiguous, pointer :: dOrbs(:, :)

    integer(I32) :: nG, nOS

    !call Timer_Timer("rest")

    dState(:) = 0.0_R64

    nG = Grid_nPoints
    nOS = Orbs_nOrbsInState

    orbs(1:nG, 1:nOS) => state(1:)
    dOrbs(1:nG, 1:nOS) => dState(1:)

    call Method_Mb_OrbBased_ApplySingleBodyOp(dOrbs, orbs, time)
    call Method_Mb_OrbBased_ApplyHartreeFockOp(dOrbs, orbs, Method_Mb_OrbBased_rdm1, time)

    dState(:) = -IU * dState(:) ! out is the time-derivative!

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
  subroutine GetCoeffs(coeffs, state)
    use M_Coeffs
    use M_Utils_UnusedVariables

    complex(R64), pointer, intent(out) :: coeffs(:)
    complex(R64), intent(in), contiguous, target :: state(:)

    if (.false.) call UnusedVariables_Mark(state)

    if (.not. allocated(staticComponents)) then
      nullify (coeffs)
      return
    end if

    coeffs(1:Coeffs_nCoeffs) => staticComponents(1:Coeffs_nCoeffs)

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
