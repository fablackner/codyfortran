! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Method_Mb_OrbBased_Tdci) S_Method_Mb_OrbBased_Tdci

  implicit none

  !=============================================================================
  ! local data
  !=============================================================================

  complex(R64), allocatable, target                  :: staticComponents(:)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Mb_OrbBased_Tdci_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Orbs
    use M_Coeffs
    use M_Method
    use M_Method_Mb_OrbBased

    call Say_Fabricate("method.mb.orbbased.tdci")

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
    use M_Coeffs
    use M_Orbs
    use M_Method
    use M_Method_Mb_OrbBased
    use M_CoeffsInit
    use M_OrbsInit

    integer(I32) :: nVals1, nVals2, nG, nOS

    call Say_Setup("method.mb.orbbased.tdci")

    nG = Grid_nPoints
    nOS = Orbs_nOrbsInState

    nVals1 = Coeffs_nCoeffs
    nVals2 = nG * nOS

    allocate (Method_state(1:nVals1))
    allocate (staticComponents(1:nVals2))

    Coeffs_coeffs(1:nVals1) => Method_state(1:nVals1)
    call CoeffsInit_Initialize(Coeffs_coeffs)

    Orbs_orbs(1:nG, 1:nOS) => staticComponents(1:nVals2)
    call OrbsInit_Initialize(Orbs_orbs)

    call Method_Mb_OrbBased_FillH1(Method_Mb_OrbBased_h1, Orbs_orbs, 0.0_R64)
    call Method_Mb_OrbBased_FillH2(Method_Mb_OrbBased_h2, Orbs_orbs, 0.0_R64)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine TimeDerivative(dCoeffs, coeffs, time)
    use M_Utils_Timer
    use M_Utils_UnusedVariables
    use M_Utils_Constants
    use M_Coeffs
    use M_Method_Mb_OrbBased
    use M_SysKinetic
    use M_SysPotential
    use M_SysInteraction
    use M_Orbs

    complex(R64), intent(out), contiguous, target :: dCoeffs(:)
    complex(R64), intent(in), contiguous, target  :: coeffs(:)
    real(R64), intent(in)             :: time

    if (.not. SysPotential_timeIndependentQ .or. .not. SysKinetic_timeIndependentQ) then
      call Method_Mb_OrbBased_FillH1(Method_Mb_OrbBased_h1, Orbs_orbs, time)
    end if

    if (.not. SysInteraction_timeIndependentQ) then
      call Method_Mb_OrbBased_FillH2(Method_Mb_OrbBased_h2, Orbs_orbs, time)
    end if

    dCoeffs(:) = 0.0_R64
    !call Timer_Timer("start")
    call Coeffs_ApplyH1FillRdm1(coeffs, dCoeffs_=dCoeffs, h1_=Method_Mb_OrbBased_h1)
    !call Timer_Timer("h1")
    call Coeffs_ApplyH2FillRdm2(coeffs, dCoeffs_=dCoeffs, h2_=Method_Mb_OrbBased_h2)
    !call Timer_Timer("h2")

    dCoeffs(:) = -IU * dCoeffs(:) ! out is the time-derivative!

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

    complex(R64), pointer, intent(out) :: coeffs(:)
    complex(R64), intent(in), contiguous, target :: state(:)

    coeffs(1:Coeffs_nCoeffs) => state(1:Coeffs_nCoeffs)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine GetOrbs(orbs, state)
    use M_Grid
    use M_Orbs
    use M_Utils_UnusedVariables

    complex(R64), pointer, intent(out) :: orbs(:, :)
    complex(R64), intent(in), contiguous, target :: state(:)

    if (.false.) call UnusedVariables_Mark(state)

    if (.not. associated(Orbs_orbs)) then
      nullify (orbs)
      return
    end if

    orbs => Orbs_orbs

  end subroutine

end submodule
