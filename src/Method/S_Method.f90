! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Method) S_Method
  !-----------------------------------------------------------------------------
  ! Dispatcher for method selection based on JSON configuration.
  !
  ! The fabricate routine inspects the JSON tree and delegates to either:
  !   - Sb (single-body): external potential only, no interactions
  !   - Mb (many-body):   multi-particle with body types and statistics
  !
  ! Also the two views of Method_TimeDerivative that diagonalizers and
  ! imaginary-time integrators take as callbacks: Method_ApplyHamiltonian
  ! and Method_ImagTimeDerivative.
  !-----------------------------------------------------------------------------

  implicit none

  !=============================================================================
  ! local data
  !=============================================================================

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Fabricate
    !---------------------------------------------------------------------------
    ! Entry point for method initialization.
    !
    ! Reads the top-level "method" key from JSON and branches to the appropriate
    ! single-body or many-body fabrication routine. Exactly one of "method.sb"
    ! or "method.mb" must be present in the configuration.
    !---------------------------------------------------------------------------
    use M_Utils_Json
    use M_Utils_Say
    use M_Method_Sb
    use M_Method_Mb

    call Say_Fabricate("method")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    ! Optional entry points: null unless the selected method binds them (no
    ! binding of an earlier fabricate or of a test survives)
    nullify (Method_PropagateExactPart, Method_ApplyExactPartPhi)

    ! Set by a Kohn-Sham method only
    Method_kohnShamQ = .false.

    !------------------------------------
    ! branch: single-body vs many-body
    !------------------------------------

    if (Json_GetExistence("method.sb")) then
      call Method_Sb_Fabricate

    else if (Json_GetExistence("method.mb")) then
      call Method_Mb_Fabricate

    else
      error stop "method is missing one of: sb, mb"
    end if

    ! A density functional (sysXc) replaces the exchange and correlation of
    ! the interaction: only a Kohn-Sham method uses it (with a correlated
    ! wavefunction it would count correlation twice)
    if (Json_GetExistence("sysXc") .and. .not. Method_kohnShamQ) &
      error stop "sysXc (an exchange-correlation functional) is only used by a Kohn-Sham method (method.mb.orbBased.tdks)"

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_ApplyHamiltonian(dState, state, time)
    use M_Utils_Constants

    complex(R64), intent(out), contiguous, target :: dState(:)
    complex(R64), intent(in), contiguous, target  :: state(:)
    real(R64), intent(in) :: time

    call Method_TimeDerivative(dState, state, time)
    dState(:) = IU * dState(:)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_ImagTimeDerivative(dState, state, time)
    use M_Utils_Constants

    complex(R64), intent(out), contiguous, target :: dState(:)
    complex(R64), intent(in), contiguous, target  :: state(:)
    real(R64), intent(in) :: time

    call Method_TimeDerivative(dState, state, time)
    dState(:) = -IU * dState(:)

  end subroutine

end submodule
