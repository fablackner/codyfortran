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

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end submodule
