! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Dispatch submodule for mixing backend selection.
!>
!> @details
!> Reads the JSON configuration under `mixing.*` and dispatches to the
!> appropriate backend. Currently supported:
!>
!>   - `mixing.linear`: Plain damped linear mixing (default if no block given)
!>   - `mixing.diis`: Pulay/Anderson (DIIS) mixing with residual history
!>
!> This is the only file that needs modification when adding new mixing
!> backends.
submodule(M_Mixing) S_Mixing

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Mixing_Fabricate()
    use M_Utils_Json
    use M_Utils_Say
    use M_Mixing_Linear
    use M_Mixing_Diis

    call Say_Fabricate("mixing")

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("mixing.diis")) then
      call Mixing_Diis_Fabricate

    else if (Json_GetExistence("mixing.linear")) then
      call Mixing_Linear_Fabricate

    else if (Json_GetExistence("mixing")) then
      error stop "mixing is missing one of: linear, diis"

    else
      ! linear is the default when no mixing block is configured
      call Mixing_Linear_Fabricate
    end if

  end subroutine

end submodule
