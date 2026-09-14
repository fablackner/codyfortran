! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Factory implementation for the Reconstruction module.
!>
!> @details
!> Reads the JSON configuration and dispatches to the appropriate family-level
!> fabricator based on which key is present:
!>
!> - `reconstruction.generic`    → Generic (spin-unresolved) schemes
!> - `reconstruction.singlet`    → Spin-singlet (S=0) schemes
!> - `reconstruction.singletTr0` → Trace-free singlet schemes
!>
!> Exactly one family must be configured; otherwise an error is raised.
submodule(M_Reconstruction) S_Reconstruction

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Reconstruction_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Reconstruction_Generic
    use M_Reconstruction_Singlet
    use M_Reconstruction_SingletTr0

    call Say_Fabricate("reconstruction")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("reconstruction.generic")) then
      call Reconstruction_Generic_Fabricate

    else if (Json_GetExistence("reconstruction.singlet")) then
      call Reconstruction_Singlet_Fabricate

    else if (Json_GetExistence("reconstruction.singletTr0")) then
      call Reconstruction_SingletTr0_Fabricate

    else
      error stop "reconstruction is missing one of: generic, singlet, singletTr0"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end submodule

