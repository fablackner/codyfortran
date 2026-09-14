! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Reconstruction_Generic) S_Reconstruction_Generic

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Reconstruction_Generic_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Reconstruction_Generic_Valdemoro

    call Say_Fabricate("reconstruction.generic")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("reconstruction.generic.genericValdemoro")) then
      call Reconstruction_Generic_Valdemoro_Fabricate

    else
      error stop "reconstruction.generic is missing one of: genericValdemoro"
    end if

  end subroutine

end submodule

