! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Purification_Generic) S_Purification_Generic

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Purification_Generic_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Purification_Generic_Alcoba

    call Say_Fabricate("purification.generic")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("purification.generic.alcoba")) then
      call Purification_Generic_Alcoba_Fabricate

    else
      error stop "purification.generic is missing one of: alcoba"
    end if

  end subroutine

end submodule

