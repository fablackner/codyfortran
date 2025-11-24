! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Absorber_Linear) S_Absorber_Linear

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Absorber_Linear_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Absorber_Linear_Cosinus

    call Say_Fabricate("absorber.linear")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("absorber.linear.cosinus")) then
      call Absorber_Linear_Cosinus_Fabricate

    else
      error stop "absorber.linear is missing one of: cosinus"
    end if

  end subroutine

end submodule
