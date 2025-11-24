! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_CoeffsInit) S_CoeffsInit

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine CoeffsInit_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_CoeffsInit_Load
    use M_CoeffsInit_Unary
    use M_CoeffsInit_Excited

    call Say_Fabricate("coeffsInit")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("coeffsInit.load")) then
      call CoeffsInit_Load_Fabricate

    else if (Json_GetExistence("coeffsInit.unary")) then
      call CoeffsInit_Unary_Fabricate

    else if (Json_GetExistence("coeffsInit.excited")) then
      call CoeffsInit_Excited_Fabricate

    else
      error stop "coeffsInit is missing one of: load, unary, excited"
    end if

  end subroutine

end submodule
