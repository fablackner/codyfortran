! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Absorber) S_Absorber

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Absorber_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Absorber_Linear

    call Say_Fabricate("absorber")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("absorber.linear")) then
      call Absorber_Linear_Fabricate

    else
      write (*, '(1X, A, A, A)') red, 'de: Absorber => noAbsorber', reset
      Absorber_ApplyAbsorber => NoOpProcedures_ApplyAbsorber

    end if

  end subroutine

end submodule
