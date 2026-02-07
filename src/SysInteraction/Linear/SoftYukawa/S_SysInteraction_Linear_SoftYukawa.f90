! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysInteraction_Linear_SoftYukawa) S_SysInteraction_Linear_SoftYukawa

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysInteraction_Linear_SoftYukawa_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysInteraction
    use M_SysInteraction_Linear_SoftYukawa_StdImpl
    use M_SysInteraction_Linear_SoftYukawa_Fftw

    call Say_Fabricate("sysInteraction.linear.softYukawa")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysInteraction_Linear_SoftYukawa_strength = Json_Get("sysInteraction.linear.softYukawa.strength", 1.0_R64)
    SysInteraction_Linear_SoftYukawa_softening1 = Json_Get("sysInteraction.linear.softYukawa.softening1", 1.0_R64)
    SysInteraction_Linear_SoftYukawa_softening2 = Json_Get("sysInteraction.linear.softYukawa.softening2", 0.0_R64)
    SysInteraction_Linear_SoftYukawa_dampening = Json_Get("sysInteraction.linear.softYukawa.dampening", 0.0_R64)

    SysInteraction_timeIndependentQ = .true.
    SysInteraction_bodyTypeIndependentQ = .true.

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("sysInteraction.linear.softYukawa.stdImpl")) then
      call SysInteraction_Linear_SoftYukawa_StdImpl_Fabricate

    else if (Json_GetExistence("sysInteraction.linear.softYukawa.fftw")) then
      call SysInteraction_Linear_SoftYukawa_Fftw_Fabricate

    else
      error stop "sysInteraction.linear.softYukawa is missing one of: stdImpl, fftw"
    end if

  end subroutine

  !--------------------------------------------------------------------
  pure module function SysInteraction_Linear_SoftYukawa_Interaction(distance) result(res)
    real(R64) :: res
    real(R64), intent(in) :: distance

    res = SysInteraction_Linear_SoftYukawa_strength * exp(-distance * SysInteraction_Linear_SoftYukawa_dampening) &
          / (sqrt(distance**2 + SysInteraction_Linear_SoftYukawa_softening1**2) + SysInteraction_Linear_SoftYukawa_softening2)
  end function

end submodule
