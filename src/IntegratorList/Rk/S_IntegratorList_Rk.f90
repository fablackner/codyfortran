! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_IntegratorList_Rk) S_IntegratorList_Rk

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine IntegratorList_Rk_Allocate(e, path)
    use M_Utils_Json
    use M_IntegratorList_Rk_O1Expl, only: IntegratorList_Rk_O1Expl_Allocate
    use M_IntegratorList_Rk_O2Expl, only: IntegratorList_Rk_O2Expl_Allocate
    use M_IntegratorList_Rk_O4Expl, only: IntegratorList_Rk_O4Expl_Allocate
    use M_IntegratorList_Rk_O2Impl, only: IntegratorList_Rk_O2Impl_Allocate

    class(T_IntegratorList_E), allocatable, intent(out) :: e
    character(len=*), intent(in) :: path

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence(path//".o1Expl")) then
      call IntegratorList_Rk_O1Expl_Allocate(e, path//".o1Expl")

    else if (Json_GetExistence(path//".o2Expl")) then
      call IntegratorList_Rk_O2Expl_Allocate(e, path//".o2Expl")

    else if (Json_GetExistence(path//".o4Expl")) then
      call IntegratorList_Rk_O4Expl_Allocate(e, path//".o4Expl")

    else if (Json_GetExistence(path//".o2Impl")) then
      call IntegratorList_Rk_O2Impl_Allocate(e, path//".o2Impl")

    else
      error stop path//". is missing one of: o1Expl, o2Expl, o4Expl, o2Impl"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Fabricate(this)
    use M_Utils_Json
    use M_Utils_Say

    class(T_IntegratorList_E_Rk), intent(inout) :: this

    call Say_Fabricate("propagator.integrator.rk")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    call this % FabricateLevel2()

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Setup(this)
    use M_Utils_Say
    use M_Utils_UnusedVariables

    class(T_IntegratorList_E_Rk), intent(inout) :: this

    call Say_Setup("propagator.integrator.rk")

    call this % SetupLevel2()

  end subroutine

end submodule
