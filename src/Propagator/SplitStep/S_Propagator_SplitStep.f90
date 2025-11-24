! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Propagator_SplitStep) S_Propagator_SplitStep

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Propagator_SplitStep_Fabricate
    use M_Utils_Say
    use M_Utils_Json
    use M_Propagator_SplitStep_Order2
    use M_Propagator_SplitStep_Order4

    call Say_Fabricate("propagator.splitStep")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("propagator.splitStep.order2")) then
      call Propagator_SplitStep_Order2_Fabricate

    else if (Json_GetExistence("propagator.splitStep.order4")) then
      call Propagator_SplitStep_Order4_Fabricate

    else
      error stop "propagator.splitStep is missing one of: order2, order4"
    end if

  end subroutine

end submodule
