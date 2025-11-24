! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Propagator) S_Propagator

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Propagator_Fabricate()
    use M_Utils_Json
    use M_Utils_Say
    use M_Propagator_Single
    use M_Propagator_SplitStep
    use M_Propagator_EigenExpansion

    call Say_Fabricate("propagator")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("propagator.single")) then
      call Propagator_Single_Fabricate

    else if (Json_GetExistence("propagator.splitStep")) then
      call Propagator_SplitStep_Fabricate

    else if (Json_GetExistence("propagator.eigenExpansion")) then
      call Propagator_EigenExpansion_Fabricate

    else
      error stop "propagator is missing one of: single, splitStep"
    end if

  end subroutine

end submodule
