! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Method_Mb) S_Method_Mb

  implicit none

  !=============================================================================
  ! local data
  !=============================================================================

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Mb_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Method_Mb_OrbBased
    use M_Method_Mb_GridBased
    use M_Method_Mb_GemBased

    integer(I32) :: i, ibt, index, counter

    call Say_Fabricate("method.mb")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    Method_Mb_nBodyTypes = Json_Get("method.mb.nBodyTypes", 1)

    allocate (Method_Mb_nBodies(Method_Mb_nBodyTypes))
    Method_Mb_nBodies = Json_Get("method.mb.nBodies", [(1, i=1, Method_Mb_nBodyTypes)])

    Method_Mb_nBodiesSum = sum(Method_Mb_nBodies)

    allocate (Method_Mb_nBodiesStart(Method_Mb_nBodyTypes))
    allocate (Method_Mb_nBodiesEnd(Method_Mb_nBodyTypes))

    Method_Mb_bodyStatistics = Json_Get("method.mb.bodyStatistics", [('f', i=1, Method_Mb_nBodyTypes)])

    allocate (Method_Mb_bodyTypeOfBody(Method_Mb_nBodiesSum))

    counter = 0
    do ibt = 1, Method_Mb_nBodyTypes
      Method_Mb_nBodiesStart(ibt) = counter + 1
      do index = 1, Method_Mb_nBodies(ibt)
        counter = counter + 1
        Method_Mb_bodyTypeOfBody(counter) = ibt
      end do
      Method_Mb_nBodiesEnd(ibt) = counter
    end do

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("method.mb.orbBased")) then
      call Method_Mb_OrbBased_Fabricate

    else if (Json_GetExistence("method.mb.gridBased")) then
      call Method_Mb_GridBased_Fabricate

    else if (Json_GetExistence("method.mb.gemBased")) then
      call Method_Mb_GemBased_Fabricate

    else
      error stop "method.mb is missing one of: orbBased, gridBased, gemBased"
    end if

  end subroutine

end submodule
