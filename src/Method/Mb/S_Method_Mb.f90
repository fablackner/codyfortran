! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Method_Mb) S_Method_Mb
  !-----------------------------------------------------------------------------
  ! Many-body (Mb) method dispatcher and metadata initialization.
  !
  ! This submodule:
  !   1. Parses body-type configuration (species, particle counts, statistics)
  !   2. Builds index mappings for packed storage layouts
  !   3. Dispatches to the selected representation: OrbBased, GridBased, or GemBased
  !-----------------------------------------------------------------------------

  implicit none

  !=============================================================================
  ! local data
  !=============================================================================

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Method_Mb_Fabricate
    !---------------------------------------------------------------------------
    ! Initializes many-body metadata from JSON and branches to representation.
    !
    ! JSON keys read:
    !   - method.mb.nBodyTypes:     Number of distinct particle species/spins
    !   - method.mb.nBodies:        Array of particle counts per body type
    !   - method.mb.bodyStatistics: Array of 'f' (fermion) or 'b' (boson)
    !
    ! Computed arrays:
    !   - bodyTypeOfBody(1:nBodiesSum): Maps global body index → body type
    !   - nBodiesStart/End(bt):         Index range for body type bt
    !---------------------------------------------------------------------------
    use M_Utils_Json
    use M_Utils_Say
    use M_Method_Mb_OrbBased
    use M_Method_Mb_GridBased
    use M_Method_Mb_GemBased

    integer(I32) :: i, ibt, index, counter

    call Say_Fabricate("method.mb")

    !------------------------------------
    ! parse body-type configuration
    !------------------------------------

    Method_Mb_nBodyTypes = Json_Get("method.mb.nBodyTypes", 1)

    allocate (Method_Mb_nBodies(Method_Mb_nBodyTypes))
    Method_Mb_nBodies = Json_Get("method.mb.nBodies", [(1, i=1, Method_Mb_nBodyTypes)])

    Method_Mb_nBodiesSum = sum(Method_Mb_nBodies)

    allocate (Method_Mb_nBodiesStart(Method_Mb_nBodyTypes))
    allocate (Method_Mb_nBodiesEnd(Method_Mb_nBodyTypes))

    Method_Mb_bodyStatistics = Json_Get("method.mb.bodyStatistics", [('f', i=1, Method_Mb_nBodyTypes)])

    !------------------------------------
    ! build body-type index mappings
    !------------------------------------

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
    ! branch: representation selection
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
