! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_OrbsInit) S_OrbsInit

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine OrbsInit_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_OrbsInit_Linear
    use M_OrbsInit_GridPoint
    use M_OrbsInit_Lattice
    use M_OrbsInit_Load
    use M_OrbsInit_Ylm

    call Say_Fabricate("orbsInit")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    OrbsInit_Initialize => Initialize

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("orbsInit.linear")) then
      call OrbsInit_Linear_Fabricate

    else if (Json_GetExistence("orbsInit.lattice")) then
      call OrbsInit_Lattice_Fabricate

    else if (Json_GetExistence("orbsInit.gridPoint")) then
      call OrbsInit_GridPoint_Fabricate

    else if (Json_GetExistence("orbsInit.load")) then
      call OrbsInit_Load_Fabricate

    else if (Json_GetExistence("orbsInit.ylm")) then
      call OrbsInit_Ylm_Fabricate

    else
      error stop "orbsInit is missing one of: linear, load, ylm"
    end if

  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Initialize(orbs)
    use M_Utils_DataStorage
    use M_Method_Mb
    use M_Method_Mb_OrbBased
    use M_Grid

    complex(R64), intent(out), contiguous  :: orbs(:, :)

    integer(I32) :: ind, ibt, i1

    i1 = 0
    do ibt = 1, Method_Mb_nBodyTypes
      do ind = 1, Method_Mb_OrbBased_nOrbs(ibt)
        i1 = i1 + 1

        call OrbsInit_InitializeOrb(orbs(:, i1), ind, ibt)

      end do
    end do

    !call Orbs_Orthonormalize(orbs)

  end subroutine

end submodule
