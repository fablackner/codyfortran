! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_ConfigList) S_ConfigList

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine ConfigList_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_ConfigList_AllActive, only: ConfigList_AllActive_Allocate

    integer(I32) :: i, nListElements
    character(len=:), allocatable :: childName

    call Say_Fabricate("configList")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    ConfigList_Setup => Setup

    !------------------------------------
    ! branch
    !------------------------------------

    nListElements = Json_GetNumChildren("configList")
    allocate (configList(nListElements))

    do i = 1, Json_GetNumChildren("configList")
      if (Json_GetChildType("configList", i) .ne. "object") then
        error stop "Each configList element must be an object"
      end if

      childName = Json_GetChildName("configList", i)

      if (index(childName, "allActive") .ne. 0) then
        call ConfigList_AllActive_Allocate(configList(i) % e, "configList."//childName)
      else
        error stop "configList element must be one of: allActive"
      end if

      configList(i) % e % path = "configList."//childName
      call configList(i) % e % Fabricate()

    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Setup()

    integer(I32) :: i

    do i = 1, size(configList)
      call configList(i) % e % Setup()
    end do

  end subroutine

end submodule
