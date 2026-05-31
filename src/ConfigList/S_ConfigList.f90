! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file S_ConfigList.f90
!> @brief Implementation of ConfigList fabrication and setup orchestration.
!>
!> @details
!> This submodule implements the factory logic for parsing JSON configuration
!> and instantiating the appropriate concrete configuration elements. It acts
!> as the top-level dispatcher, delegating to truncation-specific allocators.
!>
!> ## JSON Structure Expected
!>
!> ```json
!> {
!>   "configList": {
!>     "allActive1": { "fermionic": { "bodyTarget": 1, "nExcitations": 2 } },
!>     "allActive2": { "bosonic":   { "bodyTarget": 2, "nExcitations": 3 } }
!>   }
!> }
!> ```
!>
!> ## Extension Point
!>
!> To add new truncation schemes (e.g., "restricted", "generalized"):
!> 1. Create `M_ConfigList_Restricted` with `ConfigList_Restricted_Allocate`
!> 2. Add `else if (index(childName, "restricted") .ne. 0)` branch below
submodule(M_ConfigList) S_ConfigList

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Factory routine: parse JSON and instantiate configuration elements.
!>
!> @details
!> Iterates over children of "configList" in the JSON input, identifies the
!> truncation scheme from the key name, and delegates to the appropriate
!> allocator. Each allocated element is then fabricated (reads its own params).
!>
!> ## Supported Truncation Schemes
!>
!> | Key contains  | Allocator                         | Description          |
!> |---------------|-----------------------------------|----------------------|
!> | "allActive"   | ConfigList_AllActive_Allocate     | Full active space    |
!>
!> @post `configList(:)` allocated with one element per JSON child
!> @post Each element's `path` and `Fabricate()` have been called
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
    ! branch on truncation scheme
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
!> @brief Global setup: iterate all elements and call their Setup methods.
!>
!> @details
!> After fabrication, each element knows its nConfigurations but hasn't built
!> the codeFromConfig mapping or singles/doubles connectivity. This routine
!> triggers that second-phase initialization for all elements.
!>
!> Called via the `ConfigList_Setup` procedure pointer set during fabrication.
  subroutine Setup()

    integer(I32) :: i

    do i = 1, size(configList)
      call configList(i) % e % Setup()
    end do

  end subroutine

end submodule
