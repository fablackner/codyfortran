! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file S_IntegratorList.f90
!> @brief Submodule implementing IntegratorList fabrication and setup.
!>
!> @details
!> This submodule provides the runtime "factory" logic that:
!> 1. Reads the JSON `"integratorList"` array.
!> 2. Dispatches to backend-specific allocators based on child names.
!> 3. Attaches the user-supplied `TimeDerivative` callback.
!> 4. Invokes each backend's `Fabricate` and later `Setup` methods.
!>
!> ## Backend Dispatch
!>
!> The child name in JSON determines which allocator is called:
!>
!> | Child contains | Allocator called                   |
!> |----------------|------------------------------------|
!> | `"rk"`         | `IntegratorList_Rk_Allocate`       |
!> | `"expokit"`    | `IntegratorList_Expokit_Allocate`  |
!> | `"gslOdeiv2"`  | `IntegratorList_GslOdeiv2_Allocate`|
!> | `"sil"`        | `IntegratorList_Sil_Allocate`      |
!> | `"cn"`         | `IntegratorList_Cn_Allocate`       |
submodule(M_IntegratorList) S_IntegratorList

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Build the integrator list from JSON configuration.
  !>
  !> Reads the `"integratorList"` array from the active JSON config, allocates
  !> one concrete backend per element, attaches the supplied `TimeDerivative`
  !> callback, and calls `Fabricate` on each.
  !>
  !> @param[in] input  Array of fabrication inputs (one per integrator).
  !>                   Each entry must supply a valid `TimeDerivative` pointer.
  !>
  !> @note The size of `input` must match the number of children in
  !>   `"integratorList"`. A mismatch triggers `error stop`.
  module subroutine IntegratorList_Fabricate(input)
    use M_Utils_Json
    use M_Utils_Say
    use M_IntegratorList_Rk, only: IntegratorList_Rk_Allocate
    use M_IntegratorList_Expokit, only: IntegratorList_Expokit_Allocate
    use M_IntegratorList_GslOdeiv2, only: IntegratorList_GslOdeiv2_Allocate
    use M_IntegratorList_Sil, only: IntegratorList_Sil_Allocate
    use M_IntegratorList_Cn, only: IntegratorList_Cn_Allocate

    type(T_IntegratorList_FabricateInput), intent(in) :: input(:)

    integer(I32) :: i, nListElements
    character(len=:), allocatable :: childName

    call Say_Fabricate("integratorList")

    !------------------------------------
    ! Bind module-level setup hook
    !------------------------------------

    IntegratorList_Setup => Setup

    !------------------------------------
    ! Allocate and fabricate each backend
    !------------------------------------

    nListElements = Json_GetNumChildren("integratorList")
    allocate (integratorList(nListElements))

    if (size(input) .ne. nListElements) then
      error stop "Input size does not match number of integratorList elements"
    end if

    do i = 1, Json_GetNumChildren("integratorList")
      if (Json_GetChildType("integratorList", i) .ne. "object") then
        error stop "Each integratorList element must be an object"
      end if

      childName = Json_GetChildName("integratorList", i)

      if (index(childName, "rk") .ne. 0) then
        call IntegratorList_Rk_Allocate(integratorList(i) % e, "integratorList."//childName)
      elseif (index(childName, "expokit") .ne. 0) then
        call IntegratorList_Expokit_Allocate(integratorList(i) % e, "integratorList."//childName)
      elseif (index(childName, "gslOdeiv2") .ne. 0) then
        call IntegratorList_GslOdeiv2_Allocate(integratorList(i) % e, "integratorList."//childName)
      elseif (index(childName, "sil") .ne. 0) then
        call IntegratorList_Sil_Allocate(integratorList(i) % e, "integratorList."//childName)
      elseif (index(childName, "cn") .ne. 0) then
        call IntegratorList_Cn_Allocate(integratorList(i) % e, "integratorList."//childName)
      else
        error stop "integratorList element must be one of: rk, expokit, gslOdeiv2, sil, cn"
      end if

      integratorList(i) % e % path = "integratorList."//childName
      integratorList(i) % e % TimeDerivative => input(i) % TimeDerivative
      call integratorList(i) % e % Fabricate()

    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Iterate over all integrators and call their `Setup` method.
  !>
  !> Invoked after fabrication, once state dimensions are finalized.
  !> Each backend allocates its internal work arrays here.
  subroutine Setup()

    integer(I32) :: i

    do i = 1, size(integratorList)
      call integratorList(i) % e % Setup()
    end do

  end subroutine

end submodule
