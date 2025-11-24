! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_DiagonalizerList) S_DiagonalizerList

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine DiagonalizerList_Fabricate(input)
    use M_Utils_Json
    use M_Utils_Say
    use M_DiagonalizerList_Lapack, only: DiagonalizerList_Lapack_Allocate
    use M_DiagonalizerList_Arpack, only: DiagonalizerList_Arpack_Allocate

    type(T_DiagonalizerList_FabricateInput), intent(in) :: input(:)

    integer(I32) :: i, nListElements
    character(len=:), allocatable :: childName

    call Say_Fabricate("diagonalizerList")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    DiagonalizerList_Setup => Setup

    !------------------------------------
    ! branch
    !------------------------------------

    nListElements = Json_GetNumChildren("diagonalizerList")
    allocate (DiagonalizerList(nListElements))

    if (size(input) .ne. nListElements) then
      error stop "Input size does not match number of diagonalizerList elements"
    end if

    do i = 1, nListElements
      if (Json_GetChildType("diagonalizerList", i) .ne. "object") then
        error stop "Each diagonalizerList element must be an object"
      end if

      childName = Json_GetChildName("diagonalizerList", i)

      if (index(childName, "lapack") .ne. 0) then
        call DiagonalizerList_Lapack_Allocate(DiagonalizerList(i) % e, "diagonalizerList."//childName)
      else if (index(childName, "arpack") .ne. 0) then
        call DiagonalizerList_Arpack_Allocate(DiagonalizerList(i) % e, "diagonalizerList."//childName)
      else
        error stop "DiagonalizerList element must be one of: lapack, arpack"
      end if

      DiagonalizerList(i) % e % path = "diagonalizerList."//childName
      DiagonalizerList(i) % e % ApplyMatOnVec => input(i) % ApplyMatOnVec
      DiagonalizerList(i) % e % dim = input(i) % dim
      call DiagonalizerList(i) % e % Fabricate()

    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Setup()

    integer(I32) :: i

    do i = 1, size(DiagonalizerList)
      call DiagonalizerList(i) % e % Setup()
    end do

  end subroutine

end submodule
