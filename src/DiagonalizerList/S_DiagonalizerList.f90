! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file S_DiagonalizerList.f90
!> @brief Implementation of DiagonalizerList fabrication and setup.
!>
!> @details
!> This submodule contains the factory logic that reads the JSON configuration
!> and instantiates the appropriate eigensolver backends. The name-based
!> dispatch ("lapack" vs "arpack" in the JSON key) determines which concrete
!> implementation is allocated.
!>
!> **Fabrication Flow:**
!> 1. Read number of children under `"diagonalizerList"` JSON path
!> 2. Verify input array size matches JSON structure
!> 3. For each child: allocate backend by name, bind callback, call `Fabricate`
!> 4. Bind global `DiagonalizerList_Setup` to local `Setup` procedure
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_DiagonalizerList) S_DiagonalizerList

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Factory subroutine that creates all diagonalizer instances.
!>
!> @details
!> Reads the `"diagonalizerList"` JSON block and creates one eigensolver
!> per child object. Backend selection is based on substring matching:
!> - Name contains "lapack" → LAPACK dense solver
!> - Name contains "arpack" → ARPACK iterative solver
!>
!> **Error conditions:**
!> - `size(input)` must match number of JSON children
!> - Each JSON child must be an object (not scalar/array)
!> - Child name must contain "lapack" or "arpack"
!>
!> @param[in] input  Array of fabrication inputs with callbacks and dimensions
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
    ! Bind global setup procedure pointer
    !------------------------------------

    DiagonalizerList_Setup => Setup

    !------------------------------------
    ! Allocate container array and
    ! instantiate backends by JSON name
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

      ! Backend selection by name substring matching
      if (index(childName, "lapack") .ne. 0) then
        call DiagonalizerList_Lapack_Allocate(DiagonalizerList(i) % e, "diagonalizerList."//childName)
      else if (index(childName, "arpack") .ne. 0) then
        call DiagonalizerList_Arpack_Allocate(DiagonalizerList(i) % e, "diagonalizerList."//childName)
      else
        error stop "DiagonalizerList element must be one of: lapack, arpack"
      end if

      ! Bind common fields and invoke backend-specific fabrication
      DiagonalizerList(i) % e % path = "diagonalizerList."//childName
      DiagonalizerList(i) % e % ApplyMatOnVec => input(i) % ApplyMatOnVec
      DiagonalizerList(i) % e % dim = input(i) % dim
      call DiagonalizerList(i) % e % Fabricate()

    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Calls Setup on all allocated diagonalizer instances.
!>
!> @details
!> Iterates through `DiagonalizerList(:)` and invokes each backend's `Setup`
!> method. This allows backends to allocate working arrays sized to the
!> problem dimensions established during fabrication.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Setup()

    integer(I32) :: i

    do i = 1, size(DiagonalizerList)
      call DiagonalizerList(i) % e % Setup()
    end do

  end subroutine

end submodule
