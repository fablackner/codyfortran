! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Utilities for marking unused variables to suppress compiler warnings.
module M_Utils_UnusedVariables
  use M_Utils_Types

  implicit none

  interface UnusedVariables_MarkOld
    module procedure UnusedVariables_Mark1
    module procedure UnusedVariables_Mark2
    module procedure UnusedVariables_Mark3
    module procedure UnusedVariables_Mark4
    module procedure UnusedVariables_Mark5
  end interface

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure subroutine UnusedVariables_Mark(dummy1, dummy2, dummy3, dummy4, dummy5)
    class(*), intent(in), optional :: dummy1(..)
    class(*), intent(in), optional :: dummy2(..)
    class(*), intent(in), optional :: dummy3(..)
    class(*), intent(in), optional :: dummy4(..)
    class(*), intent(in), optional :: dummy5(..)
    integer, allocatable :: i(:)
    i = shape(dummy1) + shape(dummy2) + shape(dummy3) + shape(dummy4) + shape(dummy5)
    error stop "UnusedVariables_Mark: This function shall never be called."
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure subroutine UnusedVariables_Mark1(dummy)
    class(*), intent(in) :: dummy(..)
    integer, allocatable :: i(:)
    i = shape(dummy)
    error stop "UnusedVariables_Mark: This function shall never be called."
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure subroutine UnusedVariables_Mark2(dummy1, dummy2)
    class(*), intent(in) :: dummy1(..)
    class(*), intent(in) :: dummy2(..)
    integer, allocatable :: i(:)
    i = shape(dummy1) + shape(dummy2)
    error stop "UnusedVariables_Mark: This function shall never be called."
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure subroutine UnusedVariables_Mark3(dummy1, dummy2, dummy3)
    class(*), intent(in) :: dummy1(..)
    class(*), intent(in) :: dummy2(..)
    class(*), intent(in) :: dummy3(..)
    integer, allocatable :: i(:)
    i = shape(dummy1) + shape(dummy2) + shape(dummy3)
    error stop "UnusedVariables_Mark: This function shall never be called."
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure subroutine UnusedVariables_Mark4(dummy1, dummy2, dummy3, dummy4)
    class(*), intent(in) :: dummy1(..)
    class(*), intent(in) :: dummy2(..)
    class(*), intent(in) :: dummy3(..)
    class(*), intent(in) :: dummy4(..)
    integer, allocatable :: i(:)
    i = shape(dummy1) + shape(dummy2) + shape(dummy3) + shape(dummy4)
    error stop "UnusedVariables_Mark: This function shall never be called."
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure subroutine UnusedVariables_Mark5(dummy1, dummy2, dummy3, dummy4, dummy5)
    class(*), intent(in) :: dummy1(..)
    class(*), intent(in) :: dummy2(..)
    class(*), intent(in) :: dummy3(..)
    class(*), intent(in) :: dummy4(..)
    class(*), intent(in) :: dummy5(..)
    integer, allocatable :: i(:)
    i = shape(dummy1) + shape(dummy2) + shape(dummy3) + shape(dummy4) + shape(dummy5)
    error stop "UnusedVariables_Mark: This function shall never be called."
  end subroutine

end module
