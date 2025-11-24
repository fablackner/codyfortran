! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
module M_Utils_String
  use M_Utils_Types

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function String_InsertInteger(str1, i, str2) result(res)
    character(len=:), allocatable          :: res
    character(len=*), intent(in)           :: str1
    integer(I32), intent(in)               :: i
    character(len=*), intent(in)           :: str2

    character(range(i) + 2) :: tmp

    write (tmp, '(i0)') i

    res = str1//trim(tmp)//str2

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function String_toStr(i) result(res)
    character(len=:), allocatable          :: res
    integer(I32), intent(in)               :: i

    character(len=32)                      :: tmp

    write (tmp, '(i0)') i
    res = trim(tmp)

  end function

end module
