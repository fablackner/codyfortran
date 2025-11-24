! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> GSL-backed helpers for combinations without repetition.
!>
!> Wraps `gsl_combination_*` to enumerate combinations and export them into
!> Fortran arrays for downstream indexing tasks.
module M_Utils_CombinationGslLib
  use M_Utils_Types
  use, intrinsic :: iso_c_binding

  implicit none

  private :: GSL_SUCCESS, GSL_FAILURE
  integer(c_int), parameter :: GSL_SUCCESS = 0
  integer(c_int), parameter :: GSL_FAILURE = -1

  interface
    pure function GSL_COMBINATION_ALLOC(n, k) bind(C, name="gsl_combination_alloc") result(comb)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: comb
      integer(c_size_t), intent(in), value :: n, k
    end function

    pure subroutine GSL_COMBINATION_INIT_FIRST(comb) bind(C, name="gsl_combination_init_first")
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: comb
    end subroutine

    pure function GSL_COMBINATION_GET(comb, i) bind(C, name="gsl_combination_get") result(element)
      use, intrinsic :: iso_c_binding
      integer(c_size_t) :: element
      type(c_ptr), intent(in), value :: comb
      integer(c_size_t), intent(in), value :: i
    end function

    pure function GSL_COMBINATION_NEXT(comb) bind(C, name="gsl_combination_next") result(status)
      use, intrinsic :: iso_c_binding
      integer(c_int) :: status
      type(c_ptr), intent(in), value :: comb
    end function

    pure subroutine GSL_COMBINATION_FREE(comb) bind(C, name="gsl_combination_free")
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: comb
    end subroutine
  end interface

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure subroutine CombinationGslLib_CombiNoRepeat(co, n)
    use, intrinsic :: iso_c_binding

    integer(I32), intent(out), contiguous :: co(:, :)
    integer(I32), intent(in) :: n

    integer(I64) :: i, k
    integer(I32) :: c
    type(c_ptr) :: comb
    integer(c_int) :: status

    k = size(co, 1)

    comb = GSL_COMBINATION_ALLOC(int(n, c_size_t), int(k, c_size_t))
    if (.not. c_associated(comb)) then
      error stop "GSL_COMBINATION_ALLOC failed"
    end if

    call GSL_COMBINATION_INIT_FIRST(comb)

    c = 1
    do
      do i = 1, k
        co(i, c) = int(GSL_COMBINATION_GET(comb, i - 1), kind=I32) + 1
      end do
      c = c + 1

      status = GSL_COMBINATION_NEXT(comb)
      if (status .ne. GSL_SUCCESS) exit
    end do

    call GSL_COMBINATION_FREE(comb)

  end subroutine

end module
