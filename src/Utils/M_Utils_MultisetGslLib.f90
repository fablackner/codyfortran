! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> GSL-backed helpers for combinations with repetition (multisets).
!>
!> Wraps `gsl_multiset_*` to enumerate multisets and provide them as Fortran
!> matrices for easy downstream consumption.

module M_Utils_MultisetGslLib
  use M_Utils_Types
  use, intrinsic :: iso_c_binding

  implicit none

  private :: GSL_SUCCESS, GSL_FAILURE
  integer(c_int), parameter :: GSL_SUCCESS = 0
  integer(c_int), parameter :: GSL_FAILURE = -1

  interface
    pure function GSL_MULTISET_ALLOC(n, k) bind(C, name="gsl_multiset_alloc") result(mset)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: mset
      integer(c_size_t), intent(in), value :: n, k
    end function

    pure subroutine GSL_MULTISET_INIT_FIRST(mset) bind(C, name="gsl_multiset_init_first")
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: mset
    end subroutine

    pure function GSL_MULTISET_GET(mset, i) bind(C, name="gsl_multiset_get") result(element)
      use, intrinsic :: iso_c_binding
      integer(c_size_t) :: element
      type(c_ptr), intent(in), value :: mset
      integer(c_size_t), intent(in), value :: i
    end function

    pure function GSL_MULTISET_NEXT(mset) bind(C, name="gsl_multiset_next") result(status)
      use, intrinsic :: iso_c_binding
      integer(c_int) :: status
      type(c_ptr), intent(in), value :: mset
    end function

    pure subroutine GSL_MULTISET_FREE(mset) bind(C, name="gsl_multiset_free")
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: mset
    end subroutine
  end interface

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure subroutine MultisetGslLib_CombiWithRepeat(co, n)
    use, intrinsic :: iso_c_binding

    integer(I32), intent(out), contiguous :: co(:, :)
    integer(I32), intent(in) :: n

    integer(I64) :: i, k
    integer(I32) :: c
    type(c_ptr) :: mset
    integer(c_int) :: status

    k = size(co, 1)

    mset = GSL_MULTISET_ALLOC(int(n, c_size_t), int(k, c_size_t))
    if (.not. c_associated(mset)) then
      error stop "GSL_MULTISET_ALLOC failed"
      return
    end if

    call GSL_MULTISET_INIT_FIRST(mset)

    c = 1
    do
      do i = 1, k
        co(i, c) = int(GSL_MULTISET_GET(mset, i - 1), kind=I32) + 1
      end do
      c = c + 1

      status = GSL_MULTISET_NEXT(mset)
      if (status .ne. GSL_SUCCESS) exit
    end do

    call GSL_MULTISET_FREE(mset)

  end subroutine

end module
