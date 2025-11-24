! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> GSL-backed helpers for permutations (combinatorics basics).
!>
!> Provides thin wrappers around gsl_permutation_* to build full permutation
!> sets used for indexing and enumeration tasks.
module M_Utils_CombinatoricsBasics
  use M_Utils_Types
  use, intrinsic :: iso_c_binding

  implicit none

  private :: GSL_SUCCESS, GSL_FAILURE
  integer(c_int), parameter :: GSL_SUCCESS = 0
  integer(c_int), parameter :: GSL_FAILURE = -1

  interface
    pure function GSL_PERMUTATION_ALLOC(n) bind(C, name="gsl_permutation_alloc") result(perm)
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: perm
      integer(c_size_t), intent(in), value :: n
    end function

    pure subroutine GSL_PERMUTATION_INIT(perm) bind(C, name="gsl_permutation_init")
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: perm
    end subroutine

    pure function GSL_PERMUTATION_GET(perm, i) bind(C, name="gsl_permutation_get") result(element)
      use, intrinsic :: iso_c_binding
      integer(c_size_t) :: element
      type(c_ptr), intent(in), value :: perm
      integer(c_size_t), intent(in), value :: i
    end function

    pure function GSL_PERMUTATION_NEXT(perm) bind(C, name="gsl_permutation_next") result(status)
      use, intrinsic :: iso_c_binding
      integer(c_int) :: status
      type(c_ptr), intent(in), value :: perm
    end function

    pure subroutine GSL_PERMUTATION_FREE(perm) bind(C, name="gsl_permutation_free")
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: perm
    end subroutine
  end interface

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure subroutine PermutationGslLib_Permutations(mat)
    use, intrinsic :: iso_c_binding

    integer(I32), intent(out), contiguous :: mat(:, :)

    integer(I64) :: i, n
    integer(I32) :: c
    type(c_ptr) :: perm
    integer(c_int) :: status

    n = size(mat, 1)

    perm = GSL_PERMUTATION_ALLOC(int(n, c_size_t))
    if (.not. c_associated(perm)) then
      error stop "GSL_PERMUTATION_ALLOC failed"
      return
    end if

    call GSL_PERMUTATION_INIT(perm)

    c = 1
    do
      do i = 1, n
        mat(i, c) = int(GSL_PERMUTATION_GET(perm, i - 1), kind=I32) + 1
      end do
      c = c + 1

      status = GSL_PERMUTATION_NEXT(perm)
      if (status .ne. GSL_SUCCESS) exit
    end do

    call GSL_PERMUTATION_FREE(perm)

  end subroutine

end module
