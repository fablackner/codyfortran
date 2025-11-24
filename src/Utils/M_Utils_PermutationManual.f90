! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Pure Fortran permutation generator producing all n! orderings.
!>
!> Recursively builds permutations by fixing the last position and permuting
!> the remaining entries, writing results into a (n, n!) matrix.
module M_Utils_PermutationManual
  use M_Utils_Types

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure recursive subroutine PermutationManual_Permutation(mat)
    use M_Utils_SfManual

    integer(I32), intent(out), contiguous  :: mat(:, :)

    integer(I32), allocatable ::  seq(:)
    integer(I32) :: i, j, n, count, k, fac

    n = size(mat, 1)
    fac = SfManual_Factorial(n)

    if (n .eq. 1) then
      mat(1, 1) = 1
      return
    end if

    ! calculate the Permutations with n at the last position
    mat(n, 1:fac / n) = n
    call PermutationManual_Permutation(mat(1:n - 1, 1:fac / n))

    ! mat(1:n-1, 1:fac/n) contains now the Permutations of n-1 elements

    ! the remaining Permutations can be booted by keeping the last
    ! value fixed and permute the first n-1 elements using mat(1:n-1, 1:fac/n)
    ! we explit the index gathering capabilities of fortran

    count = fac / n + 1
    do i = n - 1, 1, -1

      seq = [(j, j=1, n)]
      seq(n) = i
      seq(i) = n

      do k = 1, fac / n
        mat(n, count) = i
        mat(1:n - 1, count) = seq(mat(1:n - 1, k))
        count = count + 1
      end do

    end do

  end subroutine

end module
