! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Pure Fortran generation of combinations with repetition (multisets).
!>
!> Recursively fills a (k, nMultiset) matrix of indices representing all
!> multisets of size k drawn from 1..n with repetition allowed.
module M_Utils_MultisetManual
  use M_Utils_Types

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure recursive subroutine MultisetManual_CombiWithRepeat(co, n)
    use M_Utils_SfManual

    integer(I32), intent(out), contiguous  :: co(:, :)
    integer(I32), intent(in) :: n

    integer(I32) :: i, k, kr, ncombi1, ncombi2, nsub

    k = size(co, 1)

    if (n .eq. 1) then

      co(:, 1) = 1

    else if (n > 1) then

      nsub = SfManual_Binomial((n - 1) + k - 1, k)
      call MultisetManual_CombiWithRepeat(co(:, 1:nsub), n - 1)

      ncombi2 = nsub
      do i = 1, k

        kr = k - i

        ncombi1 = ncombi2
        ncombi2 = ncombi1 + SfManual_Binomial((n - 1) + kr - 1, kr)

        call MultisetManual_CombiWithRepeat(co(1:kr, ncombi1 + 1:ncombi2), n - 1)

        co(kr + 1:k, ncombi1 + 1:ncombi2) = n

      end do

    end if

  end subroutine

end module
