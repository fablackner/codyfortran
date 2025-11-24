! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Pure Fortran generation of k-combinations from n without repetition.
!>
!> Implements a recursive splitter based on binomial counts to fill an output
!> matrix of shape (k, nChooseK), where each column is a combination.
module M_Utils_CombinationManual
  use M_Utils_Types

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure recursive subroutine CombinationManual_CombiNoRepeat(co, n)
    use M_Utils_SfManual

    integer(I32), intent(out), contiguous  :: co(:, :)
    integer(I32), intent(in)  :: n

    integer(I32) :: i, k, ncombi, nsub

    k = size(co, 1)

    if (n .eq. k) then

      co(:, 1) = [(i, i=1, n)]

    else if (k .eq. 1) then

      co(1, :) = [(i, i=1, n)]

    else if (n > k) then

      nsub = SfManual_Binomial(n - 1, k)
      call CombinationManual_CombiNoRepeat(co(:, 1:nsub), n - 1)

      ncombi = nsub + SfManual_Binomial(n - 1, k - 1)
      call CombinationManual_CombiNoRepeat(co(1:k - 1, nsub + 1:ncombi), n - 1)
      co(k, nsub + 1:ncombi) = n

    end if

  end subroutine

end module
