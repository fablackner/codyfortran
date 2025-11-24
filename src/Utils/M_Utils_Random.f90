! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Simple random number helpers.
module M_Utils_Random
  use M_Utils_Types
  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Random_setSeed()

    integer :: seedSize, seedClock
    integer, allocatable :: seed(:)
    integer :: k

    call random_seed(size=seedSize)
    allocate (seed(seedSize))

    call system_clock(count=seedClock)

    seed = seedClock + [(k, k=1, seedSize)] * 37

    call random_seed(put=seed)

    deallocate (seed)

  end subroutine Random_setSeed

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function Random_getNumber(a, b) result(res)
    real(R64) :: res
    real(R64), intent(in) :: a
    real(R64), intent(in) :: b

    real(R64) :: randVal

    call random_number(randVal)
    res = a + (b - a) * randVal

  end function Random_getNumber

end module
