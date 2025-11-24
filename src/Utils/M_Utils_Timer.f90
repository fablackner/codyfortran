! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Lightweight wall-clock timing helper for coarse profiling.
!>
!> Call `Timer_Timer("label")` repeatedly to print elapsed wall-clock time
!> since the previous call together with the provided label.
module M_Utils_Timer
  use M_Utils_Types

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Timer_Timer(routine)

    character(len=*), intent(in) :: routine

    integer(I64), save :: start, finish, ticksPerSec
    real(R64) :: duration

    call system_clock(finish)

    duration = dble(finish - start) / ticksPerSec

    if (ticksPerSec .ne. 0) write (*, '(A, f14.7, A)') 'wallTime: ', duration, '   routine: '//routine
    call system_clock(start, ticksPerSec)

  end subroutine

end module
