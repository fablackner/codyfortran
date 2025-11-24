! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Friendly, colorized logging helpers for consistent console output.
!>
!> Provides short routines to print headings, status lines, and emojis using
!> the tokens defined in `M_Utils_Characters`.
module M_Utils_Say
  use M_Utils_Types

  implicit none

  character(256), allocatable    :: fabricatedList(:)
  character(256), allocatable    :: setupList(:)
  integer(I64), save :: startTicks = 0_I64

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! UTF-8 byte code source: https://apps.timwhitlock.info/emoji/tables/unicode
  subroutine Say_Hello()
    use M_Utils_Characters
    character(8)  :: date
    character(10) :: time

    call SetStyleCharacters()

    call date_and_time(date=date, time=time)

    write (*, *) "STARTED on ", date(7:8), ".", date(5:6), ".", date(1:4), " at ", time(1:2), ":", time(3:4), ":", time(5:10)
    write (*, *) nl//" Compiler: "//compiler_version()
    write (*, *) "Flags: "//compiler_options()

    call PrintLogoLarge()
    call PrintGanbatte()

    call system_clock(startTicks)
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Say_Goodbye()
    use M_Utils_Characters
    character(8)  :: date
    character(10) :: time
    integer(I64) :: finishTicks, ticksPerSec

    call date_and_time(date=date, time=time)
    call system_clock(finishTicks, ticksPerSec)

    call PrintLogoSmall()
    call PrintOmedetoo()

    write (*, *) "Total elapsed wall time:", real(finishTicks - startTicks, kind=R64) / real(ticksPerSec, kind=R64), nl

    write (*, *) "FINISHED on ", date(7:8), ".", date(5:6), ".", date(1:4), " at ", time(1:2), ":", time(3:4), ":", time(5:10), nl
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Say_Section(label)
    character(len=*), intent(in) :: label
    call PrintSection(label)
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Say_Fabricate(details)
    use M_Utils_Characters
    character(len=*), intent(in) :: details

    call UpdateList(fabricatedList, details)

    if (details .eq. "coeffs") then
      if (.not. any(fabricatedList .eq. "method")) error stop details//": fabricate method first"
    end if

    if (details .eq. "twoRdm") then
      if (.not. any(fabricatedList .eq. "method")) error stop details//": fabricate method first"
    end if

    write (*, "(A)") nl//bold//" fabricate: "//details//reset
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Say_Setup(details)
    use M_Utils_Characters
    character(len=*), intent(in) :: details

    call UpdateList(setupList, details)

    if (details .eq. "absorber") then
      if (.not. any(setupList .eq. "grid")) error stop details//": setup grid first"
    end if

    if (details .eq. "method") then
      if (.not. any(setupList .eq. "coeffs")) error stop details//": setup coeffs first"
      if (.not. any(setupList .eq. "coeffsInit")) error stop details//": setup coeffsInit first"
      if (.not. any(setupList .eq. "orbs")) error stop details//": setup orbs first"
      if (.not. any(setupList .eq. "orbsInit")) error stop details//": setup orbsInit first"
    end if

    write (*, "(A)") nl//" "//purple//" setup: "//details//reset
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine PrintSection(label)
    use M_Utils_Characters
    character(len=*), intent(in) :: label
    integer :: width
    character(len=:), allocatable :: pad

    width = 76
    pad = repeat(' ', max(0, width - 1 - len_trim(label)))

    write (*, *)
    write (*, "(A)") emph//repeat('=', width)//reset
    write (*, "(A)") emph//" "//trim(label)//pad//reset
    write (*, "(A)") emph//repeat('=', width)//reset
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine PrintLogoLarge()
    use M_Utils_Characters
    write (*, *)

    write (*, "(A)") colorLogo1//"             _________       .   ____   _ _  "//reset
    write (*, "(A)") colorLogo2//" __________________  /__  __ \\ / __ \ (°_°) "//reset
    write (*, "(A)") colorLogo3//" _  ___/  __ \  __  /  / / / // \ \ \ \ \ \  "//reset
    write (*, "(A)") colorLogo4//" / /__ / /_/ / /_/ /  /_/ /  \\_/ / | |_/ /  "//reset
    write (*, "(A)") colorLogo5//" \___/ \____/\____/ \,__, /   '~~~'  '~~~~'  "//reset
    write (*, "(A)") colorLogo6//"                   /____/                    "//reset
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine PrintLogoSmall()
    use M_Utils_Characters
    write (*, *)
    write (*, "(A)") colorLogo1//"  .   cody   _ _  "//reset
    write (*, "(A)") colorLogo2//"  \\ / __ \ (°_°) "//reset
    write (*, "(A)") colorLogo3//"  // \ \ \ \ \ \  "//reset
    write (*, "(A)") colorLogo4//"  \\_/ / | |_/ /  "//reset
    write (*, "(A)") colorLogo5//"  '~~~'  '~~~~'   "//reset
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine PrintGanbatte()
    use M_Utils_Characters
    write (*, "(A)") blink//"  Ganbatte!"//emojiGanbatte//reset
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine PrintOmedetoo()
    use M_Utils_Characters
    write (*, "(A)") blink//"  Omedetoo!"//emojiOmedetto//reset
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine UpdateList(list, elementToAdd)
    character(256), intent(inout), allocatable   :: list(:)
    character(len=*), intent(in) :: elementToAdd

    if (.not. allocated(list)) then
      allocate (list(1))
      list(1) = elementToAdd//repeat(' ', 256 - len_trim(elementToAdd))
    else
      list = [list(:), elementToAdd//repeat(' ', 256 - len_trim(elementToAdd))]
    end if
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine SetStyleCharacters()
    use M_Utils_Characters
    character(32) :: stVal
    integer(I32) :: istat = 0
    character(len=:), allocatable :: emojiPrefix

    call get_environment_variable("CODY_PRETTY_PRINT", stVal, status=istat)

    if (istat .eq. 0 .and. trim(stVal) .eq. "OFF") then

      nl = new_line('A')
      bold = ""
      reset = ""
      red = ""
      green = ""
      blue = ""
      purple = ""

      emph = ""
      blink = ""

      colorLogo1 = ""
      colorLogo2 = ""
      colorLogo3 = ""
      colorLogo4 = ""
      colorLogo5 = ""
      colorLogo6 = ""

      emojiOmedetto = ""
      emojiGanbatte = ""
      emojiPresent = ""
      emojiSad = ""
      emojiConfused = ""

    else

      nl = new_line('A')
      bold = achar(27)//"[1m"
      reset = achar(27)//"[0m"
      red = achar(27)//"[41;97m"
      green = achar(27)//"[42;97m"
      blue = achar(27)//"[44;97m"
      purple = achar(27)//"[45;97m"

      emph = achar(27)//"[1;30;103m"
      blink = achar(27)//"[1;5m"

      colorLogo1 = achar(27)//"[1;48;5;226;38;5;52m"
      colorLogo2 = achar(27)//"[1;48;5;190;38;5;53m"
      colorLogo3 = achar(27)//"[1;48;5;154;38;5;54m"
      colorLogo4 = achar(27)//"[1;48;5;118;38;5;55m"
      colorLogo5 = achar(27)//"[1;48;5;82;38;5;56m"
      colorLogo6 = achar(27)//"[1;48;5;46;38;5;57m"

      emojiPrefix = char(int(z"F0"))//char(int(z"9F"))

      emojiOmedetto = emojiPrefix//char(int(z"8F"))//char(int(z"86"))
      emojiGanbatte = emojiPrefix//char(int(z"8D"))//char(int(z"80"))
      ! gift 🎁 U+1F381 = F0 9F 8E 81
      emojiPresent = emojiPrefix//char(int(z"8E"))//char(int(z"81"))
      ! pensive/sad 😔 U+1F614 = F0 9F 98 94
      emojiSad = emojiPrefix//char(int(z"98"))//char(int(z"94"))
      ! confused 😕 U+1F615 = F0 9F 98 95
      emojiConfused = emojiPrefix//char(int(z"98"))//char(int(z"95"))

    end if

  end subroutine

end module
