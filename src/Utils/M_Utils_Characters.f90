! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Centralized ANSI escape codes, emojis, and color tokens used by printers.
!>
!> Exposes lazily-allocated strings (bold, colors, resets, emojis) to keep
!> user-facing terminal output consistent across modules.
module M_Utils_Characters

  character(len=:), allocatable :: nl
  character(len=:), allocatable :: bold
  character(len=:), allocatable :: reset
  character(len=:), allocatable :: red
  character(len=:), allocatable :: green
  character(len=:), allocatable :: blue
  character(len=:), allocatable :: purple

  character(len=:), allocatable :: emph
  character(len=:), allocatable :: blink

  character(len=:), allocatable :: colorLogo1
  character(len=:), allocatable :: colorLogo2
  character(len=:), allocatable :: colorLogo3
  character(len=:), allocatable :: colorLogo4
  character(len=:), allocatable :: colorLogo5
  character(len=:), allocatable :: colorLogo6

  character(len=:), allocatable :: emojiOmedetto
  character(len=:), allocatable :: emojiGanbatte
  character(len=:), allocatable :: emojiPresent
  character(len=:), allocatable :: emojiSad
  character(len=:), allocatable :: emojiConfused

end module
