! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Lightweight binary I/O for arbitrary array shapes using modern Fortran.
!>
!> Uses unlimited polymorphic, assumed-rank dummy arguments and C interoperability
!> to load/save raw binary data to/from files without imposing a specific type
!> or rank at compile time. Intended for fast checkpoints and diagnostics.
module M_Utils_DataStorage
  use M_Utils_Types

  implicit none

contains

! This is pure magic. I wrote it during a sweaty delirium and I have no idea what is going on.
! All what matters is that it seems to work!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine LoadData(data, filename, elementBits, startPos_)
    use, intrinsic :: iso_c_binding

    type(*), intent(inout), target       :: data(..)
    character(len=*), intent(in)       :: filename
    integer(I32), intent(in)           :: elementBits
    integer(I32), intent(in), optional :: startPos_

    integer(I32) :: elementBytes
    integer(I32) :: nBytes, shift
    type(c_ptr) :: pointerToData
    integer(I08), contiguous, pointer :: bytes(:)
    integer(I08) :: testByte

    integer(I32) :: i, istat, jstat, io
    character(256) :: msg

    elementBytes = elementBits / 8
    nBytes = elementBytes * size(data)

    if (.not. present(startPos_)) shift = 0
    if (present(startPos_)) shift = (startPos_ - 1) * elementBytes

    pointerToData = c_loc(data)
    call c_f_pointer(pointerToData, bytes, [nBytes])

    open (newunit=io, form="unformatted", file=filename, access='stream', status="old", action="read", iostat=istat, iomsg=msg)
    if (istat .ne. 0) print *, trim(msg)
    if (istat .ne. 0) error stop

    do i = 1, nBytes
      read (io, iostat=jstat, pos=i + shift) bytes(i)
      if (jstat .eq. iostat_end) print *, 'more data requested than stored in '//filename
      if (jstat .eq. iostat_end) error stop
    end do

    print *
    write (*, '(7g0)') ' Loaded data from ', filename, ' ranging Byte ', 1 + shift, ' to Byte ', shift + nBytes

    read (io, iostat=jstat) testByte
    if (jstat .eq. iostat_end) print *, 'reached end of '//filename

    close (io)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine SaveData(data, filename, elementBits, startPos_)
    use, intrinsic :: iso_c_binding

    type(*), intent(in), target        :: data(..)
    character(len=*), intent(in)       :: filename
    integer(I32), intent(in)           :: elementBits
    integer(I32), intent(in), optional :: startPos_

    integer(I32) :: elementBytes
    integer(I32) :: nBytes, shift
    type(c_ptr) :: pointerToData
    integer(I08), contiguous, pointer :: bytes(:)

    integer(I32) :: i, istat, io
    character(256) :: msg

    elementBytes = elementBits / 8
    nBytes = elementBytes * size(data)

    if (.not. present(startPos_)) shift = 0
    if (present(startPos_)) shift = (startPos_ - 1) * elementBytes

    pointerToData = c_loc(data)
    call c_f_pointer(pointerToData, bytes, [nBytes])

    open (newunit=io, form="unformatted", file=filename, access='stream', status="replace", action="write", iostat=istat, iomsg=msg)

    if (istat .ne. 0) print *, trim(msg)
    if (istat .ne. 0) error stop

    do i = 1, nBytes
      write (io, pos=i + shift) bytes(i)
    end do

    print *
    write (*, '(7g0)') ' Saved data to ', filename, ' ranging byte ', 1 + shift, ' to byte ', shift + nBytes

    close (io)

  end subroutine

end module
