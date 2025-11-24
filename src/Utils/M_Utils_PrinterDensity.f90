! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Module M_Utils_PrinterDensity provides functionality for dumping density matrices
!> and density distributions to files and screen output.
!> This module handles one-body and two-body densities in both orbital basis
!> and spatial grid representations.
!>
!> Entry points project density matrices to various representations, apply
!> normalization/consistency checks, and write outputs for post-processing and
!> diagnostics.
module M_Utils_PrinterDensity
  use M_Utils_Types

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Dumps the one-body reduced density matrix diagonal elements to a file and/or screen.
!> Outputs the absolute values of the diagonal elements at a given time.
  subroutine PrinterDensity_DumpOneBodyDensity(rdm1, filename, toScreenQ, time)

    !> One-body reduced density matrix in orbital representation.
    complex(R64), intent(in), contiguous :: rdm1(:, :)
    !> Output filename. If empty, no file output is produced.
    character(len=*), intent(in) :: filename
    !> Flag controlling whether to print output to screen.
    logical, intent(in) :: toScreenQ
    !> Current time value to be written with the density data.
    real(R64), intent(in) :: time

    integer(I32) :: io, istat, nO, i
    character(256) :: msg

    nO = size(rdm1, 1)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Print to Screen
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (toScreenQ) then

      write (*, *)
      write (*, *) "------------------------"
      write (*, *) "one-body density in orbital basis:"
      write (*, *) "------------------------"
      write (*, *)

      write (*, '(*(F10.6))') time, (abs(rdm1(i, i)), i=1, nO)

    end if

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Print to file
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (filename .eq. "") return

    open (newunit=io, file=filename, status="unknown", position="append", iostat=istat, iomsg=msg)
    if (istat .ne. 0) print *, "Failed to open ", filename, "! iomsg = "//trim(msg)
    if (istat .ne. 0) error stop

    write (io, '(*(F10.6))') time, (abs(rdm1(i, i)), i=1, nO)

    close (io)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Dumps the diagonal elements of the two-body reduced density matrix to a file and/or screen.
!> Outputs the values of \(\rho_{ijij}\) for all pairs of orbitals i,j.
  subroutine PrinterDensity_DumpTwoBodyDensity(rdm2, filename, toScreenQ)

    !> Two-body reduced density matrix in orbital representation.
    complex(R64), intent(in), contiguous :: rdm2(:, :, :, :)
    !> Output filename. If empty, no file output is produced.
    character(len=*), intent(in) :: filename
    !> Flag controlling whether to print output to screen.
    logical, intent(in) :: toScreenQ

    integer(I32) :: io, istat, nO, i, j
    character(256) :: msg

    nO = size(rdm2, 1)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Print to Screen
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (toScreenQ) then

      write (*, *)
      write (*, *) "------------------------"
      write (*, *) "two-body density in orbital basis:"
      write (*, *) "------------------------"
      write (*, *)

      do j = 1, nO
        do i = 1, nO

          write (*, '(2I6, 2E20.10E3)') i, j, rdm2(i, j, i, j)

        end do
        write (*, *)
      end do

    end if

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Print to file
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (filename .eq. "") return

    open (newunit=io, file=filename, status="unknown", position="append", iostat=istat, iomsg=msg)
    if (istat .ne. 0) print *, "Failed to open ", filename, "! iomsg = "//trim(msg)
    if (istat .ne. 0) error stop

    do j = 1, nO
      do i = 1, nO

        write (io, '(2I6, 2E20.10E3)') i, j, rdm2(i, j, i, j)

      end do
      write (io, *)
    end do

    close (io)

  end subroutine

end module
