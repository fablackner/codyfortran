! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Module M_Utils_PrinterSpectrum provides functionality for calculating and visualizing
!> quantum mechanical spectrums in many-body systems.
!> It handles the output of various physical quantities such as normalization,
!> dipole moments, orthonormality, and energy matrices.
!>
!> Focused on exporting eigenvalue spectra from configured diagonalizers and
!> related summary quantities to disk and/or standard output.
module M_Utils_PrinterSpectrum
  use M_Utils_Types

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine PrinterSpectrum_DumpSpectrum(filename, toScreenQ)
    use M_DiagonalizerList

    character(len=*), intent(in) :: filename
    logical, intent(in)          :: toScreenQ

    integer(I32) :: i, io, istat
    character(256) :: msg

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Print to Screen
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (toScreenQ) then

      do i = 1, DiagonalizerList(1) % e % nFound
        write (*, '(1X,I0,E20.10)') i, DiagonalizerList(1) % e % evals(i)
      end do

    end if

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Print to file
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (filename .eq. "") return

    open (newunit=io, file=filename, status="unknown", position="append", iostat=istat, iomsg=msg)
    if (istat .ne. 0) print *, "Failed to open ", filename, "! iomsg = "//trim(msg)
    if (istat .ne. 0) error stop

    do i = 1, DiagonalizerList(1) % e % nFound
      write (io, '(I0,E20.10)') i, DiagonalizerList(1) % e % evals(i)
    end do

    write (io, *)
    close (io)

  end subroutine

end module
