! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Module M_Utils_PrinterNatOrb provides functionality for analyzing and visualizing
!> natural orbital properties in quantum many-body systems.
!> It handles the calculation and output of natural orbital occupation numbers
!> derived from one-body reduced density matrices.
!>
!> Diagonalizes the 1-RDM to obtain natural orbitals/occupations and dumps
!> them in a consistent, analysis-friendly format.
module M_Utils_PrinterNatOrb
  use M_Utils_Types

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Diagonalizes the one-body reduced density matrix to extract natural orbital
!> occupations and dumps them to a file and/or screen output.
!> Natural orbitals are eigenfunctions of the one-body RDM that satisfy:
!> \[
!> \hat{\rho}_1 |\phi_i\rangle = n_i |\phi_i\rangle
!> \]
!> where \(n_i\) are the natural orbital occupation numbers.
  subroutine PrinterNatOrb_DumpNatOrbOccupation(rdm1, filename, toScreenQ, time)
    use M_Utils_RdmBasics

    !> One-body reduced density matrix in orbital representation.
    complex(R64), intent(in), contiguous :: rdm1(:, :)
    !> Output filename. If empty, no file output is produced.
    character(len=*), intent(in) :: filename
    !> Flag controlling whether to print output to screen.
    logical, intent(in) :: toScreenQ
    !> Current time value to be written with the natural orbital data.
    real(R64), intent(in) :: time

    real(R64), allocatable    :: natocc(:)
    complex(R64), allocatable :: natorbs(:, :)

    integer(I32) :: io, istat, nO, i
    character(256) :: msg

    nO = size(rdm1, 1)

    call RdmBasics_DiagonalizeRdm1(natocc, natorbs, rdm1)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Print to Screen
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (toScreenQ) then

      write (*, *)
      write (*, *) "------------------------"
      write (*, *) "natOrb occupation:     "
      write (*, *) "------------------------"
      write (*, *)

      do i = 1, nO
        write (*, '(4g0)') ' Occupation of ', i, ' Orbital is ', natocc(i)
      end do

    end if

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Print to file
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (filename .eq. "") return

    open (newunit=io, file=filename, status="unknown", position="append", iostat=istat, iomsg=msg)
    if (istat .ne. 0) print *, "Failed to open ", filename, "! iomsg = "//trim(msg)
    if (istat .ne. 0) error stop

    write (io, '(*(F10.6))') time, (natocc(i), i=1, nO)

    close (io)

  end subroutine

end module
