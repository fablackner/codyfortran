! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Module M_Utils_PrinterRdm provides functionality for printing and visualizing
!> reduced density matrices (RDMs) in quantum many-body systems.
!> It handles the output of one-body, two-body, and three-body RDMs
!> in both orbital basis and spatial grid representations.
!>
!> Routines here assemble RDMs from the active method/state, optionally enforce
!> symmetry/trace constraints, and write results to portable files for analysis.
module M_Utils_PrinterRdm
  use M_Utils_Types

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Dumps the one-body reduced density matrix to a file and/or screen output.
!> Outputs the full matrix elements \(\rho_{ij}\) for all pairs of orbitals i,j.
  subroutine PrinterRdm_DumpRdm1(rdm1, filename, toScreenQ)

    !> One-body reduced density matrix in orbital representation.
    complex(R64), intent(in), contiguous :: rdm1(:, :)
    !> Output filename. If empty, no file output is produced.
    character(len=*), intent(in) :: filename
    !> Flag controlling whether to print output to screen.
    logical, intent(in) :: toScreenQ

    integer(I32) :: io, istat, nO, i, j
    character(256) :: msg

    nO = size(rdm1, 1)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Print to Screen
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (toScreenQ) then

      write (*, *)
      write (*, *) "------------------------"
      write (*, *) "oneRdm in orbital basis:"
      write (*, *) "------------------------"
      write (*, *)

      do j = 1, nO
        do i = 1, nO
          write (*, '(2I6, *(F10.6))') i, j, rdm1(i, j)
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

    write (io, *)
    do j = 1, nO
      do i = 1, nO
        write (io, '(2I6, *(F10.6))') i, j, rdm1(i, j)
      end do
      write (io, *)
    end do

    close (io)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Dumps the absolute values of the one-body reduced density matrix elements
!> in matrix form to a file and/or screen output at a given time.
!> Provides a visualization of the RDM magnitude across all orbital pairs.
  subroutine PrinterRdm_DumpRdm1AbsMatrixForm(rdm1, filename, toScreenQ, time)

    !> One-body reduced density matrix in orbital representation.
    complex(R64), intent(in), contiguous :: rdm1(:, :)
    !> Output filename. If empty, no file output is produced.
    character(len=*), intent(in) :: filename
    !> Flag controlling whether to print output to screen.
    logical, intent(in) :: toScreenQ
    !> Current time value to be written with the RDM data.
    real(R64), intent(in) :: time

    integer(I32) :: io, istat, nO, i, j
    character(256) :: msg

    nO = size(rdm1, 1)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Print to Screen
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (toScreenQ) then

      write (*, *)
      write (*, *) "------------------------"
      write (*, *) "absolute of oneRdm in matrix form:"
      write (*, *) "------------------------"
      write (*, *)

      write (*, *)
      do i = 1, nO
        write (*, '(*(F10.6))') (abs(rdm1(i, j)), j=1, nO)
      end do

    end if

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Print to file
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (filename .eq. "") return

    open (newunit=io, file=filename, status="unknown", position="append", iostat=istat, iomsg=msg)
    if (istat .ne. 0) print *, "Failed to open ", filename, "! iomsg = "//trim(msg)
    if (istat .ne. 0) error stop

    write (io, '(*(F10.6))') time

    do i = 1, nO

      write (io, '(*(F10.6))') (abs(rdm1(i, j)), j=1, nO)

    end do

    write (io, *)
    close (io)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Dumps the two-body reduced density matrix to a file and/or screen output.
!> Outputs the matrix elements \(\rho_{i1,i2,j1,j2}\) for all combinations
!> of orbital indices.
  subroutine PrinterRdm_DumpRdm2(rdm2, filename, toScreenQ)

    !> Two-body reduced density matrix in orbital representation.
    complex(R64), intent(in), contiguous :: rdm2(:, :, :, :)
    !> Output filename. If empty, no file output is produced.
    character(len=*), intent(in) :: filename
    !> Flag controlling whether to print output to screen.
    logical, intent(in) :: toScreenQ

    integer(I32) :: io, istat, nO, i1, j1, i2, j2
    integer(I32) :: jCount, iCount
    character(256) :: msg

    nO = size(rdm2, 1)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Print to Screen
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (toScreenQ) then

      write (*, *)
      write (*, *) "------------------------"
      write (*, *) "rdm2 in orbital basis:"
      write (*, *) "------------------------"
      write (*, *)

      jCount = 0
      do j1 = 1, nO
        do j2 = 1, nO
          jCount = jCount + 1

          iCount = 0
          do i1 = 1, nO
            do i2 = 1, nO
              iCount = iCount + 1

              write (*, '(6I4, *(F10.6))') i1, i2, j1, j2, iCount, jCount, rdm2(i1, i2, j1, j2)

            end do
          end do
          write (*, *)
        end do
      end do

    end if

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Print to file
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (filename .eq. "") return

    open (newunit=io, file=filename, status="unknown", position="append", iostat=istat, iomsg=msg)
    if (istat .ne. 0) print *, "Failed to open ", filename, "! iomsg = "//trim(msg)
    if (istat .ne. 0) error stop

    jCount = 0
    do j1 = 1, nO
      do j2 = 1, nO
        jCount = jCount + 1

        iCount = 0
        do i1 = 1, nO
          do i2 = 1, nO
            iCount = iCount + 1

            write (io, '(6I4, *(F10.6))') i1, i2, j1, j2, iCount, jCount, rdm2(i1, i2, j1, j2)
          end do
        end do
        write (io, *)
      end do
    end do

    close (io)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Computes and dumps symmetry properties of the two-body reduced density matrix
!> for fermionic systems. Evaluates norm, spin, permutation symmetry, hermiticity,
!> and singlet trace conditions.
!> These properties provide important validation checks for RDM correctness.
  subroutine PrinterRdm_DumpRdm2Symmetries(rdm2, filename, toScreenQ)
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    !> Two-body reduced density matrix in orbital representation.
    complex(R64), intent(in), contiguous :: rdm2(:, :, :, :)
    !> Output filename. If empty, no file output is produced.
    character(len=*), intent(in) :: filename
    !> Flag controlling whether to print output to screen.
    logical, intent(in) :: toScreenQ

    integer(I32) :: io, istat, nO, i1, j1, i2, j2
    character(256) :: msg

    integer(I32) :: nOa, nOb, nBodiesUP

    complex(R64) :: spinTr1, spin, norm, hermSym, permSym, norm1, norm2

    nO = size(rdm2, 1)

    nOa = Method_Mb_OrbBased_nOrbs(1)
    nOb = Method_Mb_OrbBased_nOrbs(2)
    nBodiesUP = Method_Mb_nBodies(1)

    !=======================================================================================
    !2RDMNorm
    !=======================================================================================

    norm = 0.0_R64
    do i1 = 1, nO
      do i2 = 1, nO

        norm = norm + rdm2(i1, i2, i1, i2)

      end do
    end do

    !=======================================================================================
    !2RDMSpin
    !=======================================================================================
    ! see e.g. the spin symmetries on page 41 in the book RDM-mechanics by Mazziotti

    spin = 0.0_R64
    do i1 = 1, nOa
      do i2 = 1, nOa

        spin = spin + rdm2(i1, i2 + nOa, i2, i1 + nOa)

      end do
    end do

    spin = nBodiesUP - spin

    !=======================================================================================
    !(N)Tr1(Daa)=(N-1)Tr1(Dab)
    !=======================================================================================
    ! This is the Singlet trace condition from Phys. Rev. A. 91, 023412 (2015)

    spinTr1 = 0.0_R64
    do i1 = 1, nOa
      do j1 = 1, nOa
        norm1 = 0.0_R64
        norm2 = 0.0_R64
        do i2 = 1, nOa

          norm1 = norm1 + rdm2(i1, i2 + nOa, i1, i2 + nOa)
          norm2 = norm2 + rdm2(i1, i2 + nOa, i2, i1 + nOa)

        end do

        spinTr1 = spinTr1 + abs(norm1 - Method_Mb_nBodies(1) * norm2)
      end do
    end do

    !=======================================================================================
    !Permutation-symmetry
    !=======================================================================================

    permSym = 0.0_R64
    do j2 = 1, nO
      do j1 = 1, nO
        do i2 = 1, nO
          do i1 = 1, nO

            permSym = permSym + &
                      abs(rdm2(i1, i2, j1, j2) + rdm2(i2, i1, j1, j2)) + &
                      abs(rdm2(i1, i2, j1, j2) + rdm2(i1, i2, j2, j1)) + &
                      abs(rdm2(i1, i2, j1, j2) - rdm2(i2, i1, j2, j1))

          end do
        end do
      end do
    end do

    permSym = permSym / (nOa)**4

    !=======================================================================================
    !Hermizity
    !=======================================================================================

    hermSym = 0.0_R64
    do j2 = 1, nO
      do j1 = 1, nO
        do i2 = 1, nO
          do i1 = 1, nO

            hermSym = hermSym + abs(rdm2(i1, i2, j1, j2) - conjg(rdm2(j1, j2, i1, i2)))

          end do
        end do
      end do
    end do

    hermSym = hermSym / (nOa)**4

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Print to Screen
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (toScreenQ) then

      write (*, *)
      write (*, *) "TwoRdmSymmetries:"
      write (*, "(A,E25.16)") " norm:       ", real(norm)
      write (*, "(A,E25.16)") " spin:       ", real(spin)
      write (*, "(A,E25.16)") " singletCond ", real(spinTr1)
      write (*, "(A,E25.16)") " symmetry:   ", real(permSym)
      write (*, "(A,E25.16)") " hermizity:  ", real(hermSym)

    end if

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Print to file
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (filename .eq. "") return

    open (newunit=io, file=filename, status="unknown", position="append", iostat=istat, iomsg=msg)
    if (istat .ne. 0) print *, "Failed to open ", filename, "! iomsg = "//trim(msg)
    if (istat .ne. 0) error stop

    write (io, "(5E25.16)") real(norm), real(spin), real(spinTr1), real(permSym), real(hermSym)

    close (io)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Dumps the three-body reduced density matrix to a file and/or screen output.
!> Outputs the matrix elements \(\rho_{i1,i2,i3,j1,j2,j3}\) for all combinations
!> of orbital indices, with appropriate indexing.
  subroutine PrinterRdm_DumpRdm3(rdm3, filename, toScreenQ)

    !> Three-body reduced density matrix in orbital representation.
    complex(R64), intent(in), contiguous :: rdm3(:, :, :, :, :, :)
    !> Output filename. If empty, no file output is produced.
    character(len=*), intent(in) :: filename
    !> Flag controlling whether to print output to screen.
    logical, intent(in) :: toScreenQ

    integer(I32) :: io, istat, nO, i1, j1, i2, j2, i3, j3
    integer(I32) :: jCount, iCount
    character(256) :: msg

    nO = size(rdm3, 1)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Print to Screen
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (toScreenQ) then

      write (*, *)
      write (*, *) "------------------------"
      write (*, *) "rdm3 in orbital basis:"
      write (*, *) "------------------------"
      write (*, *)

      jCount = 0
      do j1 = 1, nO
        do j2 = 1, nO
          do j3 = 1, nO
            jCount = jCount + 1

            iCount = 0
            do i1 = 1, nO
              do i2 = 1, nO
                do i3 = 1, nO
                  iCount = iCount + 1

                  write (*, '(8I4, *(F10.6))') i1, i2, i3, j1, j2, j3, iCount, jCount, rdm3(i1, i2, i3, j1, j2, j3)

                end do
              end do
            end do
            write (*, *)
          end do
        end do
      end do

    end if

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Print to file
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (filename .eq. "") return

    open (newunit=io, file=filename, status="unknown", position="append", iostat=istat, iomsg=msg)
    if (istat .ne. 0) print *, "Failed to open ", filename, "! iomsg = "//trim(msg)
    if (istat .ne. 0) error stop

    jCount = 0
    do j1 = 1, nO
      do j2 = 1, nO
        do j3 = 1, nO
          jCount = jCount + 1

          iCount = 0
          do i1 = 1, nO
            do i2 = 1, nO
              do i3 = 1, nO
                iCount = iCount + 1

                write (io, '(8I4, *(F10.6))') i1, i2, i3, j1, j2, j3, iCount, jCount, rdm3(i1, i2, i3, j1, j2, j3)
              end do
            end do
          end do
          write (io, *)
        end do
      end do
    end do

    close (io)

  end subroutine

end module
