! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Module M_Utils_PrinterRdm provides functionality for printing and visualizing
!> reduced density matrices (RDMs) in quantum many-body systems.
!> It handles the output of one-body, two-body, and three-body RDMs
!> in both orbital basis and spatial grid representations.
!>
!> Linear-geometry specialization that projects RDMs to 1D grids and writes
!> compact representations suitable for analysis and plotting.
module M_Utils_PrinterRdmLinear
  use M_Utils_Types

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Projects and dumps the one-body reduced density matrix from orbital basis to spatial grid.
!> Computes the spatial representation of the one-body RDM using:
!> \[
!> \rho(x,y) = \sum_{i,j} \rho_{ij} \phi_i(x) \phi_j^*(y)
!> \]
!> where \(\rho_{ij}\) is the one-body RDM and \(\phi_i(x)\) are the orbitals.
  subroutine PrinterRdmLinear_DumpRdm1FullExpansion(rdm1, orbs, filename, toScreenQ)
    use M_Grid
    use M_Grid_Linear
    use M_Orbs

    !> One-body reduced density matrix in orbital representation.
    complex(R64), intent(in), contiguous :: rdm1(:, :)
    !> Orbital coefficients on the grid.
    complex(R64), intent(in), contiguous :: orbs(:, :)
    !> Output filename. If empty, no file output is produced.
    character(len=*), intent(in) :: filename
    !> Flag controlling whether to print output to screen.
    logical, intent(in) :: toScreenQ

    complex(R64) :: cTmp
    integer(I32) :: iGrid, jGrid, i1, j1
    integer(I32) :: occupation
    integer(I32) :: io, istat, nO, nOS, nG
    character(256) :: msg

    nO = size(rdm1, 1)
    nG = Grid_nPoints
    nOS = Orbs_nOrbsInState
    occupation = nO / nOS

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Print to Screen
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (toScreenQ) then

      write (*, *)
      write (*, *) "------------------------"
      write (*, *) "oneRdm on grid:         "
      write (*, *) "------------------------"
      write (*, *)

      do jGrid = 1, nG
        do iGrid = 1, nG

          cTmp = 0.0_R64
          do j1 = 1, nOS
            do i1 = 1, nOS
              cTmp = cTmp + occupation * rdm1(i1, j1) * orbs(iGrid, i1) * conjg(orbs(jGrid, j1))
            end do
          end do

          write (*, '(4E20.10E3)') Grid_Linear_xCoord(iGrid), Grid_Linear_xCoord(jGrid), cTmp

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

    do jGrid = 1, nG
      do iGrid = 1, nG

        cTmp = 0.0_R64
        do j1 = 1, nOS
          do i1 = 1, nOS
            cTmp = cTmp + occupation * rdm1(i1, j1) * orbs(iGrid, i1) * conjg(orbs(jGrid, j1))
          end do
        end do

        write (io, '(4E20.10E3)') Grid_Linear_xCoord(jGrid), Grid_Linear_xCoord(iGrid), cTmp

      end do
      write (io, *)
    end do

    close (io)

  end subroutine

end module
