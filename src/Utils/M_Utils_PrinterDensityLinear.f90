! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Module M_Utils_PrinterDensity provides functionality for dumping density matrices
!> and density distributions to files and screen output.
!> This module handles one-body and two-body densities in both orbital basis
!> and spatial grid representations.
!>
!> Linear-geometry specialization of density printers with routines tailored
!> to 1D grids and their projections.
module M_Utils_PrinterDensityLinear
  use M_Utils_Types

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Projects and dumps the one-body density from orbital basis to linear grid.
!> Computes the linear representation of the density using:
!> \[
!> \rho(x) = \sum_{i,j} \rho_{ij} \phi_i(x) \phi_j^*(x)
!> \]
!> where \(\rho_{ij}\) is the one-body RDM and \(\phi_i(x)\) are the orbitals.
  subroutine PrinterDensityLinear_DumpOneBodyDensityOnGrid(rdm1, orbs, filename, toScreenQ)
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
    integer(I32) :: iGrid, i1, j1
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
      write (*, *) "one-body density on Grid:"
      write (*, *) "------------------------"
      write (*, *)

      do iGrid = 1, nG

        cTmp = 0.0_R64
        do j1 = 1, nOS
          do i1 = 1, nOS
            cTmp = cTmp + occupation * rdm1(i1, j1) * orbs(iGrid, i1) * conjg(orbs(iGrid, j1))
          end do
        end do

        write (*, '(3E20.10E3)') Grid_Linear_xCoord(iGrid), cTmp
      end do

      write (*, *)

    end if

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Print to file
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (filename .eq. "") return

    open (newunit=io, file=filename, status="unknown", position="append", iostat=istat, iomsg=msg)
    if (istat .ne. 0) print *, "Failed to open ", filename, "! iomsg = "//trim(msg)
    if (istat .ne. 0) error stop

    do iGrid = 1, nG

      cTmp = 0.0_R64
      do j1 = 1, nOS
        do i1 = 1, nOS
          cTmp = cTmp + occupation * rdm1(i1, j1) * orbs(iGrid, i1) * conjg(orbs(iGrid, j1))
        end do
      end do

      write (io, '(3E20.10E3)') Grid_Linear_xCoord(iGrid), cTmp
    end do

    write (io, *)

    close (io)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Projects and dumps the two-body density from orbital basis to spatial grid.
!> Computes the spatial representation of the two-body density using:
!> \[
!> \rho(x,y) = \sum_{i1,i2,j1,j2} \rho_{i1,i2,j1,j2} \phi_{i1}(x) \phi_{i2}(x) \phi_{j1}^*(y) \phi_{j2}^*(y)
!> \]
!> where \(\rho_{i1,i2,j1,j2}\) is the two-body RDM and \(\phi_i(x)\) are the orbitals.
  subroutine PrinterDensityLinear_DumpTwoBodyDensityOnGrid(rdm2, orbs, filename, toScreenQ)
    use M_Grid
    use M_Grid_Linear
    use M_Orbs

    !> Two-body reduced density matrix in orbital representation.
    complex(R64), intent(in), contiguous :: rdm2(:, :, :, :)
    !> Orbital coefficients on the grid.
    complex(R64), intent(in), contiguous :: orbs(:, :)
    !> Output filename. If empty, no file output is produced.
    character(len=*), intent(in) :: filename
    !> Flag controlling whether to print output to screen.
    logical, intent(in) :: toScreenQ

    complex(R64) :: cTmp
    integer(I32) :: iGrid, jGrid
    integer(I32) :: i1, j1, i2, j2
    integer(I32) :: io, istat, nO, nOS, nG
    integer(I32) :: occupation
    character(256) :: msg

    complex(R64), allocatable :: rdm2spatial(:, :, :, :)

    nO = size(rdm2, 1)
    nG = Grid_nPoints
    nOS = Orbs_nOrbsInState
    occupation = nO / nOS

    allocate (rdm2spatial(nOS, nOS, nOS, nOS))

    rdm2spatial(:, :, :, :) = rdm2(:nOS, :nOS, :nOS, :nOS)

    if (occupation .eq. 2) then
      rdm2spatial(:, :, :, :) = rdm2spatial(:nOS, :nOS, :nOS, :nOS) + &
                                rdm2(nOS + 1:, :nOS, nOS + 1:, :nOS) + &
                                rdm2(:nOS, nOS + 1:, :nOS, nOS + 1:) + &
                                rdm2(nOS + 1:, nOS + 1:, nOS + 1:, nOS + 1:)
    end if

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Print to Screen
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (toScreenQ) then

      write (*, *)
      write (*, *) "------------------------"
      write (*, *) "two-body density on Grid:"
      write (*, *) "------------------------"
      write (*, *)

      do jGrid = 1, nG
        do iGrid = 1, nG

          cTmp = 0.0_R64
          do j2 = 1, nOS
            do j1 = 1, nOS
              do i2 = 1, nOS
                do i1 = 1, nOS

                  cTmp = cTmp + &
                         rdm2spatial(i1, i2, j1, j2) * &
                         orbs(iGrid, i1) * &
                         orbs(iGrid, i2) * &
                         conjg(orbs(jGrid, j1)) * &
                         conjg(orbs(jGrid, j2))

                end do
              end do
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
        do j2 = 1, nOS
          do j1 = 1, nOS
            do i2 = 1, nOS
              do i1 = 1, nOS

                cTmp = cTmp + rdm2spatial(i1, i2, j1, j2) * orbs(iGrid, i1) * orbs(iGrid, i2) &
                       * conjg(orbs(jGrid, j1)) * conjg(orbs(jGrid, j2))

              end do
            end do
          end do
        end do

        write (io, '(4E20.10E3)') Grid_Linear_xCoord(iGrid), Grid_Linear_xCoord(jGrid), cTmp

      end do
      write (io, *)
    end do

    close (io)

  end subroutine

end module
