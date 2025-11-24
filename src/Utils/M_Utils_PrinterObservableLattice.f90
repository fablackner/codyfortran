! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Module M_Utils_PrinterObservable provides functionality for calculating and visualizing
!> quantum mechanical observables in many-body systems.
!> It handles the output of various physical quantities.
!>
!> Linear-geometry specialization for common observables along 1D grids.
module M_Utils_PrinterObservableLattice
  use M_Utils_Types

  implicit none

contains

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Calculates and dumps the center of the particle density.
  !> @details Computes the expectation value of the position operator using:
  !> \[
  !> \langle \hat{x} \rangle = \sum_{i,j} \rho_{ij} \int \phi_i(x) \, x \, \phi_j^*(x) \, dx
  !> \]
  !> where \(\rho_{ij}\) is the one-body RDM and \(\phi_i(x)\) are the orbitals.
  !>
  !> @param[in] rdm1 One-body reduced density matrix in orbital representation.
  !> @param[in] orbs Orbital coefficients on the grid.
  !> @param[in] filename Output filename. If empty, no file output is produced.
  !> @param[in] toScreenQ Flag controlling whether to print output to screen.
  !> @param[in] time Current time value to be written with the data.
  subroutine PrinterObservableLattice_DumpDensityCenter(rdm1, orbs, filename, toScreenQ, time)
    use M_Grid
    use M_Grid_Lattice
    use M_Orbs

    complex(R64), intent(in), contiguous :: rdm1(:, :)
    complex(R64), intent(in), contiguous :: orbs(:, :)
    character(len=*), intent(in) :: filename
    logical, intent(in) :: toScreenQ
    real(R64), intent(in) :: time

    integer(I32) :: i1, i2, i3, j1, j2
    integer(I32) :: iSite
    complex(R64) :: centerX, centerY, centerZ
    integer(I32) :: occupation
    integer(I32) :: io, istat, nO, nOS
    complex(R64), allocatable :: xOrbs(:, :), yOrbs(:, :), zOrbs(:, :)
    logical :: needsCenter
    character(256) :: msg

    nO = size(rdm1, 1)
    nOS = Orbs_nOrbsInState
    occupation = nO / nOS

    needsCenter = toScreenQ .or. (filename .ne. "")
    if (.not. needsCenter) return

    allocate (xOrbs(size(orbs, 1), nOS), yOrbs(size(orbs, 1), nOS), zOrbs(size(orbs, 1), nOS))
    xOrbs = 0.0_R64
    yOrbs = 0.0_R64
    zOrbs = 0.0_R64

    do j1 = 1, nOS
      do i3 = 1, Grid_Lattice_zSize
        do i2 = 1, Grid_Lattice_ySize
          do i1 = 1, Grid_Lattice_xSize
            iSite = Grid_Lattice_code(i1, i2, i3)
            xOrbs(iSite, j1) = i1 * orbs(iSite, j1)
            yOrbs(iSite, j1) = i2 * orbs(iSite, j1)
            zOrbs(iSite, j1) = i3 * orbs(iSite, j1)
          end do
        end do
      end do
    end do

    centerX = 0.0_R64
    centerY = 0.0_R64
    centerZ = 0.0_R64

    do j2 = 1, nOS
      do j1 = 1, nOS
        centerX = centerX + occupation * rdm1(j1, j2) * Grid_InnerProduct(orbs(:, j2), xOrbs(:, j1))
        centerY = centerY + occupation * rdm1(j1, j2) * Grid_InnerProduct(orbs(:, j2), yOrbs(:, j1))
        centerZ = centerZ + occupation * rdm1(j1, j2) * Grid_InnerProduct(orbs(:, j2), zOrbs(:, j1))
      end do
    end do

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Print to Screen
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (toScreenQ) then

      write (*, *)
      write (*, *) "------------------------"
      write (*, *) "Center of particle density:"
      write (*, *) "------------------------"
      write (*, *)

      write (*, '(A, E20.10E3)') "Center X: ", real(centerX, kind=R64)
      write (*, '(A, E20.10E3)') "Center Y: ", real(centerY, kind=R64)
      write (*, '(A, E20.10E3)') "Center Z: ", real(centerZ, kind=R64)

      write (*, *)

    end if

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Print to file
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (filename .eq. "") return

    open (newunit=io, file=filename, status="unknown", position="append", iostat=istat, iomsg=msg)
    if (istat .ne. 0) then
      print *, "Failed to open ", trim(filename), "! iomsg = "//trim(msg)
      error stop
    end if

    write (io, '(4E20.10E3)') time, real(centerX, kind=R64), real(centerY, kind=R64), real(centerZ, kind=R64)

    close (io)

  end subroutine

end module
