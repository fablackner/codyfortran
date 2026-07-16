! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Module M_Utils_PrinterObservableYlm provides functionality for calculating
!> and visualizing quantum mechanical observables on Ylm grids.
!>
!> Ylm-geometry specialization for common observables in a spherical
!> harmonics representation.
module M_Utils_PrinterObservableYlm
  use M_Utils_Types

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Calculates and dumps the z component of the dipole moment.
!> Computes the expectation value of the position operator z using:
!> \[
!> \langle \hat{z} \rangle = \sum_{i,j} \rho_{ij} \langle \phi_j | z | \phi_i \rangle
!> \]
!> where \(\rho_{ij}\) is the one-body RDM and \(\phi_i\) are the orbitals.
!> In the Ylm expansion z = r \sqrt{4\pi/3} Y_{10}, so z is applied as a
!> spatial product with a potential-like array holding only the (l=1, m=0)
!> component.
  subroutine PrinterObservableYlm_DumpDipole(rdm1, orbs, filename, toScreenQ, time, dipole_)
    use M_Utils_Constants
    use M_Grid
    use M_Grid_Ylm
    use M_Orbs

    !> One-body reduced density matrix in orbital representation.
    complex(R64), intent(in), contiguous :: rdm1(:, :)
    !> Orbital coefficients on the grid.
    complex(R64), intent(in), contiguous :: orbs(:, :)
    !> Output filename. If empty, no file output is produced.
    character(len=*), intent(in) :: filename
    !> Flag controlling whether to print output to screen.
    logical, intent(in) :: toScreenQ
    !> Current time value to be written with the dipole data.
    real(R64), intent(in) :: time
    !> Optional output of the computed dipole moment <z>.
    real(R64), intent(out), optional :: dipole_

    integer(I32), parameter :: lmaxOp = 1

    integer(I32) :: i1, j1
    complex(R64) :: dipole
    integer(I32) :: occupation
    integer(I32) :: io, istat, nO, nOS, nG, nRad
    complex(R64), allocatable :: zOp(:)
    complex(R64), allocatable :: zLm(:)
    complex(R64), allocatable :: zOrbs(:, :)
    logical :: needDipole
    character(256) :: msg

    nO = size(rdm1, 1)
    nG = Grid_nPoints
    nOS = Orbs_nOrbsInState
    nRad = Grid_Ylm_nRadial
    occupation = nO / nOS

    needDipole = toScreenQ .or. (filename .ne. "") .or. present(dipole_)
    if (.not. needDipole) return

    ! Assemble z = r sqrt(4 pi / 3) Y_10 as a potential-like Ylm array
    allocate (zOp((2 * lmaxOp + 1)**2 * nRad))
    allocate (zLm(nRad))
    zOp(:) = 0.0_R64

    zLm(:) = sqrt(4.0_R64 * PI / 3.0_R64) * Grid_Ylm_radialPoints(:)
    call Grid_Ylm_SetLmComponent(zOp, 1, 0, zLm)

    allocate (zOrbs(nG, nOS))

    do i1 = 1, nOS
      call Grid_Ylm_SpatialProduct(zOrbs(:, i1), zOp, orbs(:, i1), Grid_Ylm_lmax, lmaxOp, Grid_Ylm_lmax)
    end do

    dipole = 0.0_R64

    do j1 = 1, nOS
      do i1 = 1, nOS
        dipole = dipole + occupation * rdm1(i1, j1) * Grid_InnerProduct(orbs(:, j1), zOrbs(:, i1))
      end do
    end do

    if (present(dipole_)) dipole_ = real(dipole, kind=R64)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Print to Screen
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (toScreenQ) then

      write (*, *)
      write (*, '(A, E20.10E3)') " Dipole moment <z>: ", real(dipole, kind=R64)
      write (*, *)

    end if

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Print to file
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (filename .eq. "") return

    open (newunit=io, file=filename, status="unknown", position="append", iostat=istat, iomsg=msg)
    if (istat .ne. 0) print *, "Failed to open ", filename, "! iomsg = "//trim(msg)
    if (istat .ne. 0) error stop

    write (io, '(3E20.10E3)') time, dipole

    close (io)

  end subroutine

end module
