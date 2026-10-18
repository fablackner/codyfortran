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
    use M_Grid_Ylm

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

    complex(R64) :: dipole

    if (.not. (toScreenQ .or. (filename .ne. "") .or. present(dipole_))) return

    dipole = CosThetaExpectation(rdm1, orbs, Grid_Ylm_radialPoints)

    if (present(dipole_)) dipole_ = real(dipole, kind=R64)

    call Dump(" Dipole moment <z>: ", dipole, filename, toScreenQ, time)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Calculates and dumps the z component of the dipole acceleration in the
!> acceleration form (Ehrenfest theorem for a nuclear charge Z at the origin):
!> \[
!> \ddot{d}(t) = -\langle \partial_z V \rangle = -Z \langle z / r^3 \rangle
!>             = -Z \sum_{i,j} \rho_{ij} \langle \phi_j | \cos\theta / r^2 | \phi_i \rangle
!> \]
!> The force of the external laser field, -N E(t), is not included (it only
!> contributes at the driving frequency). Electron-electron forces cancel.
  subroutine PrinterObservableYlm_DumpDipoleAcceleration(rdm1, orbs, charge, filename, toScreenQ, time, &
                                                         acceleration_)
    use M_Grid_Ylm

    !> One-body reduced density matrix in orbital representation.
    complex(R64), intent(in), contiguous :: rdm1(:, :)
    !> Orbital coefficients on the grid.
    complex(R64), intent(in), contiguous :: orbs(:, :)
    !> Nuclear charge Z of the attractive -Z/r potential.
    real(R64), intent(in) :: charge
    !> Output filename. If empty, no file output is produced.
    character(len=*), intent(in) :: filename
    !> Flag controlling whether to print output to screen.
    logical, intent(in) :: toScreenQ
    !> Current time value to be written with the data.
    real(R64), intent(in) :: time
    !> Optional output of the computed dipole acceleration.
    real(R64), intent(out), optional :: acceleration_

    complex(R64) :: acceleration

    if (.not. (toScreenQ .or. (filename .ne. "") .or. present(acceleration_))) return

    acceleration = -charge * CosThetaExpectation(rdm1, orbs, 1.0_R64 / Grid_Ylm_radialPoints**2)

    if (present(acceleration_)) acceleration_ = real(acceleration, kind=R64)

    call Dump(" Dipole acceleration: ", acceleration, filename, toScreenQ, time)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Expectation value of g(r) cos(theta) for the state given by rdm1 and orbs,
!> with g given on the radial grid points.
!> cos(theta) = sqrt(4 pi / 3) Y_10, so the operator is applied as a spatial
!> product with a potential-like array holding only the (l=1, m=0) component.
  function CosThetaExpectation(rdm1, orbs, radialFunction) result(res)
    use M_Utils_Constants
    use M_Grid
    use M_Grid_Ylm
    use M_Orbs

    complex(R64), intent(in), contiguous :: rdm1(:, :)
    complex(R64), intent(in), contiguous :: orbs(:, :)
    real(R64), intent(in), contiguous :: radialFunction(:)
    complex(R64) :: res

    integer(I32), parameter :: lmaxOp = 1
    integer(I32), parameter :: mmaxOp = 0

    integer(I32) :: i1, j1, occupation, nO, nOS, nG, nRad
    complex(R64), allocatable :: op(:), opLm(:), opOrbs(:, :)

    nO = size(rdm1, 1)
    nG = Grid_nPoints
    nOS = Orbs_nOrbsInState
    nRad = Grid_Ylm_nRadial
    occupation = nO / nOS

    allocate (op(Grid_Ylm_NChannels(lmaxOp, mmaxOp) * nRad))
    allocate (opLm(nRad))
    op(:) = 0.0_R64

    opLm(:) = sqrt(4.0_R64 * PI / 3.0_R64) * radialFunction(:)
    call Grid_Ylm_SetLmComponent(op, 1, 0, opLm, mmax_=mmaxOp)

    allocate (opOrbs(nG, nOS))

    do i1 = 1, nOS
      call Grid_Ylm_SpatialProduct(opOrbs(:, i1), op, orbs(:, i1), lMaxConjg_=lmaxOp, mMaxConjg_=mmaxOp)
    end do

    res = 0.0_R64

    do j1 = 1, nOS
      do i1 = 1, nOS
        res = res + occupation * rdm1(i1, j1) * Grid_InnerProduct(orbs(:, j1), opOrbs(:, i1))
      end do
    end do

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Print a (time, value) observable to screen and/or append it to a file.
  subroutine Dump(label, value, filename, toScreenQ, time)
    character(len=*), intent(in) :: label
    complex(R64), intent(in) :: value
    character(len=*), intent(in) :: filename
    logical, intent(in) :: toScreenQ
    real(R64), intent(in) :: time

    integer(I32) :: io, istat
    character(256) :: msg

    if (toScreenQ) then
      write (*, *)
      write (*, '(A, E20.10E3)') label, real(value, kind=R64)
      write (*, *)
    end if

    if (filename .eq. "") return

    open (newunit=io, file=filename, status="unknown", position="append", iostat=istat, iomsg=msg)
    if (istat .ne. 0) print *, "Failed to open ", filename, "! iomsg = "//trim(msg)
    if (istat .ne. 0) error stop

    write (io, '(3E20.10E3)') time, value

    close (io)

  end subroutine

end module
