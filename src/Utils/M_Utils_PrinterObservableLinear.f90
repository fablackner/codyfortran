! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Module M_Utils_PrinterObservable provides functionality for calculating and visualizing
!> quantum mechanical observables in many-body systems.
!> It handles the output of various physical quantities such as normalization,
!> dipole moments, orthonormality, and energy matrices.
!>
!> Linear-geometry specialization for common observables along 1D grids.
module M_Utils_PrinterObservableLinear
  use M_Utils_Types

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Calculates and dumps the dipole moment of the quantum state.
!> Computes the expectation value of the position operator using:
!> \[
!> \langle \hat{x} \rangle = \sum_{i,j} \rho_{ij} \int \phi_i(x) \, x \, \phi_j^*(x) \, dx
!> \]
!> where \(\rho_{ij}\) is the one-body RDM and \(\phi_i(x)\) are the orbitals.
  subroutine PrinterObservableLinear_DumpDipole(rdm1, orbs, filename, toScreenQ, time)
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
    !> Current time value to be written with the dipole data.
    real(R64), intent(in) :: time

    integer(I32) :: i1, j1
    complex(R64) :: dipole
    integer(I32) :: occupation
    integer(I32) :: io, istat, nO, nOS, nG
    complex(R64), allocatable :: xOrbs(:, :)
    logical :: needDipole
    character(256) :: msg

    nO = size(rdm1, 1)
    nG = Grid_nPoints
    nOS = Orbs_nOrbsInState
    occupation = nO / nOS

    needDipole = toScreenQ .or. (filename .ne. "")
    if (.not. needDipole) return

    allocate (xOrbs(nG, nOS))

    do i1 = 1, nOS
      xOrbs(:, i1) = Grid_Linear_xCoord(:) * orbs(:, i1)
    end do

    dipole = 0.0_R64

    do j1 = 1, nOS
      do i1 = 1, nOS
        dipole = dipole + occupation * rdm1(i1, j1) * Grid_InnerProduct(orbs(:, j1), xOrbs(:, i1))
      end do
    end do

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Print to Screen
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (toScreenQ) then

      write (*, *)
      write (*, *) "------------------------"
      write (*, *) "dipole moment/velocity/acceleration:"
      write (*, *) "------------------------"
      write (*, *)

      write (*, '(A, E20.10E3)') "Dipole moment: ", real(dipole, kind=R64)

      write (*, *)

    end if

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Print to file
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (filename .eq. "") return

    open (newunit=io, file=filename, status="unknown", position="append", iostat=istat, iomsg=msg)
    if (istat .ne. 0) print *, "Failed to open ", filename, "! iomsg = "//trim(msg)
    if (istat .ne. 0) error stop

    write (io, '(2E20.10E3)') time, dipole

    close (io)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Calculates and dumps the potential energy function on the spatial grid.
!> Extracts the external potential \(V(x)\) by applying the potential operator
!> to a constant function and outputs its values at each grid point.
  subroutine PrinterObservableLinear_DumpPotentialFullExpansion(filename, toScreenQ, time)
    use M_Grid
    use M_Grid_Linear
    use M_Method_Mb
    use M_Method_Mb_OrbBased
    use M_SysPotential

    !> Output filename. If empty, no file output is produced.
    character(len=*), intent(in) :: filename
    !> Flag controlling whether to print output to screen.
    logical, intent(in) :: toScreenQ
    !> Current time value for evaluating time-dependent potentials.
    real(R64), intent(in) :: time

    integer(I32) :: iGrid
    integer(I32) :: io, istat, nG
    character(256) :: msg

    complex(R64) :: orbDummy(Grid_nPoints, 1)
    complex(R64) :: potential(Grid_nPoints, 1)

    nG = Grid_nPoints

    ! Set dummy orbital to all ones
    orbDummy(:, :) = 1.0_R64

    ! Initialize potential array to zeros
    potential(:, :) = 0.0_R64

    ! Apply potential operator to get V(x) values
    call Method_Mb_OrbBased_ApplyPotentialOp(potential, orbDummy, time)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Print to Screen
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (toScreenQ) then

      write (*, *)
      write (*, *) "------------------------"
      write (*, *) "potential on grid:      "
      write (*, *) "------------------------"
      write (*, *)

      do iGrid = 1, nG
        write (*, '(2E20.10E3)') Grid_Linear_xCoord(iGrid), dble(potential(iGrid, 1))
      end do

    end if

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Print to file
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (filename .eq. "") return

    open (newunit=io, file=filename, status="unknown", position="append", iostat=istat, iomsg=msg)
    if (istat .ne. 0) print *, "Failed to open ", filename, "! iomsg = "//trim(msg)
    if (istat .ne. 0) error stop

    do iGrid = 1, nG
      write (io, '(2E20.10E3)') Grid_Linear_xCoord(iGrid), dble(potential(iGrid, 1))
    end do

    write (io, *)
    close (io)

  end subroutine

end module
