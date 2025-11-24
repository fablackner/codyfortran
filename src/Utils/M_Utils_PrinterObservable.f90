! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Module M_Utils_PrinterObservable provides functionality for calculating and visualizing
!> quantum mechanical observables in many-body systems.
!> It handles the output of various physical quantities such as normalization,
!> dipole moments, orthonormality, and energy matrices.
!>
!> This module computes selected observables from the current method/state and
!> dumps them to files and/or screen with consistent formatting helpers.
module M_Utils_PrinterObservable
  use M_Utils_Types

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Calculates and dumps normalization quantities for the quantum state.
!> Computes the trace of one-body and two-body RDMs, and orbital norms.
  subroutine PrinterObservable_DumpNorm(rdm1, rdm2, orbs, filename, toScreenQ)
    use M_Grid
    use M_Orbs

    !> One-body reduced density matrix in orbital representation.
    complex(R64), intent(in), contiguous :: rdm1(:, :)
    !> Two-body reduced density matrix in orbital representation.
    complex(R64), intent(in), contiguous :: rdm2(:, :, :, :)
    !> Orbital coefficients on the grid.
    complex(R64), intent(in), contiguous :: orbs(:, :)
    !> Output filename. If empty, no file output is produced.
    character(len=*), intent(in) :: filename
    !> Flag controlling whether to print output to screen.
    logical, intent(in) :: toScreenQ

    integer(I32) :: iGrid
    integer(I32) :: io, istat, nO, nOS, nG
    integer(I32) :: i1, i2
    character(256) :: msg

    complex(R64) :: normRdm1
    complex(R64) :: normRdm2
    complex(R64), allocatable :: normOrbs(:)

    nO = size(rdm1, 1)
    nG = Grid_nPoints
    nOS = Orbs_nOrbsInState

    allocate (normOrbs(nOS))

    normRdm1 = 0.0_R64
    do i1 = 1, nO
      normRdm1 = normRdm1 + rdm1(i1, i1)
    end do

    normRdm2 = 0.0_R64
    do i1 = 1, nO
      do i2 = 1, nO
        normRdm2 = normRdm2 + rdm2(i1, i2, i1, i2)
      end do
    end do

    normOrbs = 0.0_R64
    do i1 = 1, nOS
      do iGrid = 1, nG
        normOrbs(i1) = Grid_InnerProduct(orbs(:, i1), orbs(:, i1))
      end do
    end do

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Print to Screen
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (toScreenQ) then

      write (*, *)
      write (*, '(A, 2E20.10E3)') " Norm of 1RDM: ", normRdm1

      write (*, '(A, 2E20.10E3)') " Norm of 2RDM: ", normRdm2

      do i1 = 1, nOS
        write (*, '(A, I4, A, 2E20.10E3)') " Norm of Orbital ", i1, " is: ", normOrbs(i1)
      end do

    end if

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Print to file
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (filename .eq. "") return

    open (newunit=io, file=filename, status="unknown", position="append", iostat=istat, iomsg=msg)
    if (istat .ne. 0) print *, "Failed to open ", filename, "! iomsg = "//trim(msg)
    if (istat .ne. 0) error stop

    write (io, '(*(F10.6))') normRdm1, normRdm2, (normOrbs(i1), i1=1, nOS)

    close (io)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Calculates and dumps the orbital orthonormality matrix.
!> Computes the inner product between all pairs of orbitals to verify:
!> \[
!> \langle \phi_i | \phi_j \rangle = \delta_{ij}
!> \]
!> where \(\delta_{ij}\) is the Kronecker delta.
  subroutine PrinterObservable_DumpOrthonormality(orbs, filename, toScreenQ)
    use M_Grid
    use M_Orbs

    !> Orbital coefficients on the grid.
    complex(R64), intent(in), contiguous :: orbs(:, :)
    !> Output filename. If empty, no file output is produced.
    character(len=*), intent(in) :: filename
    !> Flag controlling whether to print output to screen.
    logical, intent(in) :: toScreenQ

    integer(I32) :: iGrid, i1, j1
    integer(I32) :: io, istat, nOS, nG
    character(256) :: msg

    complex(R64), allocatable :: innerP(:, :)

    nG = Grid_nPoints
    nOS = Orbs_nOrbsInState

    allocate (innerP(nOS, nOS))

    innerP(:, :) = 0.0_R64

    do i1 = 1, nOS
      do j1 = 1, nOS
        do iGrid = 1, nG
          innerP(i1, j1) = innerP(i1, j1) + Grid_InnerProduct(orbs(:, i1), orbs(:, j1))
        end do
      end do
    end do

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Print to Screen
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (toScreenQ) then

      do i1 = 1, nOS
        do j1 = 1, nOS
          write (*, '(A, I4, I4, A, 2E20.10E3)') " Inner product of orbital pair ", i1, j1, " is: ", innerP(i1, j1)
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

    do i1 = 1, nOS
      do j1 = 1, nOS
        write (io, '(I4, I4, 2E20.10E3)') i1, j1, innerP(i1, j1)
      end do
    end do

    close (io)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Calculates and dumps the orbital energy matrix.
!> Computes matrix elements of the system Hamiltonian in the orbital basis:
!> \[
!> E_{ij} = \langle \phi_i | \hat{H} | \phi_j \rangle
!> \]
!> where \(\hat{H}\) is the full Hamiltonian including kinetic, potential, and correlation terms.
  subroutine PrinterObservable_DumpOrbitalEnergyMatrix(rdm1, rdm2, orbs, filename, toScreenQ, time)
    use M_Grid
    use M_Orbs
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    !> Orbital coefficients on the grid.
    complex(R64), intent(in), contiguous :: orbs(:, :)
    !> One-body reduced density matrix in orbital representation.
    complex(R64), intent(in), contiguous :: rdm1(:, :)
    !> Two-body reduced density matrix in orbital representation.
    complex(R64), intent(in), contiguous :: rdm2(:, :, :, :)
    !> Output filename. If empty, no file output is produced.
    character(len=*), intent(in) :: filename
    !> Flag controlling whether to print output to screen.
    logical, intent(in) :: toScreenQ
    !> Current time value to be written with the energy data.
    real(R64), intent(in) :: time

    complex(R64), allocatable :: dOrbs(:, :)
    complex(R64), allocatable :: energyMatrix(:, :)

    integer(I32) :: nG, nOS

    integer(I32) :: io, istat, i, j
    character(256) :: msg

    nG = Grid_nPoints
    nOS = Orbs_nOrbsInState

    allocate (dOrbs(nG, nOS), energyMatrix(nOS, nOS))

    dOrbs(:, :) = 0.0_R64
    call Method_Mb_OrbBased_ApplySingleBodyOp(dOrbs, orbs, time)
    call Method_Mb_OrbBased_ApplyCorrelationOp(dOrbs, orbs, rdm1, rdm2, time)

    energyMatrix(:, :) = 0.0_R64
    do j = 1, nOS
      do i = 1, nOS
        energyMatrix(i, j) = energyMatrix(i, j) + Grid_InnerProduct(orbs(:, i), dOrbs(:, j))
      end do
    end do

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Print to Screen
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (toScreenQ) then

      write (*, *)
      write (*, *) "------------------------"
      write (*, *) "absolute of orbital energy matrix:  "
      write (*, *) "------------------------"
      write (*, *)

      do i = 1, nOS
        write (*, '(*(F10.6))') (abs(energyMatrix(i, j)), j=1, nOS)
      end do

    end if

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Print to file
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (filename .eq. "") return

    open (newunit=io, file=filename, status="unknown", position="append", iostat=istat, iomsg=msg)
    if (istat .ne. 0) print *, "Failed to open ", filename, "! iomsg = "//trim(msg)
    if (istat .ne. 0) error stop

    do i = 1, nOS
      write (io, '(*(F10.6))') (abs(energyMatrix(i, j)), j=1, nOS)
    end do

    write (io, *)
    close (io)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Calculates and dumps the total energy of the system.
!> Computes the total energy at a given time `t` by calling `Method_GetEnergy`.
!> \[
!> E = \langle \Psi | \hat{H} | \Psi \rangle
!> \]
  subroutine PrinterObservable_DumpEnergy(filename, toScreenQ, time)
    use M_Method
    use M_Orbs
    use M_Utils_UnusedVariables

    !> Output filename. If empty, no file output is produced.
    character(len=*), intent(in) :: filename
    !> Flag controlling whether to print output to screen.
    logical, intent(in) :: toScreenQ
    !> Current time value.
    real(R64), intent(in) :: time

    integer(I32) :: io, istat
    character(256) :: msg
    complex(R64) :: energy

    energy = Method_GetEnergy(time)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Print to Screen
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (toScreenQ) then
      write (*, *)
      write (*, '(A, 2E20.10E3)') " Energy: ", energy
    end if

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Print to file
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (filename .eq. "") return

    open (newunit=io, file=filename, status="unknown", position="append", iostat=istat, iomsg=msg)
    if (istat .ne. 0) print *, "Failed to open ", filename, "! iomsg = "//trim(msg)
    if (istat .ne. 0) error stop

    write (io, '(F10.6, 2E20.10E3)') time, energy

    close (io)

  end subroutine

end module
