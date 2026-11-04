! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Module M_Utils_PrinterObservableFem provides observables of the state on the
!> 3D finite-element grid (grid.fem).
!>
!> Fem-geometry specialization: the one-body density
!>   rho(r) = sum_ij rho_s(i, j) phi_j^*(r) phi_i(r)
!> of the spatial 1-RDM rho_s and the orbitals, integrated with the nodal
!> quadrature weights Grid_Fem_weights (the grid metric of Grid_InnerProduct),
!> gives
!>   N = int rho = Tr(rho_s S),  S_ji = <phi_j|phi_i>   (electron number),
!>   <r_a> = int r_a rho = sum_ij rho_s(i, j) <phi_j| r_a |phi_i>  (dipole,
!>           a = x, y, z, positions of the nodes relative to the origin),
!> and the dipole velocity is v_a = i Tr(rho_s [T, r_a]) with the kinetic
!> operator T of SysKinetic. None assumes orthonormal orbitals (absorber,
!> Cayley tolerance).
!>
!> Method-aware (as M_Utils_PrinterObservable): rdm1(i, j) = <a+_j a_i> is the
!> 1-RDM of Method over Method_Mb_OrbBased_nOrbsSum spin orbitals, orbs holds
!> the Orbs_nOrbsInState stored orbitals, and rho_s is
!> Method_Mb_OrbBased_FillSpatialRdm1 of rdm1: with restricted storage
!> (nOrbsSum = 2 nOrbsInState) the sum of the two diagonal spin blocks (for a
!> closed shell twice the first block, the occupation factor nO/nOS of the
!> Linear printer), otherwise rdm1. Stored orbital i has the body type
!> Method_Mb_OrbBased_bodyTypeOfOrb(i).
module M_Utils_PrinterObservableFem
  use M_Utils_Types

  implicit none

  private
  public :: PrinterObservableFem_DumpDipole, PrinterObservableFem_DumpNorm, PrinterObservableFem_DumpDipoleVelocity

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Calculates and dumps the dipole moment (<x>, <y>, <z>) of the state:
!> \[
!> \langle \hat{r}_a \rangle = \sum_{i,j} \rho_{ij} \langle \phi_j | r_a | \phi_i \rangle
!> \]
!> One line per call appended to filename: "time <x> <y> <z>" (4 columns).
  subroutine PrinterObservableFem_DumpDipole(rdm1, orbs, filename, toScreenQ, time, dipole_)
    !> One-body reduced density matrix in orbital representation.
    complex(R64), intent(in), contiguous :: rdm1(:, :)
    !> Orbital coefficients on the grid (Grid_nPoints x Orbs_nOrbsInState).
    complex(R64), intent(in), contiguous :: orbs(:, :)
    !> Output filename. If empty, no file output is produced.
    character(len=*), intent(in) :: filename
    !> Flag controlling whether to print output to screen.
    logical, intent(in) :: toScreenQ
    !> Current time value to be written with the dipole data.
    real(R64), intent(in) :: time
    !> Optional output of the dipole moment (<x>, <y>, <z>).
    real(R64), intent(out), optional :: dipole_(3)

    real(R64) :: moments(0:3)

    if (.not. (toScreenQ .or. (filename .ne. "") .or. present(dipole_))) return

    call FillDensityMoments(moments, rdm1, orbs)

    if (present(dipole_)) dipole_(:) = moments(1:3)

    call Dump(" Dipole moment <x>, <y>, <z>: ", moments(1:3), filename, toScreenQ, time)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Calculates and dumps the particle number N = Tr(rho_s S) of the state (the
!> integral of the one-body density; equals the number of bodies for
!> orthonormal orbitals and a normalized many-body part).
!> One line per call appended to filename: "time N" (2 columns).
  subroutine PrinterObservableFem_DumpNorm(rdm1, orbs, filename, toScreenQ, time, norm_)
    !> One-body reduced density matrix in orbital representation.
    complex(R64), intent(in), contiguous :: rdm1(:, :)
    !> Orbital coefficients on the grid (Grid_nPoints x Orbs_nOrbsInState).
    complex(R64), intent(in), contiguous :: orbs(:, :)
    !> Output filename. If empty, no file output is produced.
    character(len=*), intent(in) :: filename
    !> Flag controlling whether to print output to screen.
    logical, intent(in) :: toScreenQ
    !> Current time value to be written with the norm.
    real(R64), intent(in) :: time
    !> Optional output of the particle number N.
    real(R64), intent(out), optional :: norm_

    real(R64) :: moments(0:3)

    if (.not. (toScreenQ .or. (filename .ne. "") .or. present(norm_))) return

    call FillDensityMoments(moments, rdm1, orbs)

    if (present(norm_)) norm_ = moments(0)

    call Dump(" Particle number Tr(rho S): ", moments(0:0), filename, toScreenQ, time)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Calculates and dumps the dipole velocity (v_x, v_y, v_z), the velocity form
!> of the dipole on the grid:
!> \[
!> v_a = i\, \mathrm{Tr}(\rho [\hat{T}, r_a])
!>   = \mathrm{Re}\, i \sum_{i,j} \rho_{ij} \left( \langle T\phi_j | r_a \phi_i \rangle
!>     - \langle r_a \phi_j | T\phi_i \rangle \right)
!> \]
!> with the kinetic operator T of SysKinetic (mass of the body type of each
!> orbital). Potential, length gauge and interaction commute with the total r_a,
!> so v_a = d<r_a>/dt wherever the equations of motion obey Ehrenfest's theorem
!> for r_a (exact dynamics, TDHF). After a kick exp(i kappa r_a) of every orbital
!> it is N kappa/m up to the discretization of T (TRK sum rule on the grid).
!> One line per call appended to filename: "time v_x v_y v_z" (4 columns).
  subroutine PrinterObservableFem_DumpDipoleVelocity(rdm1, orbs, filename, toScreenQ, time, velocity_)
    !> One-body reduced density matrix in orbital representation.
    complex(R64), intent(in), contiguous :: rdm1(:, :)
    !> Orbital coefficients on the grid (Grid_nPoints x Orbs_nOrbsInState).
    complex(R64), intent(in), contiguous :: orbs(:, :)
    !> Output filename. If empty, no file output is produced.
    character(len=*), intent(in) :: filename
    !> Flag controlling whether to print output to screen.
    logical, intent(in) :: toScreenQ
    !> Current time value to be written with the data (and passed to the
    !> kinetic operator).
    real(R64), intent(in) :: time
    !> Optional output of the dipole velocity (d<x>/dt, d<y>/dt, d<z>/dt).
    real(R64), intent(out), optional :: velocity_(3)

    real(R64) :: velocity(3)

    if (.not. (toScreenQ .or. (filename .ne. "") .or. present(velocity_))) return

    call FillDipoleVelocity(velocity, rdm1, orbs, time)

    if (present(velocity_)) velocity_(:) = velocity(:)

    call Dump(" Dipole velocity d<x>/dt, d<y>/dt, d<z>/dt: ", velocity, filename, toScreenQ, time)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Moments of the one-body density: moments(0) = N = int rho,
!> moments(a) = int r_a rho (a = 1, 2, 3 for x, y, z), with the nodal weights.
!> One pass over the nodes, O(nOS^2) per node; the node coordinates come from
!> Grid_Fem_FillNodeCoordinate (three real arrays of Grid_nPoints).
  subroutine FillDensityMoments(moments, rdm1, orbs)
    !$ use omp_lib, only: omp_in_parallel
    use M_Grid_Fem
    use M_Orbs
    use M_Method_Mb_OrbBased, only: Method_Mb_OrbBased_FillSpatialRdm1

    real(R64), intent(out) :: moments(0:3)
    complex(R64), intent(in), contiguous :: rdm1(:, :)
    complex(R64), intent(in), contiguous :: orbs(:, :)

    complex(R64), allocatable :: rdm1Spatial(:, :)
    real(R64), allocatable :: xNodes(:), yNodes(:), zNodes(:)
    complex(R64) :: density
    real(R64) :: weightedDensity, m0, m1, m2, m3
    integer(I32) :: nOS, ix, iy, iz, node, i, j

    call CheckArguments(rdm1, orbs)

    nOS = Orbs_nOrbsInState
    call Method_Mb_OrbBased_FillSpatialRdm1(rdm1Spatial, rdm1)

    call Grid_Fem_FillNodeCoordinate(xNodes, [1.0_R64, 0.0_R64, 0.0_R64])
    call Grid_Fem_FillNodeCoordinate(yNodes, [0.0_R64, 1.0_R64, 0.0_R64])
    call Grid_Fem_FillNodeCoordinate(zNodes, [0.0_R64, 0.0_R64, 1.0_R64])

    m0 = 0.0_R64
    m1 = 0.0_R64
    m2 = 0.0_R64
    m3 = 0.0_R64

    ! Nodes in storage order, OpenMP over the z planes (static schedule)
    !$omp parallel do default(shared) private(ix, iy, iz, node, i, j, density, weightedDensity) &
    !$omp reduction(+:m0, m1, m2, m3) schedule(static) if (.not. omp_in_parallel())
    do iz = 1, Grid_Fem_nz
      do iy = 1, Grid_Fem_ny
        do ix = 1, Grid_Fem_nx
          node = Grid_Fem_Index(ix, iy, iz)

          ! rho = sum_ij rho_s(i, j) conjg(phi_j) phi_i (real for Hermitian rho_s)
          density = (0.0_R64, 0.0_R64)
          do j = 1, nOS
            do i = 1, nOS
              density = density + rdm1Spatial(i, j) * conjg(orbs(node, j)) * orbs(node, i)
            end do
          end do
          weightedDensity = Grid_Fem_weights(node) * real(density, kind=R64)

          m0 = m0 + weightedDensity
          m1 = m1 + weightedDensity * xNodes(node)
          m2 = m2 + weightedDensity * yNodes(node)
          m3 = m3 + weightedDensity * zNodes(node)
        end do
      end do
    end do
    !$omp end parallel do

    moments(:) = [m0, m1, m2, m3]

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> velocity(a) = Re(i sum_ij rho_s(i, j) (<T phi_j|r_a phi_i> - <r_a phi_j|T phi_i>))
!> for a = 1, 2, 3 (x, y, z): T phi_i once per orbital (SysKinetic at time,
!> body type Method_Mb_OrbBased_bodyTypeOfOrb(i)), r_a phi_i per axis (node
!> coordinates of Grid_Fem_FillNodeCoordinate), 2 nOS^2 inner products per
!> axis.
  subroutine FillDipoleVelocity(velocity, rdm1, orbs, time)
    use M_Utils_Constants, only: IU
    use M_Grid
    use M_Grid_Fem
    use M_SysKinetic
    use M_Orbs
    use M_Method_Mb_OrbBased, only: Method_Mb_OrbBased_FillSpatialRdm1, Method_Mb_OrbBased_bodyTypeOfOrb

    real(R64), intent(out) :: velocity(3)
    complex(R64), intent(in), contiguous :: rdm1(:, :)
    complex(R64), intent(in), contiguous :: orbs(:, :)
    real(R64), intent(in) :: time

    complex(R64), allocatable :: rdm1Spatial(:, :), rOrbs(:, :), tOrbs(:, :)
    real(R64), allocatable :: coordinate(:)
    real(R64) :: direction(3)
    complex(R64) :: trace
    integer(I32) :: nOS, iAxis, i, j

    call CheckArguments(rdm1, orbs)

    nOS = Orbs_nOrbsInState
    call Method_Mb_OrbBased_FillSpatialRdm1(rdm1Spatial, rdm1)

    allocate (rOrbs(Grid_nPoints, nOS), tOrbs(Grid_nPoints, nOS))
    do i = 1, nOS
      call SysKinetic_MultiplyWithKineticOp(tOrbs(:, i), orbs(:, i), time, bt_=Method_Mb_OrbBased_bodyTypeOfOrb(i))
    end do

    do iAxis = 1, 3
      direction(:) = 0.0_R64
      direction(iAxis) = 1.0_R64
      call Grid_Fem_FillNodeCoordinate(coordinate, direction)
      do i = 1, nOS
        rOrbs(:, i) = coordinate(:) * orbs(:, i)
      end do

      trace = (0.0_R64, 0.0_R64)
      do j = 1, nOS
        do i = 1, nOS
          trace = trace + rdm1Spatial(i, j) * (Grid_InnerProduct(tOrbs(:, j), rOrbs(:, i)) &
                                               - Grid_InnerProduct(rOrbs(:, j), tOrbs(:, i)))
        end do
      end do
      velocity(iAxis) = real(IU * trace, kind=R64)
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Stops unless the grid is grid.fem (after Grid_Setup), orbs holds the
!> stored orbitals and rdm1 is the 1-RDM of Method.
  subroutine CheckArguments(rdm1, orbs)
    use M_Grid
    use M_Grid_Fem
    use M_Orbs
    use M_Method_Mb_OrbBased, only: Method_Mb_OrbBased_nOrbsSum

    complex(R64), intent(in), contiguous :: rdm1(:, :)
    complex(R64), intent(in), contiguous :: orbs(:, :)

    if (.not. allocated(Grid_Fem_weights)) &
      error stop "PrinterObservableFem: requires grid.fem (after Grid_Setup)"
    if (size(orbs, 1) .ne. Grid_nPoints .or. size(orbs, 2) < Orbs_nOrbsInState) &
      error stop "PrinterObservableFem: orbs must hold Orbs_nOrbsInState orbitals of Grid_nPoints nodes"
    if (Orbs_nOrbsInState < 1 .or. size(rdm1, 1) .ne. Method_Mb_OrbBased_nOrbsSum &
        .or. size(rdm1, 2) .ne. Method_Mb_OrbBased_nOrbsSum) &
      error stop "PrinterObservableFem: rdm1 must be the nOrbsSum x nOrbsSum 1-RDM of Method_Mb_OrbBased_FillRdm1"

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Print (time, values) to screen and/or append them to a file as one line.
  subroutine Dump(label, values, filename, toScreenQ, time)
    character(len=*), intent(in) :: label
    real(R64), intent(in) :: values(:)
    character(len=*), intent(in) :: filename
    logical, intent(in) :: toScreenQ
    real(R64), intent(in) :: time

    integer(I32) :: io, istat
    character(256) :: msg

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Print to Screen
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (toScreenQ) then
      write (*, '(A, *(ES24.15E3))') label, values
    end if

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Print to file
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (filename .eq. "") return

    open (newunit=io, file=filename, status="unknown", position="append", iostat=istat, iomsg=msg)
    if (istat .ne. 0) print *, "Failed to open ", filename, "! iomsg = "//trim(msg)
    if (istat .ne. 0) error stop

    write (io, '(*(ES24.15E3))') time, values

    close (io)

  end subroutine

end module
