! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file S_SysKinetic_Lattice_NearestNeighbor.f90
!> @brief Implementation of nearest-neighbor tight-binding kinetic operator.
submodule(M_SysKinetic_Lattice_NearestNeighbor) S_SysKinetic_Lattice_NearestNeighbor

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Read hopping parameters from JSON and bind the kinetic operator.
  module subroutine SysKinetic_Lattice_NearestNeighbor_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysKinetic

    call Say_Fabricate("sysKinetic.lattice.nearestneighbor")

    !------------------------------------
    ! read configuration values
    !------------------------------------

    SysKinetic_Lattice_NearestNeighbor_hoppX = Json_Get("sysKinetic.lattice.nearestNeighbor.hoppX", 1.0_R64)
    SysKinetic_Lattice_NearestNeighbor_hoppY = Json_Get("sysKinetic.lattice.nearestNeighbor.hoppY", 1.0_R64)
    SysKinetic_Lattice_NearestNeighbor_hoppZ = Json_Get("sysKinetic.lattice.nearestNeighbor.hoppZ", 1.0_R64)

    !------------------------------------
    ! set global flags and bind operator
    !------------------------------------

    SysKinetic_timeIndependentQ = .true.
    SysKinetic_bodyTypeIndependentQ = .true.

    SysKinetic_MultiplyWithKineticOp => MultiplyWithKineticOp

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Apply T̂ψ via nearest-neighbor hopping on the 3D lattice.
  !>
  !> @details
  !> Implements the tight-binding kinetic operator by summing hopping terms
  !> from all 6 neighboring sites. For each lattice site (ix, iy, iz), the
  !> operator adds contributions:
  !>
  !>     dOrb(neighbor) ← dOrb(neighbor) − t_d × orb(site)
  !>
  !> where d ∈ {±x, ±y, ±z} and t_d is the corresponding hopping amplitude.
  !> Hard-wall boundaries skip hops that would exit the lattice; periodic
  !> boundaries wrap via modulo arithmetic.
  !>
  !> Complexity: O(6 × nSites) = O(N) per orbital application.
  subroutine MultiplyWithKineticOp(dOrb, orb, time, bt_)
    use M_Utils_UnusedVariables
    use M_Grid
    use M_Grid_Lattice

    complex(R64), intent(out), contiguous   :: dOrb(:)
    complex(R64), intent(in), contiguous    :: orb(:)
    real(R64), intent(in)                   :: time
    integer(I32), intent(in), optional      :: bt_

    integer(I32) :: iDir
    integer(I32) :: ix, iy, iz
    integer(I32) :: tix, tiy, tiz
    integer(I32) :: pos, tpos

    if (.false.) call UnusedVariables_Mark(time, bt_)

    dOrb(:) = 0.0_R64

    ! Loop over all 6 hopping directions (right, left, top, bottom, front, back)
    Do iDir = 1, 6

      Do iz = 1, Grid_Lattice_zSize
        Do iy = 1, Grid_Lattice_ySize
          Do ix = 1, Grid_Lattice_xSize

            select case (iDir)

            case (1) ! hop +x (right)
              if (.not. Grid_Lattice_xPeriodicQ .and. ix .eq. Grid_Lattice_xSize) cycle

              tix = modulo(ix, Grid_Lattice_xSize) + 1
              tiy = iy
              tiz = iz

              pos = Grid_Lattice_code(ix, iy, iz)
              tpos = Grid_Lattice_code(tix, tiy, tiz)

              dOrb(tpos) = dOrb(tpos) - SysKinetic_Lattice_NearestNeighbor_hoppX * orb(pos)

            case (2) ! hop −x (left)
              if (.not. Grid_Lattice_xPeriodicQ .and. ix .eq. 1) cycle

              tix = modulo(ix - 2, Grid_Lattice_xSize) + 1
              tiy = iy
              tiz = iz

              pos = Grid_Lattice_code(ix, iy, iz)
              tpos = Grid_Lattice_code(tix, tiy, tiz)

              dOrb(tpos) = dOrb(tpos) - SysKinetic_Lattice_NearestNeighbor_hoppX * orb(pos)

            case (3) ! hop +y (top)
              if (.not. Grid_Lattice_yPeriodicQ .and. iy .eq. Grid_Lattice_ySize) cycle

              tix = ix
              tiy = modulo(iy, Grid_Lattice_ySize) + 1
              tiz = iz

              pos = Grid_Lattice_code(ix, iy, iz)
              tpos = Grid_Lattice_code(tix, tiy, tiz)

              dOrb(tpos) = dOrb(tpos) - SysKinetic_Lattice_NearestNeighbor_hoppY * orb(pos)

            case (4) ! hop −y (bottom)
              if (.not. Grid_Lattice_yPeriodicQ .and. iy .eq. 1) cycle

              tix = ix
              tiy = modulo(iy - 2, Grid_Lattice_ySize) + 1
              tiz = iz

              pos = Grid_Lattice_code(ix, iy, iz)
              tpos = Grid_Lattice_code(tix, tiy, tiz)

              dOrb(tpos) = dOrb(tpos) - SysKinetic_Lattice_NearestNeighbor_hoppY * orb(pos)

            case (5) ! hop +z (front)
              if (.not. Grid_Lattice_zPeriodicQ .and. iz .eq. Grid_Lattice_zSize) cycle

              tix = ix
              tiy = iy
              tiz = modulo(iz, Grid_Lattice_zSize) + 1

              pos = Grid_Lattice_code(ix, iy, iz)
              tpos = Grid_Lattice_code(tix, tiy, tiz)

              dOrb(tpos) = dOrb(tpos) - SysKinetic_Lattice_NearestNeighbor_hoppZ * orb(pos)

            case (6) ! hop −z (back)
              if (.not. Grid_Lattice_zPeriodicQ .and. iz .eq. 1) cycle

              tix = ix
              tiy = iy
              tiz = modulo(iz - 2, Grid_Lattice_zSize) + 1

              pos = Grid_Lattice_code(ix, iy, iz)
              tpos = Grid_Lattice_code(tix, tiy, tiz)

              dOrb(tpos) = dOrb(tpos) - SysKinetic_Lattice_NearestNeighbor_hoppZ * orb(pos)

            end select

          end do
        end do
      end do

    end do

  end subroutine

end submodule
