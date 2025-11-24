! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysKinetic_Lattice_NearestNeighbor) S_SysKinetic_Lattice_NearestNeighbor

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysKinetic_Lattice_NearestNeighbor_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysKinetic

    call Say_Fabricate("sysKinetic.lattice.nearestneighbor")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysKinetic_Lattice_NearestNeighbor_hoppX = Json_Get("sysKinetic.lattice.nearestNeighbor.hoppX", 1.0_R64)
    SysKinetic_Lattice_NearestNeighbor_hoppY = Json_Get("sysKinetic.lattice.nearestNeighbor.hoppY", 1.0_R64)
    SysKinetic_Lattice_NearestNeighbor_hoppZ = Json_Get("sysKinetic.lattice.nearestNeighbor.hoppZ", 1.0_R64)

    SysKinetic_timeIndependentQ = .true.
    SysKinetic_bodyTypeIndependentQ = .true.

    SysKinetic_MultiplyWithKineticOp => MultiplyWithKineticOp

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine MultiplyWithKineticOp(dOrb, orb, time, bt_)
    use M_Utils_UnusedVariables
    use M_Grid
    use M_Grid_Lattice

    complex(R64), intent(out), contiguous   :: dOrb(:)
    complex(R64), intent(in), contiguous    :: orb(:)
    real(R64), intent(in)                   :: time
    integer(I32), intent(in), optional      :: bt_

    integer(I32) :: rltbfb
    integer(I32) :: ix, iy, iz
    integer(I32) :: tix, tiy, tiz
    integer(I32) :: pos, tpos

    if (.false.) call UnusedVariables_Mark(time, bt_)

    dOrb(:) = 0.0_R64

    ! hopping is defined with negative sign!

    Do rltbfb = 1, 6 ! sum over all possible directions to jump (right left top bottom front back)

      Do iz = 1, Grid_Lattice_zSize
        Do iy = 1, Grid_Lattice_ySize
          Do ix = 1, Grid_Lattice_xSize

            select case (rltbfb)

            case (1)! hop to the right
              if (.not. Grid_Lattice_xPeriodicQ .and. ix .eq. Grid_Lattice_xSize) cycle ! for hard walls no jump

              tix = modulo(ix, Grid_Lattice_xSize) + 1
              tiy = iy
              tiz = iz

              pos = Grid_Lattice_code(ix, iy, iz)
              tpos = Grid_Lattice_code(tix, tiy, tiz)

              dOrb(tpos) = dOrb(tpos) - SysKinetic_Lattice_NearestNeighbor_hoppX * orb(pos)

            case (2) ! hop to the left
              if (.not. Grid_Lattice_xPeriodicQ .and. ix .eq. 1) cycle ! for hard walls no jump

              tix = modulo(ix - 2, Grid_Lattice_xSize) + 1
              tiy = iy
              tiz = iz

              pos = Grid_Lattice_code(ix, iy, iz)
              tpos = Grid_Lattice_code(tix, tiy, tiz)

              dOrb(tpos) = dOrb(tpos) - SysKinetic_Lattice_NearestNeighbor_hoppX * orb(pos)

            case (3) ! hop to the top
              if (.not. Grid_Lattice_yPeriodicQ .and. iy .eq. Grid_Lattice_ySize) cycle ! for hard walls no jump

              tix = ix
              tiy = modulo(iy, Grid_Lattice_ySize) + 1
              tiz = iz

              pos = Grid_Lattice_code(ix, iy, iz)
              tpos = Grid_Lattice_code(tix, tiy, tiz)

              dOrb(tpos) = dOrb(tpos) - SysKinetic_Lattice_NearestNeighbor_hoppY * orb(pos)

            case (4) ! hop to the bottom
              if (.not. Grid_Lattice_yPeriodicQ .and. iy .eq. 1) cycle ! for hard walls no jump

              tix = ix
              tiy = modulo(iy - 2, Grid_Lattice_ySize) + 1
              tiz = iz

              pos = Grid_Lattice_code(ix, iy, iz)
              tpos = Grid_Lattice_code(tix, tiy, tiz)

              dOrb(tpos) = dOrb(tpos) - SysKinetic_Lattice_NearestNeighbor_hoppY * orb(pos)

            case (5)! hop to the front
              if (.not. Grid_Lattice_zPeriodicQ .and. iz .eq. Grid_Lattice_zSize) cycle ! for hard walls no jump

              tix = ix
              tiy = iy
              tiz = modulo(iz, Grid_Lattice_zSize) + 1

              pos = Grid_Lattice_code(ix, iy, iz)
              tpos = Grid_Lattice_code(tix, tiy, tiz)

              dOrb(tpos) = dOrb(tpos) - SysKinetic_Lattice_NearestNeighbor_hoppZ * orb(pos)

            case (6) ! hop to the back
              if (.not. Grid_Lattice_zPeriodicQ .and. iz .eq. 1) cycle ! for hard walls no jump

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
