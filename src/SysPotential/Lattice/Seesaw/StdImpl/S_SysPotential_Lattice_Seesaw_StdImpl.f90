! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysPotential_Lattice_Seesaw_StdImpl) S_SysPotential_Lattice_Seesaw_StdImpl

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysPotential_Lattice_Seesaw_StdImpl_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysPotential
    use M_SysPotential_Lattice

    call Say_Fabricate("cosinusLattice")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysPotential_FillExternalPotential => FillExternalPotential

  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillExternalPotential(externalPotential, time, bt_)
    use M_Utils_UnusedVariables
    use M_SysPotential_Lattice_Seesaw
    use M_Grid
    use M_Grid_Lattice
    use M_Utils_Constants, only: PI

    complex(R64), intent(out), allocatable :: externalPotential(:)
    real(R64), intent(in) :: time
    integer(I32), intent(in), optional :: bt_

    integer(I32) :: ix, iy, iz, iGrid
    real(R64) :: slopeX, slopeY, slopeZ
    real(R64) :: slopeMinX, slopeMinY, slopeMinZ
    real(R64) :: slopeMaxX, slopeMaxY, slopeMaxZ
    real(R64) :: frequencyX, frequencyY, frequencyZ
    real(R64) :: centerX, centerY, centerZ

    if (.false.) call UnusedVariables_Mark(bt_)

    if (.not. allocated(externalPotential)) allocate (externalPotential(Grid_nPoints))

    slopeMinX = SysPotential_Lattice_Seesaw_slopeMinX
    slopeMinY = SysPotential_Lattice_Seesaw_slopeMinY
    slopeMinZ = SysPotential_Lattice_Seesaw_slopeMinZ
    slopeMaxX = SysPotential_Lattice_Seesaw_slopeMaxX
    slopeMaxY = SysPotential_Lattice_Seesaw_slopeMaxY
    slopeMaxZ = SysPotential_Lattice_Seesaw_slopeMaxZ
    frequencyX = SysPotential_Lattice_Seesaw_frequencyX
    frequencyY = SysPotential_Lattice_Seesaw_frequencyY
    frequencyZ = SysPotential_Lattice_Seesaw_frequencyZ

    slopeX = 0.5_R64 * (slopeMaxX + slopeMinX) + 0.5_R64 * (slopeMaxX - slopeMinX) * sin(2.0_R64 * PI * frequencyX * time)
    slopeY = 0.5_R64 * (slopeMaxY + slopeMinY) + 0.5_R64 * (slopeMaxY - slopeMinY) * sin(2.0_R64 * PI * frequencyY * time)
    slopeZ = 0.5_R64 * (slopeMaxZ + slopeMinZ) + 0.5_R64 * (slopeMaxZ - slopeMinZ) * sin(2.0_R64 * PI * frequencyZ * time)

    centerX = 0.5_R64 * Grid_Lattice_xSize
    centerY = 0.5_R64 * Grid_Lattice_ySize
    centerZ = 0.5_R64 * Grid_Lattice_zSize

    do iz = 1, Grid_Lattice_zSize
      do iy = 1, Grid_Lattice_ySize
        do ix = 1, Grid_Lattice_xSize
          iGrid = Grid_Lattice_code(ix, iy, iz)
          externalPotential(iGrid) = slopeX * (ix - centerX) + &
                                     slopeY * (iy - centerY) + &
                                     slopeZ * (iz - centerZ)
        end do
      end do
    end do

  end subroutine

end submodule
