! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_OrbsInit_Lattice) S_OrbsInit_Lattice

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine OrbsInit_Lattice_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_OrbsInit
    use M_OrbsInit_Lattice_OnSite

    call Say_Fabricate("orbsInit.lattice")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    OrbsInit_InitializeOrb => InitializeOrb

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("orbsInit.lattice.onSite")) then
      call OrbsInit_Lattice_OnSite_Fabricate

    else
      error stop "orbsInit.lattice is missing one of: onSite"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine InitializeOrb(orb, ind, bt_)
    use M_Grid
    use M_Grid_Lattice

    complex(R64), intent(out), contiguous :: orb(:)
    integer(I32), intent(in)              :: ind
    integer(I32), intent(in), optional    :: bt_

    integer(I32) :: iGrid
    real(R64)    :: norm
    integer(I32) :: ix, iy, iz

    ! Initialize to zero
    orb(:) = 0.0_R64

    ! Calculate orbital values for each grid point
    do iz = 1, Grid_Lattice_zSize
      do iy = 1, Grid_Lattice_ySize
        do ix = 1, Grid_Lattice_xSize
          iGrid = Grid_Lattice_code(ix, iy, iz)
          orb(iGrid) = OrbsInit_Lattice_InitFunction(ix, iy, iz, ind, bt_)
        end do
      end do
    end do

    ! Now verify normalization directly with Grid_InnerProduct
    norm = real(Grid_InnerProduct(orb, orb), kind=R64)

    ! Apply normalization correction
    orb(:) = orb(:) / sqrt(norm)

  end subroutine

end submodule
