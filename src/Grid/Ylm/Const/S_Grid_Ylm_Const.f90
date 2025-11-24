! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Grid_Ylm_Const) S_Grid_Ylm_Const

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Grid_Ylm_Const_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Grid
    use M_Grid_Ylm

    call Say_Fabricate("grid.ylm.const")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    ! Read radial grid parameters from JSON
    Grid_Ylm_nRadial = Json_Get("grid.ylm.const.nRad", 100)

    ! Update total grid size to include the radial component
    Grid_nPoints = Grid_Ylm_nRadial * (Grid_Ylm_lmax + 1)**2

    ! Set the procedure pointer for radial grid setup
    Grid_Ylm_Setup => Setup

    !------------------------------------
    ! branch
    !------------------------------------

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Setup
    use M_Utils_Say
    use M_Grid
    use M_Grid_Ylm

    integer(I32) :: iRad
    real(R64) :: dr, r

    call Say_Setup("grid.ylm.const")

    ! Allocate arrays for the radial grid
    allocate (Grid_Ylm_radialPoints(Grid_Ylm_nRadial))
    allocate (Grid_Ylm_radialWeights(Grid_Ylm_nRadial))

    ! Set up the radial grid with constant spacing
    dr = Grid_Ylm_rmax / (Grid_Ylm_nRadial + 1)  ! +1 to avoid r=0 and r=rmax

    do iRad = 1, Grid_Ylm_nRadial
      r = iRad * dr
      Grid_Ylm_radialPoints(iRad) = r
      Grid_Ylm_radialWeights(iRad) = dr * r * r
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end submodule
