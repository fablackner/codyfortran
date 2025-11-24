! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Grid_Ylm_Fedvr) S_Grid_Ylm_Fedvr

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Grid_Ylm_Fedvr_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Grid
    use M_Grid_Ylm

    call Say_Fabricate("grid.ylm.fedvr")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    ! Read FEDVR-specific parameters
    Grid_Ylm_Fedvr_nElements = Json_Get("grid.ylm.fedvr.nElements", 10)
    Grid_Ylm_Fedvr_nLocals = Json_Get("grid.ylm.fedvr.nLocals", 7)

    ! Calculate total grid points from FEDVR parameters (without rmin and rmax)
    Grid_Ylm_nRadial = Grid_Ylm_Fedvr_nElements * (Grid_Ylm_Fedvr_nLocals - 1)

    ! Update total grid size to include the radial component
    Grid_nPoints = Grid_Ylm_nRadial * (Grid_Ylm_lmax + 1)**2

    ! Set the setup procedure
    Grid_Ylm_Setup => Setup

    !------------------------------------
    ! branch
    !------------------------------------

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Setup
    use M_Utils_Say
    use M_Utils_Fedvr
    use M_Utils_DerivativeFedvr
    use M_Grid
    use M_Grid_Ylm

    integer(I32) :: iRad
    real(R64) :: r, w

    call Say_Setup("grid.ylm.fedvr")

    ! Allocate arrays for the radial grid
    allocate (Grid_Ylm_radialPoints(Grid_Ylm_nRadial))
    allocate (Grid_Ylm_radialWeights(Grid_Ylm_nRadial))

    ! Initialize the FEDVR grid using the utility module
    call Fedvr_CreateCtx(Grid_Ylm_Fedvr_fedvrCtx, &
                         0.0_R64, &
                         Grid_Ylm_rmax, &
                         Grid_Ylm_Fedvr_nElements, &
                         Grid_Ylm_Fedvr_nLocals, &
                         .true., &
                         .false.)

    call DerivativeFedvr_CreateCtx(Grid_Ylm_Fedvr_derivativeCtx, Grid_Ylm_Fedvr_fedvrCtx)

    do iRad = 1, Grid_Ylm_nRadial
      r = Grid_Ylm_Fedvr_fedvrCtx % points(iRad) ! exclude r=0 and rmax
      w = Grid_Ylm_Fedvr_fedvrCtx % weights(iRad) ! exclude r=0 and rmax

      Grid_Ylm_radialPoints(iRad) = r
      Grid_Ylm_radialWeights(iRad) = w * r * r
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end submodule
