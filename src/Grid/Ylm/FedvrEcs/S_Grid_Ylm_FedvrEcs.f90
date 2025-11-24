! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Grid_Ylm_FedvrEcs) S_Grid_Ylm_FedvrEcs

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Grid_Ylm_FedvrEcs_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Grid
    use M_Grid_Ylm

    real(R64) :: thetaDeg, pi

    integer(I32) :: defaultStart

    call Say_Fabricate("grid.ylm.FedvrEcs")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    ! Read FedvrEcs-specific parameters
    Grid_Ylm_FedvrEcs_nElements = Json_Get("grid.ylm.FedvrEcs.nElements", 10)
    Grid_Ylm_FedvrEcs_nLocals = Json_Get("grid.ylm.FedvrEcs.nLocals", 7)

    defaultStart = Grid_Ylm_FedvrEcs_nElements + 1
    Grid_Ylm_FedvrEcs_firstEcsElement = Json_Get("grid.ylm.FedvrEcs.firstEcsElement", defaultStart)
    Grid_Ylm_FedvrEcs_transitionElements = Json_Get("grid.ylm.FedvrEcs.transitionElements", 0)

    thetaDeg = Json_Get("grid.ylm.FedvrEcs.thetaDeg", 0.0_R64)
    pi = acos(-1.0_R64)
    Grid_Ylm_FedvrEcs_theta = thetaDeg * (pi / 180.0_R64)

    ! Calculate total grid points from FedvrEcs parameters (without rmin and rmax)
    Grid_Ylm_nRadial = Grid_Ylm_FedvrEcs_nElements * (Grid_Ylm_FedvrEcs_nLocals - 1)

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
    use M_Utils_FedvrEcs
    use M_Utils_DerivativeFedvrEcs
    use M_Grid
    use M_Grid_Ylm

    integer(I32) :: iRad
    real(R64) :: rReal
    complex(R64) :: zContour, wMass, wVolume

    call Say_Setup("grid.ylm.FedvrEcs")

    ! Allocate arrays for the radial grid
    allocate (Grid_Ylm_radialPoints(Grid_Ylm_nRadial))
    allocate (Grid_Ylm_radialWeights(Grid_Ylm_nRadial))

    if (allocated(Grid_Ylm_FedvrEcs_contourPoints)) deallocate (Grid_Ylm_FedvrEcs_contourPoints)
    if (allocated(Grid_Ylm_FedvrEcs_massWeights)) deallocate (Grid_Ylm_FedvrEcs_massWeights)
    if (allocated(Grid_Ylm_FedvrEcs_volumeWeights)) deallocate (Grid_Ylm_FedvrEcs_volumeWeights)

    allocate (Grid_Ylm_FedvrEcs_contourPoints(Grid_Ylm_nRadial))
    allocate (Grid_Ylm_FedvrEcs_massWeights(Grid_Ylm_nRadial))
    allocate (Grid_Ylm_FedvrEcs_volumeWeights(Grid_Ylm_nRadial))

    ! Initialize the FedvrEcs grid using the utility module
    call FedvrEcs_CreateCtx(Grid_Ylm_FedvrEcs_FedvrEcsCtx, &
                            0.0_R64, &
                            Grid_Ylm_rmax, &
                            Grid_Ylm_FedvrEcs_nElements, &
                            Grid_Ylm_FedvrEcs_nLocals, &
                            .true., &
                            .false., &
                            Grid_Ylm_FedvrEcs_firstEcsElement, &
                            Grid_Ylm_FedvrEcs_theta, &
                            Grid_Ylm_FedvrEcs_transitionElements)

    Grid_Ylm_FedvrEcs_firstEcsElement = Grid_Ylm_FedvrEcs_FedvrEcsCtx % firstEcsElement
    Grid_Ylm_FedvrEcs_transitionElements = Grid_Ylm_FedvrEcs_FedvrEcsCtx % transitionElements
    Grid_Ylm_FedvrEcs_theta = Grid_Ylm_FedvrEcs_FedvrEcsCtx % ecsAngle
    Grid_Ylm_FedvrEcs_ecsRadius = Grid_Ylm_FedvrEcs_FedvrEcsCtx % ecsRadius

    call DerivativeFedvrEcs_CreateCtx(Grid_Ylm_FedvrEcs_derivativeCtx, Grid_Ylm_FedvrEcs_FedvrEcsCtx)

    do iRad = 1, Grid_Ylm_nRadial
      zContour = Grid_Ylm_FedvrEcs_FedvrEcsCtx % points(iRad)
      rReal = Grid_Ylm_FedvrEcs_FedvrEcsCtx % pointsReal(iRad)
      wMass = Grid_Ylm_FedvrEcs_FedvrEcsCtx % weights(iRad)
      wVolume = wMass * zContour * zContour

      Grid_Ylm_FedvrEcs_contourPoints(iRad) = zContour
      Grid_Ylm_FedvrEcs_massWeights(iRad) = wMass
      Grid_Ylm_FedvrEcs_volumeWeights(iRad) = wVolume

      Grid_Ylm_radialPoints(iRad) = rReal
      Grid_Ylm_radialWeights(iRad) = real(wMass) * rReal * rReal
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end submodule
