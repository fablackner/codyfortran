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
    logical :: ecsActiveQ

    call Say_Fabricate("grid.ylm.fedvrEcs")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    ! Read FedvrEcs-specific parameters
    Grid_Ylm_FedvrEcs_nElements = Json_Get("grid.ylm.fedvrEcs.nElements", 10)
    Grid_Ylm_FedvrEcs_nLocals = Json_Get("grid.ylm.fedvrEcs.nLocals", 7)

    defaultStart = Grid_Ylm_FedvrEcs_nElements + 1
    Grid_Ylm_FedvrEcs_firstEcsElement = Json_Get("grid.ylm.fedvrEcs.firstEcsElement", defaultStart)
    Grid_Ylm_FedvrEcs_transitionElements = Json_Get("grid.ylm.fedvrEcs.transitionElements", 0)

    thetaDeg = Json_Get("grid.ylm.fedvrEcs.thetaDeg", 0.0_R64)
    pi = acos(-1.0_R64)
    Grid_Ylm_FedvrEcs_theta = thetaDeg * (pi / 180.0_R64)

    ! Calculate total grid points from FedvrEcs parameters (without rmin and rmax)
    Grid_Ylm_nRadial = Grid_Ylm_FedvrEcs_nElements * (Grid_Ylm_FedvrEcs_nLocals - 1)

    ! Update total grid size to include the radial component
    Grid_nPoints = Grid_Ylm_nRadial * (Grid_Ylm_lmax + 1)**2

    ! Set the setup procedure
    Grid_Ylm_Setup => Setup

    ! With an actively rotated contour the Hamiltonian is complex-symmetric,
    ! not Hermitian: the consistent metric is the c-product (no conjugation)
    ! along the contour, so declare it and override the default Ylm inner
    ! product; Setup fills Grid_Ylm_radialMetricWeights accordingly
    ecsActiveQ = (abs(Grid_Ylm_FedvrEcs_theta) > 0.0_R64) .and. &
                 (Grid_Ylm_FedvrEcs_firstEcsElement <= Grid_Ylm_FedvrEcs_nElements)
    Grid_Ylm_cProductQ = ecsActiveQ
    if (ecsActiveQ) Grid_InnerProduct => InnerProduct

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

    integer(I32) :: iRad, iChannel, iGrid
    real(R64) :: rReal
    complex(R64) :: zContour, wMass, wVolume

    call Say_Setup("grid.ylm.fedvrEcs")

    ! Allocate arrays for the radial grid
    allocate (Grid_Ylm_radialPoints(Grid_Ylm_nRadial))
    allocate (Grid_Ylm_radialWeights(Grid_Ylm_nRadial))

    if (allocated(Grid_Ylm_FedvrEcs_contourPoints)) deallocate (Grid_Ylm_FedvrEcs_contourPoints)
    if (allocated(Grid_Ylm_FedvrEcs_massWeights)) deallocate (Grid_Ylm_FedvrEcs_massWeights)
    if (allocated(Grid_Ylm_FedvrEcs_volumeWeights)) deallocate (Grid_Ylm_FedvrEcs_volumeWeights)
    if (allocated(Grid_Ylm_FedvrEcs_weights)) deallocate (Grid_Ylm_FedvrEcs_weights)

    allocate (Grid_Ylm_FedvrEcs_contourPoints(Grid_Ylm_nRadial))
    allocate (Grid_Ylm_FedvrEcs_massWeights(Grid_Ylm_nRadial))
    allocate (Grid_Ylm_FedvrEcs_volumeWeights(Grid_Ylm_nRadial))
    allocate (Grid_Ylm_FedvrEcs_weights(Grid_nPoints))

    ! Initialize the FedvrEcs grid using the utility module
    call FedvrEcs_CreateCtx(Grid_Ylm_FedvrEcs_fedvrEcsCtx, &
                            0.0_R64, &
                            Grid_Ylm_rmax, &
                            Grid_Ylm_FedvrEcs_nElements, &
                            Grid_Ylm_FedvrEcs_nLocals, &
                            .true., &
                            .false., &
                            ecsStartElement_=Grid_Ylm_FedvrEcs_firstEcsElement, &
                            ecsAngle_=Grid_Ylm_FedvrEcs_theta, &
                            ecsTransitionElements_=Grid_Ylm_FedvrEcs_transitionElements)

    Grid_Ylm_FedvrEcs_firstEcsElement = Grid_Ylm_FedvrEcs_fedvrEcsCtx % firstEcsElement
    Grid_Ylm_FedvrEcs_transitionElements = Grid_Ylm_FedvrEcs_fedvrEcsCtx % transitionElements
    Grid_Ylm_FedvrEcs_theta = Grid_Ylm_FedvrEcs_fedvrEcsCtx % ecsAngle
    Grid_Ylm_FedvrEcs_ecsRadius = Grid_Ylm_FedvrEcs_fedvrEcsCtx % ecsRadius

    call DerivativeFedvrEcs_CreateCtx(Grid_Ylm_FedvrEcs_derivativeCtx, Grid_Ylm_FedvrEcs_fedvrEcsCtx)

    do iRad = 1, Grid_Ylm_nRadial
      zContour = Grid_Ylm_FedvrEcs_fedvrEcsCtx % points(iRad)
      rReal = Grid_Ylm_FedvrEcs_fedvrEcsCtx % pointsReal(iRad)
      wMass = Grid_Ylm_FedvrEcs_fedvrEcsCtx % weights(iRad)
      wVolume = wMass * zContour * zContour

      Grid_Ylm_FedvrEcs_contourPoints(iRad) = zContour
      Grid_Ylm_FedvrEcs_massWeights(iRad) = wMass
      Grid_Ylm_FedvrEcs_volumeWeights(iRad) = wVolume

      Grid_Ylm_radialPoints(iRad) = rReal
      Grid_Ylm_radialWeights(iRad) = real(wMass) * rReal * rReal
    end do

    ! Repeat the complex volume weights over all (l,m) channels; the flattened
    ! layout is (l,m) outer, radial inner (see S_Grid_Ylm Setup)
    iGrid = 0
    do iChannel = 1, (Grid_Ylm_lmax + 1)**2
      do iRad = 1, Grid_Ylm_nRadial
        iGrid = iGrid + 1
        Grid_Ylm_FedvrEcs_weights(iGrid) = Grid_Ylm_FedvrEcs_volumeWeights(iRad)
      end do
    end do

    ! On an active contour the radial metric is the complex volume weights;
    ! otherwise the generic Ylm Setup builds the default Hermitian metric
    if (Grid_Ylm_cProductQ) then
      if (allocated(Grid_Ylm_radialMetricWeights)) deallocate (Grid_Ylm_radialMetricWeights)
      allocate (Grid_Ylm_radialMetricWeights(Grid_Ylm_nRadial))
      Grid_Ylm_radialMetricWeights(:) = Grid_Ylm_FedvrEcs_volumeWeights(:)
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Computes the c-product of two fields along the ECS contour.
  !>
  !> On a rotated contour the Hamiltonian is complex-symmetric; orthogonality
  !> and matrix elements are defined by the bilinear c-product with complex
  !> volume weights, without conjugation of the bra.
  pure function InnerProduct(fConjg, f) result(res)

    complex(R64), intent(in), contiguous, dimension(:) :: fConjg
    complex(R64), intent(in), contiguous, dimension(:) :: f
    complex(R64) :: res

    res = sum(fConjg(:) * f(:) * Grid_Ylm_FedvrEcs_weights(:))

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end submodule
