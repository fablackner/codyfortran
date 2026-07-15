! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Grid_Prolate) S_Grid_Prolate

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Grid_Prolate_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Grid
    use M_Grid_Prolate_Fedvr

    call Say_Fabricate("grid.prolate")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    Grid_Prolate_R = Json_Get("grid.prolate.R", 2.0_R64)
    Grid_Prolate_a = 0.5_R64 * Grid_Prolate_R
    Grid_Prolate_ximax = Json_Get("grid.prolate.ximax", 10.0_R64)
    Grid_Prolate_mmax = Json_Get("grid.prolate.mmax", 0)
    Grid_Prolate_nEta = Json_Get("grid.prolate.nEta", 8)

    Grid_InnerProduct => InnerProduct
    Grid_Prolate_SpatialProduct => SpatialProduct
    Grid_Prolate_GetMComponent => GetMComponent
    Grid_Prolate_SetMComponent => SetMComponent
    Grid_Prolate_AddMComponent => AddMComponent
    Grid_Setup => Setup

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("grid.prolate.fedvr")) then
      call Grid_Prolate_Fedvr_Fabricate

    else
      error stop "grid.prolate is missing one of: fedvr"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Setup
    use stdlib_quadrature, only: gauss_legendre
    use M_Utils_Say
    use M_Utils_DerivativeFedvr, only: LagrangeDerivative
    use M_Grid
    use M_Grid_Prolate

    integer(I32) :: iXi, jEta, kEta, iEta1, iEta2, mChannel, iGrid
    real(R64) :: jacobian

    call Say_Setup("grid.prolate")

    ! Call the specific setup procedure of the xi back-end
    ! (fills xiPoints, xiWeights, xiKinMatrix)
    call Grid_Prolate_Setup

    !------------------------------------
    ! eta direction: single-element Gauss-Legendre DVR (interior nodes only)
    !------------------------------------

    allocate (Grid_Prolate_etaPoints(Grid_Prolate_nEta))
    allocate (Grid_Prolate_etaWeights(Grid_Prolate_nEta))

    call gauss_legendre(Grid_Prolate_etaPoints, Grid_Prolate_etaWeights, [-1.0_R64, 1.0_R64])

    ! Sturm-Liouville matrix S_eta(i,j) = sum_k w_k (1-eta_k^2) D(k,i) D(k,j);
    ! the quadrature is exact (integrand degree 2*nEta-2 <= 2*nEta-1) and the
    ! boundary terms of the integration by parts vanish since (1-eta^2) = 0
    ! at eta = +-1
    allocate (Grid_Prolate_etaKinMatrix(Grid_Prolate_nEta, Grid_Prolate_nEta))
    Grid_Prolate_etaKinMatrix = 0.0_R64
    do iEta2 = 1, Grid_Prolate_nEta
      do iEta1 = 1, Grid_Prolate_nEta
        do kEta = 1, Grid_Prolate_nEta
          Grid_Prolate_etaKinMatrix(iEta1, iEta2) = Grid_Prolate_etaKinMatrix(iEta1, iEta2) + &
            Grid_Prolate_etaWeights(kEta) * (1.0_R64 - Grid_Prolate_etaPoints(kEta)**2) * &
            LagrangeDerivative(Grid_Prolate_etaPoints, iEta1, kEta, Grid_Prolate_nEta) * &
            LagrangeDerivative(Grid_Prolate_etaPoints, iEta2, kEta, Grid_Prolate_nEta)
        end do
      end do
    end do

    !------------------------------------
    ! flattened arrays (channel-outer layout, xi fastest)
    !------------------------------------

    Grid_Prolate_nSpatial = Grid_Prolate_nXi * Grid_Prolate_nEta

    allocate (Grid_Prolate_spatialWeights(Grid_Prolate_nSpatial))
    allocate (Grid_Prolate_weights(Grid_nPoints))
    allocate (Grid_Prolate_xiCoord(Grid_nPoints))
    allocate (Grid_Prolate_etaCoord(Grid_nPoints))
    allocate (Grid_Prolate_mCoord(Grid_nPoints))

    iGrid = 0
    do jEta = 1, Grid_Prolate_nEta
      do iXi = 1, Grid_Prolate_nXi
        iGrid = iGrid + 1

        ! Metric weight: quadrature x prolate-spheroidal Jacobian
        jacobian = Grid_Prolate_a**3 * (Grid_Prolate_xiPoints(iXi)**2 - Grid_Prolate_etaPoints(jEta)**2)
        Grid_Prolate_spatialWeights(iGrid) = Grid_Prolate_xiWeights(iXi) * Grid_Prolate_etaWeights(jEta) * jacobian
      end do
    end do

    ! Channels in interleaved order m = 0, -1, +1, -2, +2, ...
    do mChannel = 0, 2 * Grid_Prolate_mmax
      associate (base => mChannel * Grid_Prolate_nSpatial)
        Grid_Prolate_weights(base + 1:base + Grid_Prolate_nSpatial) = Grid_Prolate_spatialWeights(:)
        Grid_Prolate_mCoord(base + 1:base + Grid_Prolate_nSpatial) = MOfChannel(mChannel)

        iGrid = base
        do jEta = 1, Grid_Prolate_nEta
          do iXi = 1, Grid_Prolate_nXi
            iGrid = iGrid + 1
            Grid_Prolate_xiCoord(iGrid) = Grid_Prolate_xiPoints(iXi)
            Grid_Prolate_etaCoord(iGrid) = Grid_Prolate_etaPoints(jEta)
          end do
        end do
      end associate
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Position of the m channel in the interleaved order 0, -1, +1, -2, +2, ...
  !> (independent of the field's own mmax)
  pure function ChannelOfM(m) result(res)
    integer(I32), intent(in) :: m
    integer(I32) :: res

    if (m .eq. 0) then
      res = 0
    else if (m < 0) then
      res = 2 * abs(m) - 1
    else
      res = 2 * m
    end if

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Inverse of ChannelOfM: azimuthal quantum number of channel position
  pure function MOfChannel(channel) result(res)
    integer(I32), intent(in) :: channel
    integer(I32) :: res

    if (mod(channel, 2) .eq. 0) then
      res = channel / 2
    else
      res = -(channel + 1) / 2
    end if

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Computes the inner product of two fields
  pure function InnerProduct(fConjg, f) result(res)

    complex(R64), intent(in), contiguous, dimension(:) :: fConjg
    complex(R64), intent(in), contiguous, dimension(:) :: f
    complex(R64) :: res

    res = sum(conjg(fConjg(:)) * f(:) * Grid_Prolate_weights(:))

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Computes the spatial product of two fields (pointwise in (xi, eta),
  !> convolution in the azimuthal channels)
  pure subroutine SpatialProduct(fOut, fConjg, f, mMaxOut, mMaxConjg, mMax, conjgQ_, withWeightsQ_)
    use M_Utils_Constants

    complex(R64), intent(out), contiguous :: fOut(:)
    complex(R64), intent(in), contiguous  :: fConjg(:)
    complex(R64), intent(in), contiguous  :: f(:)
    integer(I32), intent(in)              :: mMaxOut
    integer(I32), intent(in)              :: mMaxConjg
    integer(I32), intent(in)              :: mMax
    logical, intent(in), optional         :: conjgQ_
    logical, intent(in), optional         :: withWeightsQ_

    integer(I32) :: m1, m2, mOut, mChannel
    logical :: conjgQ, withWeightsQ
    complex(R64), allocatable :: aTmp(:), bTmp(:), outTmp(:)

    conjgQ = .false.
    if (present(conjgQ_)) conjgQ = conjgQ_

    withWeightsQ = .false.
    if (present(withWeightsQ_)) withWeightsQ = withWeightsQ_

    allocate (aTmp(Grid_Prolate_nSpatial), bTmp(Grid_Prolate_nSpatial), outTmp(Grid_Prolate_nSpatial))

    fOut = (0.0_R64, 0.0_R64)

    do m1 = -mMaxConjg, mMaxConjg
      call GetMComponent(aTmp, m1, fConjg)
      if (all(abs(aTmp) < 1.0e-14_R64)) cycle
      if (conjgQ) aTmp = conjg(aTmp)

      do m2 = -mMax, mMax
        ! conjg(e^(i m1 phi)) e^(i m2 phi) = e^(i (m2-m1) phi)
        if (conjgQ) then
          mOut = m2 - m1
        else
          mOut = m2 + m1
        end if
        if (abs(mOut) > mMaxOut) cycle

        call GetMComponent(bTmp, m2, f)
        if (all(abs(bTmp) < 1.0e-14_R64)) cycle

        ! One factor 1/sqrt(2 pi) from the product of two channel basis
        ! functions e^(i m phi)/sqrt(2 pi)
        outTmp(:) = aTmp(:) * bTmp(:) / sqrt(TWOPI)

        call AddMComponent(fOut, mOut, outTmp)
      end do
    end do

    ! Apply weights if requested
    if (withWeightsQ) then
      do mChannel = 0, 2 * mMaxOut
        associate (base => mChannel * Grid_Prolate_nSpatial)
          fOut(base + 1:base + Grid_Prolate_nSpatial) = &
            fOut(base + 1:base + Grid_Prolate_nSpatial) * Grid_Prolate_spatialWeights(:)
        end associate
      end do
    end if

    deallocate (aTmp, bTmp, outTmp)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure subroutine GetMComponent(fMOut, m, f)

    complex(R64), intent(out) :: fMOut(:)
    integer(I32), intent(in) :: m
    complex(R64), intent(in) :: f(:)

    integer(I32) :: base

    base = ChannelOfM(m) * Grid_Prolate_nSpatial
    fMOut(1:Grid_Prolate_nSpatial) = f(base + 1:base + Grid_Prolate_nSpatial)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure subroutine SetMComponent(fOut, m, fM)

    complex(R64), intent(out) :: fOut(:)
    integer(I32), intent(in) :: m
    complex(R64), intent(in) :: fM(:)

    integer(I32) :: base

    base = ChannelOfM(m) * Grid_Prolate_nSpatial
    fOut(base + 1:base + Grid_Prolate_nSpatial) = fM(1:Grid_Prolate_nSpatial)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure subroutine AddMComponent(fInOut, m, fM)

    complex(R64), intent(inout) :: fInOut(:)
    integer(I32), intent(in) :: m
    complex(R64), intent(in) :: fM(:)

    integer(I32) :: base

    base = ChannelOfM(m) * Grid_Prolate_nSpatial
    fInOut(base + 1:base + Grid_Prolate_nSpatial) = &
      fInOut(base + 1:base + Grid_Prolate_nSpatial) + fM(1:Grid_Prolate_nSpatial)

  end subroutine

end submodule
