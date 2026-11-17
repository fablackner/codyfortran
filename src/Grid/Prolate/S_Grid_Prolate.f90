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
    real(R64), allocatable :: dEta(:, :), gOdd(:, :)

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

    ! Sturm-Liouville matrices S_eta(i,j,parity) = sum_k w_k (1-eta_k^2) chi_i'(eta_k) chi_j'(eta_k);
    ! the boundary terms of the integration by parts vanish since (1-eta^2) = 0
    ! at eta = +-1.
    !   parity 0 (even |m|): chi_i = L_i; the quadrature is exact (integrand
    !     degree 2*nEta-2 <= 2*nEta-1).
    !   parity 1 (odd |m|): chi_i = s(eta) L_i / s(eta_i) with s = sqrt(1-eta^2),
    !     the axis factor of odd-m fields; s chi_i' = g_i / s_i at the nodes with
    !     g_i(eta_k) = -eta_k delta_ki + (1-eta_k^2) L_i'(eta_k) (DVR quadrature).
    allocate (Grid_Prolate_etaKinMatrix(Grid_Prolate_nEta, Grid_Prolate_nEta, 0:1))
    allocate (dEta(Grid_Prolate_nEta, Grid_Prolate_nEta), gOdd(Grid_Prolate_nEta, Grid_Prolate_nEta))
    do iEta1 = 1, Grid_Prolate_nEta
      do kEta = 1, Grid_Prolate_nEta
        dEta(kEta, iEta1) = LagrangeDerivative(Grid_Prolate_etaPoints, iEta1, kEta, Grid_Prolate_nEta)
        gOdd(kEta, iEta1) = (1.0_R64 - Grid_Prolate_etaPoints(kEta)**2) * dEta(kEta, iEta1)
        if (kEta .eq. iEta1) gOdd(kEta, iEta1) = gOdd(kEta, iEta1) - Grid_Prolate_etaPoints(kEta)
        gOdd(kEta, iEta1) = gOdd(kEta, iEta1) / sqrt(1.0_R64 - Grid_Prolate_etaPoints(iEta1)**2)
      end do
    end do
    Grid_Prolate_etaKinMatrix = 0.0_R64
    do iEta2 = 1, Grid_Prolate_nEta
      do iEta1 = 1, Grid_Prolate_nEta
        do kEta = 1, Grid_Prolate_nEta
          Grid_Prolate_etaKinMatrix(iEta1, iEta2, 0) = &
            Grid_Prolate_etaKinMatrix(iEta1, iEta2, 0) + &
            Grid_Prolate_etaWeights(kEta) * (1.0_R64 - Grid_Prolate_etaPoints(kEta)**2) * &
            dEta(kEta, iEta1) * dEta(kEta, iEta2)
          Grid_Prolate_etaKinMatrix(iEta1, iEta2, 1) = Grid_Prolate_etaKinMatrix(iEta1, iEta2, 1) + &
                                                       Grid_Prolate_etaWeights(kEta) * gOdd(kEta, iEta1) * gOdd(kEta, iEta2)
        end do
      end do
    end do
    deallocate (dEta, gOdd)

    !------------------------------------
    ! flattened arrays (channel-outer layout, xi fastest)
    !------------------------------------

    Grid_Prolate_nSpatial = Grid_Prolate_nXi * Grid_Prolate_nEta

    allocate (Grid_Prolate_spatialWeights(Grid_Prolate_nSpatial))
    allocate (Grid_weights(Grid_nPoints))
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
        Grid_weights(base + 1:base + Grid_Prolate_nSpatial) = Grid_Prolate_spatialWeights(:)
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

    integer(I32) :: m1, m2, mOut, mChannel, base1, base2, baseOut, nSpatial
    logical :: conjgQ, withWeightsQ
    logical, allocatable :: nonzeroQ(:)
    complex(R64), allocatable :: aTmp(:)

    conjgQ = .false.
    if (present(conjgQ_)) conjgQ = conjgQ_

    withWeightsQ = .false.
    if (present(withWeightsQ_)) withWeightsQ = withWeightsQ_

    nSpatial = Grid_Prolate_nSpatial
    allocate (aTmp(nSpatial))

    ! channels of f with a non-negligible value (checked once, not per m1;
    ! real and imaginary parts separately, which avoids a complex abs)
    allocate (nonzeroQ(-mMax:mMax))
    do m2 = -mMax, mMax
      base2 = ChannelOfM(m2) * nSpatial
      nonzeroQ(m2) = .not. NegligibleQ(f(base2 + 1:base2 + nSpatial))
    end do

    fOut = (0.0_R64, 0.0_R64)

    do m1 = -mMaxConjg, mMaxConjg
      base1 = ChannelOfM(m1) * nSpatial
      if (NegligibleQ(fConjg(base1 + 1:base1 + nSpatial))) cycle

      ! One factor 1/sqrt(2 pi) from the product of two channel basis
      ! functions e^(i m phi)/sqrt(2 pi)
      if (conjgQ) then
        aTmp(:) = conjg(fConjg(base1 + 1:base1 + nSpatial)) / sqrt(TWOPI)
      else
        aTmp(:) = fConjg(base1 + 1:base1 + nSpatial) / sqrt(TWOPI)
      end if

      do m2 = -mMax, mMax
        ! conjg(e^(i m1 phi)) e^(i m2 phi) = e^(i (m2-m1) phi)
        if (conjgQ) then
          mOut = m2 - m1
        else
          mOut = m2 + m1
        end if
        if (abs(mOut) > mMaxOut) cycle
        if (.not. nonzeroQ(m2)) cycle

        base2 = ChannelOfM(m2) * nSpatial
        baseOut = ChannelOfM(mOut) * nSpatial
        fOut(baseOut + 1:baseOut + nSpatial) = fOut(baseOut + 1:baseOut + nSpatial) + &
                                               aTmp(:) * f(base2 + 1:base2 + nSpatial)
      end do
    end do

    ! Channel Dirichlet point: an m /= 0 field vanishes on the axis (xi = 1),
    ! as in the kinetic operator and the Poisson solver. Without it a
    ! potential times orbital product keeps the axis value of the m /= 0
    ! channels, which the kinetic operator does not couple: every such point
    ! is then a spurious eigenvector with eigenvalue V(1, eta) (about -1e3 Ha
    ! near the nuclei of N2)
    do mChannel = 1, 2 * mMaxOut
      associate (base => mChannel * Grid_Prolate_nSpatial)
        fOut(base + 1:base + Grid_Prolate_nSpatial:Grid_Prolate_nXi) = (0.0_R64, 0.0_R64)
      end associate
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

    deallocate (aTmp, nonzeroQ)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> True if every real and imaginary part of a channel is below 1e-14.
  pure function NegligibleQ(fM) result(res)
    complex(R64), intent(in) :: fM(:)
    logical :: res

    res = all(abs(real(fM, kind=R64)) < 1.0e-14_R64) .and. all(abs(aimag(fM)) < 1.0e-14_R64)

  end function

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

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Grid_Prolate_FillNodeCoordinate(coord, direction)
    real(R64), allocatable, intent(out) :: coord(:)
    real(R64), intent(in)               :: direction(3)

    if (.not. allocated(Grid_Prolate_xiCoord)) error stop "Grid_Prolate_FillNodeCoordinate: call Grid_Setup first"
    if (any(abs(direction(1:2)) > 0.0_R64)) &
      error stop "Grid_Prolate_FillNodeCoordinate: only along z (x and y couple the m channels)"

    coord = direction(3) * Grid_Prolate_a * Grid_Prolate_xiCoord(:) * Grid_Prolate_etaCoord(:)

  end subroutine

end submodule
