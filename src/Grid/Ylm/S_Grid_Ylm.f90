! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Grid_Ylm) S_Grid_Ylm

  implicit none

  !=============================================================================
  ! local data
  !=============================================================================

  !> Gauss-Legendre tables of Grid_Ylm_ToNodes/FromNodes and of the m = 0
  !> fast path of SpatialProduct (Setup): quadY(l, q, m) = Theta_lm(x_q) and
  !> quadWY(q, l, m) = 2 pi w_q Theta_lm(x_q) for 0 <= m <= Grid_Ylm_mmaxNodes,
  !> m <= l <= Grid_Ylm_lmaxNodes (zero for l < m) on the Grid_Ylm_nNodes
  !> Gauss-Legendre nodes x_q in cos(theta); Theta_l,-m = (-1)^m Theta_lm.
  real(R64), allocatable :: quadY(:, :, :)
  real(R64), allocatable :: quadWY(:, :, :)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Grid_Ylm_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Grid
    use M_Grid_Ylm_Const
    use M_Grid_Ylm_Fedvr
    use M_Grid_Ylm_FedvrEcs

    call Say_Fabricate("grid.ylm")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    ! Set maximum l value from JSON configuration
    Grid_Ylm_lmax = Json_Get("grid.ylm.lmax", 1)

    ! Set maximum |m| of the orbitals (default: full expansion, mmax = lmax)
    Grid_Ylm_mmax = Json_Get("grid.ylm.mmax", Grid_Ylm_lmax)
    if (Grid_Ylm_mmax < 0) error stop "grid.ylm.mmax must be >= 0"
    Grid_Ylm_mmax = min(Grid_Ylm_mmax, Grid_Ylm_lmax)

    ! Read radial grid parameters from JSON
    Grid_Ylm_rmax = Json_Get("grid.ylm.rmax", 20.0_R64)

    Grid_InnerProduct => InnerProduct
    Grid_Ylm_SpatialProduct => SpatialProduct
    Grid_Ylm_GetLmComponent => GetLmComponent
    Grid_Ylm_SetLmComponent => SetLmComponent
    Grid_Ylm_AddLmComponent => AddLmComponent
    Grid_Setup => Setup

    ! Hermitian metric by default; a c-product back-end overrides this below
    Grid_Ylm_cProductQ = .false.

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("grid.ylm.const")) then
      call Grid_Ylm_Const_Fabricate

    else if (Json_GetExistence("grid.ylm.fedvr")) then
      call Grid_Ylm_Fedvr_Fabricate

    else if (Json_GetExistence("grid.ylm.fedvrEcs")) then
      call Grid_Ylm_FedvrEcs_Fabricate

    else
      error stop "grid.ylm is missing one of: const, fedvr, fedvrEcs"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Setup
    use M_Utils_Say
    use M_Grid
    use M_Grid_Ylm

    integer(I32) :: iRad, iL, iM, iGrid

    call Say_Setup("grid.ylm")

    ! Call the specific setup procedure for the chosen grid type
    call Grid_Ylm_Setup

    ! Default Hermitian metric: real radial weights; c-product back-ends have
    ! already filled Grid_Ylm_radialMetricWeights in their own Setup
    if (.not. Grid_Ylm_cProductQ) then
      if (allocated(Grid_Ylm_radialMetricWeights)) deallocate (Grid_Ylm_radialMetricWeights)
      allocate (Grid_Ylm_radialMetricWeights(Grid_Ylm_nRadial))
      Grid_Ylm_radialMetricWeights(:) = cmplx(Grid_Ylm_radialWeights(:), 0.0_R64, kind=R64)
    end if

    ! Allocate and initialize weights for the whole grid
    allocate (Grid_Ylm_weights(Grid_nPoints))

    ! Allocate grid coordinate and indexing arrays
    allocate (Grid_Ylm_rCoord(Grid_nPoints))  ! Store (r, l, m) coordinates
    allocate (Grid_Ylm_lCoord(Grid_nPoints))  ! Store (r, l, m) coordinates
    allocate (Grid_Ylm_mCoord(Grid_nPoints))  ! Store (r, l, m) coordinates

    ! Set up grid structure - we need to handle conversion between (l,m) quantum numbers
    ! and array indices, considering that m can be negative
    iGrid = 0
    do iL = 0, Grid_Ylm_lmax
      do iM = -min(iL, Grid_Ylm_mmax), min(iL, Grid_Ylm_mmax)
        do iRad = 1, Grid_Ylm_nRadial
          iGrid = iGrid + 1

          ! Store actual quantum numbers and radius in Grid_Ylm_ylmcoord
          Grid_Ylm_rCoord(iGrid) = Grid_Ylm_radialPoints(iRad)  ! r value
          Grid_Ylm_lCoord(iGrid) = iL                           ! l value
          Grid_Ylm_mCoord(iGrid) = iM                           ! m value

          ! Calculate weights
          Grid_Ylm_weights(iGrid) = Grid_Ylm_radialWeights(iRad)
        end do
      end do
    end do

    call SetupQuadrature

  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Tabulate Theta_lm (Y_lm = Theta_lm(theta) e^(i m phi), Condon-Shortley
  !> phase) on Gauss-Legendre nodes for Grid_Ylm_ToNodes/FromNodes and the
  !> m = 0 fast path of SpatialProduct. The tables cover l <= 2 lmax + 2 and
  !> |m| <= 2 mmax (orbitals, densities and potentials of the default
  !> interaction layout, multiplicative operators). Products of fixed-m fields
  !> whose m add up to zero are polynomials in cos(theta) of degree l1 + l2 + ...
  !> (the sin(theta)^|m| factors pair up), so products of total degree
  !> <= 4 lmax + 5 are integrated exactly.
  subroutine SetupQuadrature
    use M_Utils_Constants
    use M_Grid_Ylm
    use stdlib_quadrature, only: gauss_legendre

    real(R64), allocatable :: x(:), w(:), sinTheta(:)
    integer(I32) :: l, m, lmaxNodes, mmaxNodes

    lmaxNodes = 2 * Grid_Ylm_lmax + 2
    mmaxNodes = min(2 * Grid_Ylm_mmax, lmaxNodes)
    Grid_Ylm_lmaxNodes = lmaxNodes
    Grid_Ylm_mmaxNodes = mmaxNodes
    Grid_Ylm_nNodes = 2 * Grid_Ylm_lmax + 3

    allocate (x(Grid_Ylm_nNodes), w(Grid_Ylm_nNodes))
    call gauss_legendre(x, w)
    sinTheta = sqrt(1.0_R64 - x(:)**2)

    if (allocated(quadY)) deallocate (quadY)
    if (allocated(quadWY)) deallocate (quadWY)
    allocate (quadY(0:lmaxNodes, Grid_Ylm_nNodes, 0:mmaxNodes), quadWY(Grid_Ylm_nNodes, 0:lmaxNodes, 0:mmaxNodes))
    quadY(:, :, :) = 0.0_R64

    ! Normalized associated Legendre functions by the standard stable
    ! recurrences: Theta_00 = 1/sqrt(4 pi),
    !   Theta_mm     = -sqrt((2m + 1)/(2m)) sin(theta) Theta_m-1,m-1,
    !   Theta_m+1,m  = sqrt(2m + 3) x Theta_mm,
    !   Theta_lm     = a_lm (x Theta_l-1,m - b_lm Theta_l-2,m),
    ! a_lm = sqrt((4l^2 - 1)/(l^2 - m^2)), b_lm = sqrt(((l-1)^2 - m^2)/(4(l-1)^2 - 1))
    quadY(0, :, 0) = 1.0_R64 / sqrt(4.0_R64 * PI)
    do m = 1, mmaxNodes
      quadY(m, :, m) = -sqrt((2 * m + 1) / (2.0_R64 * m)) * sinTheta(:) * quadY(m - 1, :, m - 1)
    end do
    do m = 0, mmaxNodes
      if (m + 1 <= lmaxNodes) quadY(m + 1, :, m) = sqrt(2.0_R64 * m + 3.0_R64) * x(:) * quadY(m, :, m)
      do l = m + 2, lmaxNodes
        quadY(l, :, m) = sqrt((4.0_R64 * l**2 - 1.0_R64) / real(l**2 - m**2, kind=R64)) * &
                         (x(:) * quadY(l - 1, :, m) - &
                          sqrt(real((l - 1)**2 - m**2, kind=R64) / (4.0_R64 * (l - 1)**2 - 1.0_R64)) * &
                          quadY(l - 2, :, m))
      end do
    end do

    do m = 0, mmaxNodes
      do l = 0, lmaxNodes
        quadWY(:, l, m) = 2.0_R64 * PI * w(:) * quadY(l, :, m)
      end do
    end do

    Grid_Ylm_nodeWeights = 2.0_R64 * PI * w(:)

  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Computes the inner product of two fields
  pure function InnerProduct(fConjg, f) result(res)

    complex(R64), intent(in), contiguous, dimension(:) :: fConjg
    complex(R64), intent(in), contiguous, dimension(:) :: f
    complex(R64) :: res

    res = sum(conjg(fConjg(:)) * f(:) * Grid_Ylm_weights(:))

  end function

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Computes the spatial product of two fields
  pure subroutine SpatialProduct(fOut, fConjg, f, conjgQ_, withWeightsQ_, &
                                 lMaxOut_, lMaxConjg_, lMax_, mMaxOut_, mMaxConjg_, mMax_)

    complex(R64), intent(out), contiguous :: fOut(:)
    complex(R64), intent(in), contiguous  :: fConjg(:)
    complex(R64), intent(in), contiguous  :: f(:)
    logical, intent(in), optional         :: conjgQ_
    logical, intent(in), optional         :: withWeightsQ_
    integer(I32), intent(in), optional    :: lMaxOut_
    integer(I32), intent(in), optional    :: lMaxConjg_
    integer(I32), intent(in), optional    :: lMax_
    integer(I32), intent(in), optional    :: mMaxOut_
    integer(I32), intent(in), optional    :: mMaxConjg_
    integer(I32), intent(in), optional    :: mMax_

    integer(I32) :: lMaxOut, lMaxConjg, lMax
    integer(I32) :: mMaxOut, mMaxConjg, mMax
    logical :: conjgQ, withWeightsQ

    if (.not. present(conjgQ_)) conjgQ = .false.
    if (present(conjgQ_)) conjgQ = conjgQ_

    if (.not. present(withWeightsQ_)) withWeightsQ = .false.
    if (present(withWeightsQ_)) withWeightsQ = withWeightsQ_

    ! Channel layouts default to the orbital layout
    lMaxOut = Grid_Ylm_lmax
    if (present(lMaxOut_)) lMaxOut = lMaxOut_

    lMaxConjg = Grid_Ylm_lmax
    if (present(lMaxConjg_)) lMaxConjg = lMaxConjg_

    lMax = Grid_Ylm_lmax
    if (present(lMax_)) lMax = lMax_

    mMaxOut = Grid_Ylm_mmax
    if (present(mMaxOut_)) mMaxOut = mMaxOut_

    mMaxConjg = Grid_Ylm_mmax
    if (present(mMaxConjg_)) mMaxConjg = mMaxConjg_

    mMax = Grid_Ylm_mmax
    if (present(mMax_)) mMax = mMax_

    if (lMaxConjg .eq. 0) then
      ! Spherically symmetric first factor (e.g. a central potential): the Gaunt
      ! sum collapses to Y_00 Y_lm = Y_lm / sqrt(4 pi), a radial scaling of f
      call SpatialProductL0(fOut, fConjg, f, conjgQ, lMaxOut, lMax, mMaxOut, mMax)

    else if (mMaxOut .eq. 0 .and. mMaxConjg .eq. 0 .and. mMax .eq. 0 .and. &
             max(lMaxOut, lMaxConjg, lMax) <= Grid_Ylm_lmaxNodes .and. &
             lMaxOut + lMaxConjg + lMax < 2 * Grid_Ylm_nNodes) then
      ! Axially symmetric operands: exact Gauss-Legendre quadrature in cos(theta)
      ! (O(lmax^2) per radial point instead of the O(lmax^3) Gaunt sum)
      call SpatialProductM0(fOut, fConjg, f, conjgQ, lMaxOut, lMaxConjg, lMax)

    else
      ! General case: Gaunt-coefficient expansion
      call SpatialProductGaunt(fOut, fConjg, f, conjgQ, lMaxOut, lMaxConjg, lMax, mMaxOut, mMaxConjg, mMax)
    end if

    ! Apply metric weights if requested (complex on c-product back-ends)
    if (withWeightsQ) call ApplyRadialMetricWeights(fOut, Grid_Ylm_NChannels(lMaxOut, mMaxOut))

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Multiplies the first nChannels radial blocks of f by the radial metric
  !> weights (complex on c-product back-ends).
  pure subroutine ApplyRadialMetricWeights(f, nChannels)
    complex(R64), intent(inout), contiguous :: f(:)
    integer(I32), intent(in)                :: nChannels

    integer(I32) :: nRad, iCh, i

    nRad = Grid_Ylm_nRadial

    do iCh = 1, nChannels
      i = (iCh - 1) * nRad + 1
      f(i:i + nRad - 1) = f(i:i + nRad - 1) * Grid_Ylm_radialMetricWeights(:)
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Spatial product of general fields by the Gaunt-coefficient expansion.
  !>
  !> Every product Y_{l1,m1} Y_{l2,m2} is expanded into the Y_{l,m} of the
  !> output layout; negligible channels of either factor are skipped.
  pure subroutine SpatialProductGaunt(fOut, fConjg, f, conjgQ, lMaxOut, lMaxConjg, lMax, mMaxOut, mMaxConjg, mMax)
    use M_Utils_SphericalHarmonics

    complex(R64), intent(out), contiguous :: fOut(:)
    complex(R64), intent(in), contiguous  :: fConjg(:)
    complex(R64), intent(in), contiguous  :: f(:)
    logical, intent(in)                   :: conjgQ
    integer(I32), intent(in)              :: lMaxOut
    integer(I32), intent(in)              :: lMaxConjg
    integer(I32), intent(in)              :: lMax
    integer(I32), intent(in)              :: mMaxOut
    integer(I32), intent(in)              :: mMaxConjg
    integer(I32), intent(in)              :: mMax

    integer(I32) :: l, m
    integer(I32) :: l1, m1, l2, m2
    integer(I32) :: nRad, iCh, i1, i2, iOut
    real(R64)    :: gVal
    logical, allocatable :: nonZeroQ(:)
    complex(R64), allocatable :: lmTmpConjg(:)

    nRad = Grid_Ylm_nRadial

    allocate (lmTmpConjg(nRad))

    ! Flag the non-negligible channels of f once; reused for every fConjg channel
    allocate (nonZeroQ(Grid_Ylm_NChannels(lMax, mMax)))
    do iCh = 1, size(nonZeroQ)
      i2 = (iCh - 1) * nRad + 1
      nonZeroQ(iCh) = .not. all(abs(f(i2:i2 + nRad - 1)) < 1.0e-14_R64)
    end do

    ! Initialize output array
    fOut = (0.0_R64, 0.0_R64)

    ! Build up expansion using Gaunt coefficients
    ! Loop over all (l₁, m₁) pairs
    do l1 = 0, lMaxConjg
      do m1 = -min(l1, mMaxConjg), min(l1, mMaxConjg)
        ! Complex conjugation: conjg(f)_{l1,m1} = (-1)^m1 conjg(f_{l1,-m1})
        if (conjgQ) then
          i1 = BaseIndex(l1, -m1, mMaxConjg)
        else
          i1 = BaseIndex(l1, m1, mMaxConjg)
        end if

        ! Check if this component is negligibly small
        if (all(abs(fConjg(i1:i1 + nRad - 1)) < 1.0e-14_R64)) cycle

        if (conjgQ) then
          lmTmpConjg(:) = conjg(fConjg(i1:i1 + nRad - 1))
          if (mod(abs(m1), 2) .ne. 0) lmTmpConjg(:) = -lmTmpConjg(:)
        else
          lmTmpConjg(:) = fConjg(i1:i1 + nRad - 1)
        end if

        ! Loop over all (l₂, m₂) pairs; channels of f are stored consecutively
        iCh = 0
        do l2 = 0, lMax
          do m2 = -min(l2, mMax), min(l2, mMax)
            iCh = iCh + 1
            if (.not. nonZeroQ(iCh)) cycle

            ! Calculate m value for output
            m = m1 + m2

            ! No output channel can hold this m value
            if (abs(m) > min(lMaxOut, mMaxOut)) cycle

            i2 = (iCh - 1) * nRad + 1

            ! Expand Y_{l1,m1}(Ω)·Y_{l2,m2}(Ω) into Y_{l,m}(Ω) using Gaunt coefficients
            do l = max(abs(l1 - l2), abs(m)), min(l1 + l2, lMaxOut)
              ! Gaunt coefficient vanishes for odd l1+l2+l (parity selection rule);
              ! skip before the (expensive) Wigner-3j evaluation
              if (mod(l1 + l2 + l, 2) .ne. 0) cycle

              ! Look up Gaunt coefficient (precomputed table with fallback)
              gVal = SphericalHarmonics_GauntTabulated(l1, m1, l2, m2, l)

              ! Skip if Gaunt coefficient is zero
              if (abs(gVal) < 1.0e-14_R64) cycle

              ! Add contribution directly to the (l,m) output channel
              iOut = BaseIndex(l, m, mMaxOut)
              fOut(iOut:iOut + nRad - 1) = fOut(iOut:iOut + nRad - 1) &
                                           + lmTmpConjg(:) * f(i2:i2 + nRad - 1) * gVal
            end do
          end do
        end do
      end do
    end do

    deallocate (lmTmpConjg, nonZeroQ)
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Spatial product with a spherically symmetric first factor g(r) Y_00.
  !>
  !> Since Y_00 Y_lm = Y_lm / sqrt(4 pi), every (l,m) channel of f that fits
  !> into the output layout is scaled by g(r) / sqrt(4 pi); the result equals
  !> the Gaunt-coefficient expansion exactly.
  pure subroutine SpatialProductL0(fOut, fConjg, f, conjgQ, lMaxOut, lMax, mMaxOut, mMax)
    use M_Utils_Constants, only: PI

    complex(R64), intent(out), contiguous :: fOut(:)
    complex(R64), intent(in), contiguous  :: fConjg(:)
    complex(R64), intent(in), contiguous  :: f(:)
    logical, intent(in)                   :: conjgQ
    integer(I32), intent(in)              :: lMaxOut
    integer(I32), intent(in)              :: lMax
    integer(I32), intent(in)              :: mMaxOut
    integer(I32), intent(in)              :: mMax

    integer(I32) :: nRad, l, m, iIn, iOut
    complex(R64), allocatable :: g(:)

    nRad = Grid_Ylm_nRadial
    allocate (g(nRad))

    ! conjg(Y_00) = Y_00, so only the radial factor is conjugated
    if (conjgQ) then
      g(:) = conjg(fConjg(1:nRad)) / sqrt(4.0_R64 * PI)
    else
      g(:) = fConjg(1:nRad) / sqrt(4.0_R64 * PI)
    end if

    fOut = (0.0_R64, 0.0_R64)

    do l = 0, min(lMax, lMaxOut)
      do m = -min(l, mMax, mMaxOut), min(l, mMax, mMaxOut)
        iIn = BaseIndex(l, m, mMax)
        iOut = BaseIndex(l, m, mMaxOut)
        fOut(iOut:iOut + nRad - 1) = g(:) * f(iIn:iIn + nRad - 1)
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Spatial product of axially symmetric (m = 0) fields by quadrature.
  !>
  !> With m = 0 the layout is a matrix f(r, l). Both factors are evaluated on
  !> the Gauss-Legendre nodes in cos(theta), multiplied pointwise and projected
  !> back onto Y_l0, l <= lMaxOut. The node count makes the angular integral
  !> exact, so the result equals the Gaunt-coefficient expansion.
  pure subroutine SpatialProductM0(fOut, fConjg, f, conjgQ, lMaxOut, lMaxConjg, lMax)
    complex(R64), intent(out), contiguous :: fOut(:)
    complex(R64), intent(in), contiguous  :: fConjg(:)
    complex(R64), intent(in), contiguous  :: f(:)
    logical, intent(in)                   :: conjgQ
    integer(I32), intent(in)              :: lMaxOut
    integer(I32), intent(in)              :: lMaxConjg
    integer(I32), intent(in)              :: lMax

    integer(I32) :: nRad
    real(R64), allocatable :: aRe(:, :), aIm(:, :), bRe(:, :), bIm(:, :), pRe(:, :), pIm(:, :)

    nRad = Grid_Ylm_nRadial

    fOut = (0.0_R64, 0.0_R64)

    ! First factor on the angular nodes (conjugation commutes with real Y_l0)
    aRe = matmul(reshape(real(fConjg(1:nRad * (lMaxConjg + 1))), [nRad, lMaxConjg + 1]), quadY(0:lMaxConjg, :, 0))
    aIm = matmul(reshape(aimag(fConjg(1:nRad * (lMaxConjg + 1))), [nRad, lMaxConjg + 1]), quadY(0:lMaxConjg, :, 0))
    if (conjgQ) aIm = -aIm

    ! Second factor on the angular nodes
    bRe = matmul(reshape(real(f(1:nRad * (lMax + 1))), [nRad, lMax + 1]), quadY(0:lMax, :, 0))
    bIm = matmul(reshape(aimag(f(1:nRad * (lMax + 1))), [nRad, lMax + 1]), quadY(0:lMax, :, 0))

    ! Pointwise product, then projection onto Y_l0
    pRe = aRe * bRe - aIm * bIm
    pIm = aRe * bIm + aIm * bRe

    fOut(1:nRad * (lMaxOut + 1)) = reshape(cmplx(matmul(pRe, quadWY(:, 0:lMaxOut, 0)), &
                                                 matmul(pIm, quadWY(:, 0:lMaxOut, 0)), kind=R64), &
                                           [nRad * (lMaxOut + 1)])

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Grid_Ylm_ToNodes(fNodes, f, lmax, m, mmax_)
    complex(R64), intent(out), contiguous :: fNodes(:, :)
    complex(R64), intent(in), contiguous  :: f(:)
    integer(I32), intent(in)              :: lmax
    integer(I32), intent(in)              :: m
    integer(I32), intent(in), optional    :: mmax_

    integer(I32) :: nRad, mmax, mAbs, l, i
    complex(R64), allocatable :: fM(:, :)

    mmax = Grid_Ylm_mmax
    if (present(mmax_)) mmax = mmax_
    mAbs = abs(m)

    if (mAbs > mmax) error stop "Grid_Ylm_ToNodes: |m| exceeds the layout"

    nRad = Grid_Ylm_nRadial

    ! Channels (l, m), l = |m| .. lmax, as a matrix fM(r, l)
    allocate (fM(nRad, max(0, lmax - mAbs + 1)))
    do l = mAbs, lmax
      i = BaseIndex(l, m, mmax)
      fM(:, l - mAbs + 1) = f(i:i + nRad - 1)
    end do

    call Grid_Ylm_ChannelsToNodes(fNodes, fM, m)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Grid_Ylm_FromNodes(f, fNodes, lmax, m, mmax_, withWeightsQ_)
    complex(R64), intent(out), contiguous :: f(:)
    complex(R64), intent(in), contiguous  :: fNodes(:, :)
    integer(I32), intent(in)              :: lmax
    integer(I32), intent(in)              :: m
    integer(I32), intent(in), optional    :: mmax_
    logical, intent(in), optional         :: withWeightsQ_

    integer(I32) :: nRad, mmax, mAbs, l, i
    complex(R64), allocatable :: fM(:, :)
    logical :: withWeightsQ

    if (.not. present(withWeightsQ_)) withWeightsQ = .false.
    if (present(withWeightsQ_)) withWeightsQ = withWeightsQ_

    mmax = Grid_Ylm_mmax
    if (present(mmax_)) mmax = mmax_
    mAbs = abs(m)

    if (mAbs > mmax) error stop "Grid_Ylm_FromNodes: |m| exceeds the layout"

    nRad = Grid_Ylm_nRadial

    f = (0.0_R64, 0.0_R64)
    if (mAbs > lmax) return

    allocate (fM(nRad, lmax - mAbs + 1))
    call Grid_Ylm_NodesToChannels(fM, fNodes, m)

    ! Metric weights if requested (complex on c-product back-ends)
    if (withWeightsQ) fM(:, :) = fM(:, :) * spread(Grid_Ylm_radialMetricWeights(:), 2, size(fM, 2))

    do l = mAbs, lmax
      i = BaseIndex(l, m, mmax)
      f(i:i + nRad - 1) = fM(:, l - mAbs + 1)
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> The real GEMMs act on the complex arrays viewed as real (2 nRad, ...)
  !> matrices: the tables are real, so real and imaginary parts transform
  !> alike and stay interleaved.
  module subroutine Grid_Ylm_ChannelsToNodes(fNodes, fM, m)
    use, intrinsic :: iso_c_binding, only: c_loc, c_f_pointer
    use M_Utils_BlasLib, only: DGEMM

    complex(R64), intent(out), contiguous, target :: fNodes(:, :)
    complex(R64), intent(in), contiguous, target  :: fM(:, :)
    integer(I32), intent(in)                      :: m

    real(R64), contiguous, pointer :: fNodesR(:, :), fMR(:, :)
    integer(I32) :: nRad2, nL, mAbs

    mAbs = abs(m)
    nL = size(fM, 2)
    nRad2 = 2 * size(fNodes, 1)

    if (mAbs + nL - 1 > Grid_Ylm_lmaxNodes) error stop "Grid_Ylm_ChannelsToNodes: lmax exceeds the quadrature tables"
    if (mAbs > Grid_Ylm_mmaxNodes) error stop "Grid_Ylm_ChannelsToNodes: |m| exceeds the quadrature tables"

    if (nL .eq. 0) then
      fNodes(:, :) = (0.0_R64, 0.0_R64)
      return
    end if

    call c_f_pointer(c_loc(fNodes), fNodesR, [nRad2, Grid_Ylm_nNodes])
    call c_f_pointer(c_loc(fM), fMR, [nRad2, nL])
    call DGEMM('N', 'N', nRad2, Grid_Ylm_nNodes, nL, 1.0_R64, fMR, nRad2, quadY(mAbs, 1, mAbs), &
               Grid_Ylm_lmaxNodes + 1, 0.0_R64, fNodesR, nRad2)
    if (m < 0 .and. mod(mAbs, 2) .eq. 1) fNodes(:, :) = -fNodes(:, :)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Grid_Ylm_NodesToChannels(fM, fNodes, m)
    use, intrinsic :: iso_c_binding, only: c_loc, c_f_pointer
    use M_Utils_BlasLib, only: DGEMM

    complex(R64), intent(out), contiguous, target :: fM(:, :)
    complex(R64), intent(in), contiguous, target  :: fNodes(:, :)
    integer(I32), intent(in)                      :: m

    real(R64), contiguous, pointer :: fNodesR(:, :), fMR(:, :)
    integer(I32) :: nRad2, nL, mAbs

    mAbs = abs(m)
    nL = size(fM, 2)
    nRad2 = 2 * size(fNodes, 1)

    if (mAbs + nL - 1 > Grid_Ylm_lmaxNodes) error stop "Grid_Ylm_NodesToChannels: lmax exceeds the quadrature tables"
    if (mAbs > Grid_Ylm_mmaxNodes) error stop "Grid_Ylm_NodesToChannels: |m| exceeds the quadrature tables"

    if (nL .eq. 0) return

    call c_f_pointer(c_loc(fNodes), fNodesR, [nRad2, Grid_Ylm_nNodes])
    call c_f_pointer(c_loc(fM), fMR, [nRad2, nL])
    call DGEMM('N', 'N', nRad2, nL, Grid_Ylm_nNodes, 1.0_R64, fNodesR, nRad2, quadWY(1, mAbs, mAbs), &
               Grid_Ylm_nNodes, 0.0_R64, fMR, nRad2)
    if (m < 0 .and. mod(mAbs, 2) .eq. 1) fM(:, :) = -fM(:, :)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure module function Grid_Ylm_NChannels(lmax, mmax) result(res)
    use M_Utils_SphericalHarmonics

    integer(I32), intent(in) :: lmax
    integer(I32), intent(in) :: mmax
    integer(I32) :: res

    res = SphericalHarmonics_NChannels(lmax, mmax)

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> First flattened index of channel (l,m); channels are l-major with m
  !> running over -min(l,mmax)..min(l,mmax), each holding nRadial points.
  pure function BaseIndex(l, m, mmax) result(res)
    use M_Utils_SphericalHarmonics

    integer(I32), intent(in) :: l
    integer(I32), intent(in) :: m
    integer(I32), intent(in) :: mmax
    integer(I32) :: res

    res = 1 + (SphericalHarmonics_ChannelIndex(l, m, mmax) - 1) * Grid_Ylm_nRadial

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure subroutine GetLmComponent(fLmOut, l, m, f, mmax_)

    complex(R64), intent(out) :: fLmOut(:)  ! Output radial function for (l,m)
    integer(I32), intent(in) :: l         ! Angular momentum quantum number
    integer(I32), intent(in) :: m         ! Magnetic quantum number
    complex(R64), intent(in) :: f(:)      ! Full function on grid
    integer(I32), intent(in), optional :: mmax_

    integer(I32) :: baseIdx, mmax

    mmax = Grid_Ylm_mmax
    if (present(mmax_)) mmax = mmax_

    baseIdx = BaseIndex(l, m, mmax)
    fLmOut(1:Grid_Ylm_nRadial) = f(baseIdx:baseIdx + Grid_Ylm_nRadial - 1)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure subroutine SetLmComponent(fOut, l, m, fLm, mmax_)

    complex(R64), intent(inout) :: fOut(:)   ! Full function on grid (target)
    integer(I32), intent(in) :: l         ! Angular momentum quantum number
    integer(I32), intent(in) :: m         ! Magnetic quantum number
    complex(R64), intent(in) :: fLm(:)   ! Input radial function for (l,m)
    integer(I32), intent(in), optional :: mmax_

    integer(I32) :: baseIdx, mmax

    mmax = Grid_Ylm_mmax
    if (present(mmax_)) mmax = mmax_

    baseIdx = BaseIndex(l, m, mmax)
    fOut(baseIdx:baseIdx + Grid_Ylm_nRadial - 1) = fLm(1:Grid_Ylm_nRadial)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure subroutine AddLmComponent(fInOut, l, m, fLm, mmax_)

    complex(R64), intent(inout) :: fInOut(:)   ! Full function on grid (target)
    integer(I32), intent(in) :: l         ! Angular momentum quantum number
    integer(I32), intent(in) :: m         ! Magnetic quantum number
    complex(R64), intent(in) :: fLm(:)   ! Input radial function for (l,m)
    integer(I32), intent(in), optional :: mmax_

    integer(I32) :: baseIdx, mmax

    mmax = Grid_Ylm_mmax
    if (present(mmax_)) mmax = mmax_

    baseIdx = BaseIndex(l, m, mmax)
    fInOut(baseIdx:baseIdx + Grid_Ylm_nRadial - 1) = fInOut(baseIdx:baseIdx + Grid_Ylm_nRadial - 1) &
                                                     + fLm(1:Grid_Ylm_nRadial)

  end subroutine

end submodule
