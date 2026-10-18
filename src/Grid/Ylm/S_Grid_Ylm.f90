! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Grid_Ylm) S_Grid_Ylm

  implicit none

  !=============================================================================
  ! local data
  !=============================================================================

  !> Gauss-Legendre tables for the m = 0 fast path of SpatialProduct (Setup):
  !> quadY(l, q) = Y_l0(x_q) and quadWY(q, l) = 2 pi w_q Y_l0(x_q) for
  !> 0 <= l <= quadLcap on the quadNq Gauss-Legendre nodes x_q in cos(theta).
  integer(I32) :: quadLcap = -1
  integer(I32) :: quadNq = 0
  real(R64), allocatable :: quadY(:, :)
  real(R64), allocatable :: quadWY(:, :)

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

    call SetupQuadratureM0

  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Tabulate Y_l0 on Gauss-Legendre nodes for the m = 0 fast path of
  !> SpatialProduct. The tables cover l <= 2 lmax + 2 (orbitals, densities and
  !> potentials of the default interaction layout, multiplicative operators)
  !> and integrate products of total degree <= 4 lmax + 5 exactly.
  subroutine SetupQuadratureM0
    use M_Utils_Constants
    use M_Grid_Ylm
    use stdlib_quadrature, only: gauss_legendre

    real(R64), allocatable :: x(:), w(:), p(:, :)
    integer(I32) :: l

    quadLcap = 2 * Grid_Ylm_lmax + 2
    quadNq = 2 * Grid_Ylm_lmax + 3

    allocate (x(quadNq), w(quadNq), p(0:quadLcap, quadNq))
    call gauss_legendre(x, w)

    ! Legendre polynomials by the Bonnet recurrence
    p(0, :) = 1.0_R64
    if (quadLcap >= 1) p(1, :) = x(:)
    do l = 1, quadLcap - 1
      p(l + 1, :) = ((2 * l + 1) * x(:) * p(l, :) - l * p(l - 1, :)) / (l + 1)
    end do

    if (allocated(quadY)) deallocate (quadY)
    if (allocated(quadWY)) deallocate (quadWY)
    allocate (quadY(0:quadLcap, quadNq), quadWY(quadNq, 0:quadLcap))

    do l = 0, quadLcap
      quadY(l, :) = sqrt((2 * l + 1) / (4.0_R64 * PI)) * p(l, :)
      quadWY(:, l) = 2.0_R64 * PI * w(:) * quadY(l, :)
    end do

    Grid_Ylm_nNodesM0 = quadNq
    Grid_Ylm_lmaxNodesM0 = quadLcap
    Grid_Ylm_nodeWeightsM0 = 2.0_R64 * PI * w(:)

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
    use M_Utils_SphericalHarmonics

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

    integer(I32) :: l, m
    integer(I32) :: l1, m1, l2, m2
    integer(I32) :: lMaxOut, lMaxConjg, lMax
    integer(I32) :: mMaxOut, mMaxConjg, mMax
    integer(I32) :: nRad, iCh, i1, i2, iOut
    real(R64)    :: gVal
    logical :: conjgQ, withWeightsQ
    logical, allocatable :: nonZeroQ(:)
    complex(R64), allocatable :: lmTmpConjg(:)

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

    nRad = Grid_Ylm_nRadial

    ! Axially symmetric operands: exact Gauss-Legendre quadrature in cos(theta)
    ! (O(lmax^2) per radial point instead of the O(lmax^3) Gaunt sum)
    if (mMaxOut .eq. 0 .and. mMaxConjg .eq. 0 .and. mMax .eq. 0 .and. &
        max(lMaxOut, lMaxConjg, lMax) <= quadLcap .and. lMaxOut + lMaxConjg + lMax < 2 * quadNq) then
      call SpatialProductM0(fOut, fConjg, f, conjgQ, lMaxOut, lMaxConjg, lMax)
      if (withWeightsQ) then
        do iCh = 1, lMaxOut + 1
          iOut = (iCh - 1) * nRad + 1
          fOut(iOut:iOut + nRad - 1) = fOut(iOut:iOut + nRad - 1) * Grid_Ylm_radialMetricWeights(:)
        end do
      end if
      return
    end if

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

    ! Apply metric weights if requested (complex on c-product back-ends)
    if (withWeightsQ) then
      do iCh = 1, Grid_Ylm_NChannels(lMaxOut, mMaxOut)
        iOut = (iCh - 1) * nRad + 1
        fOut(iOut:iOut + nRad - 1) = fOut(iOut:iOut + nRad - 1) * Grid_Ylm_radialMetricWeights(:)
      end do
    end if

    deallocate (lmTmpConjg, nonZeroQ)
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
    aRe = matmul(reshape(real(fConjg(1:nRad * (lMaxConjg + 1))), [nRad, lMaxConjg + 1]), quadY(0:lMaxConjg, :))
    aIm = matmul(reshape(aimag(fConjg(1:nRad * (lMaxConjg + 1))), [nRad, lMaxConjg + 1]), quadY(0:lMaxConjg, :))
    if (conjgQ) aIm = -aIm

    ! Second factor on the angular nodes
    bRe = matmul(reshape(real(f(1:nRad * (lMax + 1))), [nRad, lMax + 1]), quadY(0:lMax, :))
    bIm = matmul(reshape(aimag(f(1:nRad * (lMax + 1))), [nRad, lMax + 1]), quadY(0:lMax, :))

    ! Pointwise product, then projection onto Y_l0
    pRe = aRe * bRe - aIm * bIm
    pIm = aRe * bIm + aIm * bRe

    fOut(1:nRad * (lMaxOut + 1)) = reshape(cmplx(matmul(pRe, quadWY(:, 0:lMaxOut)), &
                                                 matmul(pIm, quadWY(:, 0:lMaxOut)), kind=R64), &
                                           [nRad * (lMaxOut + 1)])

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Grid_Ylm_ToNodesM0(fNodes, f, lmax)
    complex(R64), intent(out), contiguous :: fNodes(:, :)
    complex(R64), intent(in), contiguous  :: f(:)
    integer(I32), intent(in)              :: lmax

    integer(I32) :: nRad

    if (lmax > quadLcap) error stop "Grid_Ylm_ToNodesM0: lmax exceeds the quadrature tables"

    nRad = Grid_Ylm_nRadial

    fNodes(:, :) = cmplx(matmul(reshape(real(f(1:nRad * (lmax + 1))), [nRad, lmax + 1]), quadY(0:lmax, :)), &
                         matmul(reshape(aimag(f(1:nRad * (lmax + 1))), [nRad, lmax + 1]), quadY(0:lmax, :)), &
                         kind=R64)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Grid_Ylm_FromNodesM0(f, fNodes, lmax)
    complex(R64), intent(out), contiguous :: f(:)
    complex(R64), intent(in), contiguous  :: fNodes(:, :)
    integer(I32), intent(in)              :: lmax

    integer(I32) :: nRad

    if (lmax > quadLcap) error stop "Grid_Ylm_FromNodesM0: lmax exceeds the quadrature tables"

    nRad = Grid_Ylm_nRadial

    f = (0.0_R64, 0.0_R64)
    f(1:nRad * (lmax + 1)) = reshape(cmplx(matmul(real(fNodes), quadWY(:, 0:lmax)), &
                                           matmul(aimag(fNodes), quadWY(:, 0:lmax)), kind=R64), &
                                     [nRad * (lmax + 1)])

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure module function Grid_Ylm_NChannels(lmax, mmax) result(res)

    integer(I32), intent(in) :: lmax
    integer(I32), intent(in) :: mmax
    integer(I32) :: res

    res = ChannelOffset(lmax + 1, mmax)

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Number of channels preceding shell l in a layout with |m| <= min(l, mmax).
  !>
  !> Shells l' < l hold (2*min(l',mmax)+1) channels: l**2 while l <= mmax,
  !> afterwards each further shell adds (2*mmax+1).
  pure function ChannelOffset(l, mmax) result(res)

    integer(I32), intent(in) :: l
    integer(I32), intent(in) :: mmax
    integer(I32) :: res

    if (l <= mmax) then
      res = l * l
    else
      res = mmax * mmax + (l - mmax) * (2 * mmax + 1)
    end if

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> First flattened index of channel (l,m); channels are l-major with m
  !> running over -min(l,mmax)..min(l,mmax), each holding nRadial points.
  pure function BaseIndex(l, m, mmax) result(res)

    integer(I32), intent(in) :: l
    integer(I32), intent(in) :: m
    integer(I32), intent(in) :: mmax
    integer(I32) :: res

    res = 1 + (ChannelOffset(l, mmax) + m + min(l, mmax)) * Grid_Ylm_nRadial

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

    complex(R64), intent(out) :: fOut(:)   ! Full function on grid (target)
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
