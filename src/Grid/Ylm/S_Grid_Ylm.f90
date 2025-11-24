! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Grid_Ylm) S_Grid_Ylm

  implicit none

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

    ! Read radial grid parameters from JSON
    Grid_Ylm_rmax = Json_Get("grid.ylm.rmax", 20.0_R64)

    Grid_InnerProduct => InnerProduct
    Grid_Ylm_SpatialProduct => SpatialProduct
    Grid_Ylm_GetLmComponent => GetLmComponent
    Grid_Ylm_SetLmComponent => SetLmComponent
    Grid_Ylm_AddLmComponent => AddLmComponent
    Grid_Setup => Setup

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
      do iM = -iL, iL
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
  pure subroutine SpatialProduct(fOut, fConjg, f, lMaxOut, lMaxConjg, lMax, conjgQ_, withWeightsQ_)
    use M_Utils_SphericalHarmonics

    complex(R64), intent(out), contiguous :: fOut(:)
    complex(R64), intent(in), contiguous  :: fConjg(:)
    complex(R64), intent(in), contiguous  :: f(:)
    integer(I32), intent(in)              :: lMaxOut
    integer(I32), intent(in)              :: lMaxConjg
    integer(I32), intent(in)              :: lMax
    logical, intent(in), optional         :: conjgQ_
    logical, intent(in), optional         :: withWeightsQ_

    integer(I32) :: l, m
    integer(I32) :: phase
    integer(I32) :: l1, m1, l2, m2
    integer(I32) :: nRad
    real(R64)    :: gVal
    logical :: conjgQ, withWeightsQ

    complex(R64), allocatable :: lmTmp(:), lmTmpConjg(:), fOutLm(:)

    ! Get number of radial points from the module
    nRad = Grid_Ylm_nRadial

    if (.not. present(conjgQ_)) conjgQ = .false.
    if (present(conjgQ_)) conjgQ = conjgQ_

    if (.not. present(withWeightsQ_)) withWeightsQ = .false.
    if (present(withWeightsQ_)) withWeightsQ = withWeightsQ_

    ! Allocate arrays for temporary radial functions
    allocate (lmTmp(nRad), lmTmpConjg(nRad), fOutLm(nRad))

    ! Initialize output array
    fOut = (0.0_R64, 0.0_R64)

    ! Build up expansion using Gaunt coefficients
    ! Loop over all (l₁, m₁) pairs
    do l1 = 0, lMaxConjg
      do m1 = -l1, l1
        if (conjgQ) then
          call Grid_Ylm_GetLmComponent(lmTmpConjg, l1, -m1, fConjg)
          ! Apply phase for complex conjugation
          if (mod(abs(m1), 2) .eq. 0) then
            phase = 1
          else
            phase = -1
          end if
          lmTmpConjg = phase * conjg(lmTmpConjg)
        else
          call Grid_Ylm_GetLmComponent(lmTmpConjg, l1, m1, fConjg)
        end if

        ! Check if this component is negligibly small
        if (all(abs(lmTmpConjg) < 1.0e-14_R64)) cycle

        ! Loop over all (l₂, m₂) pairs
        do l2 = 0, lMax
          do m2 = -l2, l2
            call Grid_Ylm_GetLmComponent(lmTmp, l2, m2, f)

            ! Check if this component is negligibly small
            if (all(abs(lmTmp) < 1.0e-14_R64)) cycle

            ! Calculate m value for output
            m = m1 + m2

            ! Expand Y_{l1,m1}(Ω)·Y_{l2,m2}(Ω) into Y_{l,m}(Ω) using Gaunt coefficients
            do l = abs(l1 - l2), min(l1 + l2, lMaxOut)
              if (abs(m) > l) cycle

              ! Calculate Gaunt coefficient
              gVal = SphericalHarmonics_GauntCoefficient( &
                     l1, m1, &
                     l2, m2, &
                     l, m)

              ! Skip if Gaunt coefficient is zero
              if (abs(gVal) < 1.0e-14_R64) cycle

              ! Calculate contribution to this (l,m) component
              fOutLm(:) = lmTmpConjg(:) * lmTmp(:) * gVal

              ! Add contribution directly to output
              call Grid_Ylm_AddLmComponent(fOut, l, m, fOutLm)
            end do
          end do
        end do
      end do
    end do

    ! Apply weights if requested
    if (withWeightsQ) then
      do l = 0, lMaxOut
        do m = -l, l
          call Grid_Ylm_GetLmComponent(lmTmp, l, m, fOut)
          lmTmp(:) = lmTmp(:) * Grid_Ylm_radialWeights(:)
          call Grid_Ylm_SetLmComponent(fOut, l, m, lmTmp)
        end do
      end do
    end if

    deallocate (lmTmp, lmTmpConjg, fOutLm)
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure subroutine GetLmComponent(fLmOut, l, m, f)

    complex(R64), intent(out) :: fLmOut(:)  ! Output radial function for (l,m)
    integer(I32), intent(in) :: l         ! Angular momentum quantum number
    integer(I32), intent(in) :: m         ! Magnetic quantum number
    complex(R64), intent(in) :: f(:)      ! Full function on grid

    integer(I32) :: iRad, iGrid, baseIdx, iL

    ! Initialize output
    fLmOut = (0.0_R64, 0.0_R64)

    ! Calculate base index for this (l,m) component
    ! For each l' < l, there are (2*l'+1) values of m, each with nRad points
    ! Within the current l shell, we offset by (m+l) * nRad to get to our m value
    baseIdx = 1
    ! Add all lower l shells
    do iL = 0, l - 1
      baseIdx = baseIdx + (2 * iL + 1) * Grid_Ylm_nRadial
    end do
    ! Add offset for m within current l shell
    baseIdx = baseIdx + (m + l) * Grid_Ylm_nRadial

    ! Extract the radial function for this (l,m)
    do iRad = 1, Grid_Ylm_nRadial
      iGrid = baseIdx + (iRad - 1)
      fLmOut(iRad) = f(iGrid)
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure subroutine SetLmComponent(fOut, l, m, fLm)

    complex(R64), intent(out) :: fOut(:)   ! Full function on grid (target)
    integer(I32), intent(in) :: l         ! Angular momentum quantum number
    integer(I32), intent(in) :: m         ! Magnetic quantum number
    complex(R64), intent(in) :: fLm(:)   ! Input radial function for (l,m)

    integer(I32) :: iRad, iGrid, baseIdx, iL

    ! Calculate base index for this (l,m) component
    ! For each l' < l, there are (2*l'+1) values of m, each with nRad points
    ! Within the current l shell, we offset by (m+l) * nRad to get to our m value
    baseIdx = 1
    ! Add all lower l shells
    do iL = 0, l - 1
      baseIdx = baseIdx + (2 * iL + 1) * Grid_Ylm_nRadial
    end do
    ! Add offset for m within current l shell
    baseIdx = baseIdx + (m + l) * Grid_Ylm_nRadial

    ! Add the radial function for this (l,m) to existing values
    do iRad = 1, Grid_Ylm_nRadial
      iGrid = baseIdx + (iRad - 1)
      fOut(iGrid) = fLm(iRad)  ! Add instead of replace
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure subroutine AddLmComponent(fInOut, l, m, fLm)

    complex(R64), intent(inout) :: fInOut(:)   ! Full function on grid (target)
    integer(I32), intent(in) :: l         ! Angular momentum quantum number
    integer(I32), intent(in) :: m         ! Magnetic quantum number
    complex(R64), intent(in) :: fLm(:)   ! Input radial function for (l,m)

    integer(I32) :: iRad, iGrid, baseIdx, iL

    ! Calculate base index for this (l,m) component
    ! For each l' < l, there are (2*l'+1) values of m, each with nRad points
    ! Within the current l shell, we offset by (m+l) * nRad to get to our m value
    baseIdx = 1
    ! Add all lower l shells
    do iL = 0, l - 1
      baseIdx = baseIdx + (2 * iL + 1) * Grid_Ylm_nRadial
    end do
    ! Add offset for m within current l shell
    baseIdx = baseIdx + (m + l) * Grid_Ylm_nRadial

    ! Add the radial function for this (l,m) to existing values
    do iRad = 1, Grid_Ylm_nRadial
      iGrid = baseIdx + (iRad - 1)
      fInOut(iGrid) = fInOut(iGrid) + fLm(iRad)  ! Add instead of replace
    end do

  end subroutine

end submodule
