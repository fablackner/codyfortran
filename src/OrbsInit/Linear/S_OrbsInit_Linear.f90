! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Implementation submodule for M_OrbsInit_Linear.
!>
!> @details
!> Provides the fabrication logic and the generic InitializeOrb routine that
!> samples the backend-specific InitFunction over the 1D grid and normalizes.
submodule(M_OrbsInit_Linear) S_OrbsInit_Linear

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Select and wire the linear initialization sub-backend from JSON config.
!>
!> @details
!> Currently supports: harmonic (quantum HO eigenstates).
!> Binds `OrbsInit_InitializeOrb` to the local InitializeOrb routine.
  module subroutine OrbsInit_Linear_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_OrbsInit
    use M_OrbsInit_Linear_Harmonic

    call Say_Fabricate("orbsInit.linear")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    OrbsInit_InitializeOrb => InitializeOrb

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("orbsInit.linear.harmonic")) then
      call OrbsInit_Linear_Harmonic_Fabricate

    else
      error stop "orbsInit.linear is missing one of: harmonic"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Initialize a single orbital on the 1D linear grid.
!>
!> @details
!> Samples the backend-specific `OrbsInit_Linear_InitFunction` at each grid
!> point, then normalizes using the grid's inner product (which accounts for
!> the spatial metric/quadrature weights).
!>
!> @param[out] orb  Complex orbital vector of length nGrid, filled and normalized.
!> @param[in]  ind  Orbital index (1-based; maps to quantum number n = ind - 1).
!> @param[in]  bt_  Optional body type for species-dependent initialization.
  subroutine InitializeOrb(orb, ind, bt_)
    use M_Grid
    use M_Grid_Linear

    complex(R64), intent(out), contiguous :: orb(:)
    integer(I32), intent(in)              :: ind
    integer(I32), intent(in), optional    :: bt_

    integer(I32) :: iGrid
    real(R64)    :: x, norm

    orb(:) = 0.0_R64

    ! Sample the initialization function at each grid point
    do iGrid = 1, Grid_nPoints
      x = Grid_Linear_xCoord(iGrid)
      orb(iGrid) = OrbsInit_Linear_InitFunction(x, ind, bt_)
    end do

    ! Normalize using grid-aware inner product (handles quadrature weights)
    norm = real(Grid_InnerProduct(orb, orb), kind=R64)
    orb(:) = orb(:) / sqrt(norm)

  end subroutine

end submodule
