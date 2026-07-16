! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Implementation of the cosine-profile radial absorber for Ylm grids.
!>
!> @details
!> Contains the setup and application routines for the radial cosine mask
!> absorber. The mask is precomputed once during setup on the flattened
!> (r,l,m) grid and applied via element-wise multiplication during propagation.
submodule(M_Absorber_Ylm_Cosinus) S_Absorber_Ylm_Cosinus

  implicit none

  !-----------------------------------------------------------------------------
  ! Module-level configuration and state
  !-----------------------------------------------------------------------------

  !> Radius threshold where absorption begins (r ≥ onset triggers damping)
  real(R64) :: onset

  !> Exponent denominator for the cosine taper: M = cos^(1/order)(...)
  !> Larger values produce sharper transitions.
  integer(I32) :: order

  !> Precomputed mask function M(r) evaluated at each flattened (r,l,m) grid
  !> point. Shape: (Grid_nPoints). Values in [0, 1].
  real(R64), allocatable :: maskFunction(:)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Fabrication routine: read configuration and bind procedure pointers.
  !>
  !> @details
  !> Reads `onset` and `order` from JSON configuration and wires the
  !> `Absorber_Setup` and `Absorber_ApplyAbsorber` pointers to this module's
  !> implementations. The onset defaults to 0.8 · `Grid_Ylm_rmax`, so
  !> `Grid_Fabricate` must have run beforehand.
  module subroutine Absorber_Ylm_Cosinus_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Absorber
    use M_Grid_Ylm

    call Say_Fabricate("absorber.ylm.cosinus")

    !------------------------------------
    ! read configuration parameters
    !------------------------------------

    onset = Json_Get("absorber.ylm.cosinus.onset", 0.8_R64 * Grid_Ylm_rmax)
    order = Json_Get("absorber.ylm.cosinus.order", 6)

    !------------------------------------
    ! bind procedure pointers
    !------------------------------------

    Absorber_Setup => Setup
    Absorber_ApplyAbsorber => ApplyAbsorber

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Precompute the radial mask function on the flattened Ylm grid.
  !>
  !> @details
  !> Allocates `maskFunction` and evaluates:
  !>
  !>   M(r_i) = 1                              if r_i < onset
  !>          = cos^(1/order)(π/2 · ξ_i)      otherwise
  !>
  !> where ξ_i = (r_i - onset) / (r_max - onset) ∈ [0, 1] and r_max is the
  !> outermost radial grid point. The mask is stored per flattened (r,l,m)
  !> point via `Grid_Ylm_rCoord`, so every (l,m) channel sees the same
  !> radial profile.
  !>
  !> @pre  Grid must be set up (Grid_nPoints, Grid_Ylm_rCoord available).
  !> @post maskFunction(:) is allocated and populated.
  subroutine Setup
    use M_Utils_Constants
    use M_Utils_Say
    use M_Grid
    use M_Grid_Ylm

    integer(I32) :: iGrid
    real(R64)    :: ri      ! radial coordinate r
    real(R64)    :: riMin   ! = onset (start of absorbing layer)
    real(R64)    :: riMax   ! = r_max (outermost radial grid point)
    real(R64)    :: t       ! normalized position ξ ∈ [0, 1]

    call Say_Setup("absorber.ylm.cosinus")

    allocate (maskFunction(Grid_nPoints))

    riMin = onset
    riMax = Grid_Ylm_radialPoints(Grid_Ylm_nRadial)

    if (riMax <= riMin) error stop "absorber.ylm.cosinus: onset must lie inside the radial grid"

    do iGrid = 1, Grid_nPoints

      ri = Grid_Ylm_rCoord(iGrid)

      t = (ri - riMin) / (riMax - riMin)

      if (ri < riMin) then
        maskFunction(iGrid) = 1.0_R64
      else
        maskFunction(iGrid) = cos(0.5_R64 * PI * t)**(1.0_R64 / order)
      end if

    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Apply the precomputed mask to all orbitals in place.
  !>
  !> @details
  !> Performs element-wise multiplication: orbs(:, i) ← orbs(:, i) * maskFunction(:)
  !> for each orbital i = 1, ..., nOrbsInState.
  !>
  !> @param[in,out] orbs  Orbital coefficient matrix (nGrid, nOrbitals).
  !>
  !> @pre  Setup must have been called (maskFunction allocated).
  subroutine ApplyAbsorber(orbs)
    use M_Orbs

    complex(R64), intent(inout), contiguous :: orbs(:, :)

    integer(I32) :: iOrb

    do iOrb = 1, Orbs_nOrbsInState
      orbs(:, iOrb) = orbs(:, iOrb) * maskFunction(:)
    end do

  end subroutine

end submodule
