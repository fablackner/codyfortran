! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Implementation of the cosine-profile absorber for linear grids.
!>
!> @details
!> Contains the setup and application routines for the cosine mask absorber.
!> The mask is precomputed once during setup and applied via element-wise
!> multiplication during propagation.
submodule(M_Absorber_Linear_Cosinus) S_Absorber_Linear_Cosinus

  implicit none

  !-----------------------------------------------------------------------------
  ! Module-level configuration and state
  !-----------------------------------------------------------------------------

  !> Coordinate threshold where absorption begins (|x| ≥ onset triggers damping)
  real(R64) :: onset

  !> Exponent denominator for the cosine taper: M = cos^(1/order)(...)
  !> Larger values produce sharper transitions.
  integer(I32) :: order

  !> Precomputed mask function M(x) evaluated at each grid point.
  !> Shape: (Grid_nPoints). Values in [0, 1].
  real(R64), allocatable :: maskFunction(:)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Fabrication routine: read configuration and bind procedure pointers.
  !>
  !> @details
  !> Reads `onset` and `order` from JSON configuration and wires the
  !> `Absorber_Setup` and `Absorber_ApplyAbsorber` pointers to this module's
  !> implementations.
  module subroutine Absorber_Linear_Cosinus_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Absorber

    call Say_Fabricate("absorber.linear.cosinus")

    !------------------------------------
    ! read configuration parameters
    !------------------------------------

    onset = Json_Get("absorber.linear.cosinus.onset", 100.0_R64)
    order = Json_Get("absorber.linear.cosinus.order", 6)

    !------------------------------------
    ! bind procedure pointers
    !------------------------------------

    Absorber_Setup => Setup
    Absorber_ApplyAbsorber => ApplyAbsorber

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Precompute the mask function on the spatial grid.
  !>
  !> @details
  !> Allocates `maskFunction` and evaluates:
  !>
  !>   M(x_i) = 1                              if |x_i| < onset
  !>          = cos^(1/order)(π/2 · ξ_i)      otherwise
  !>
  !> where ξ_i = (|x_i| - onset) / (x_max - onset) ∈ [0, 1].
  !>
  !> @pre  Grid must be set up (Grid_nPoints, Grid_Linear_xCoord available).
  !> @post maskFunction(:) is allocated and populated.
  subroutine Setup
    use M_Utils_Constants
    use M_Utils_Say
    use M_Grid
    use M_Grid_Linear

    integer(I32) :: iGrid
    real(R64)    :: xi      ! absolute coordinate |x|
    real(R64)    :: xiMin   ! = onset (start of absorbing layer)
    real(R64)    :: xiMax   ! = x_max (domain boundary)
    real(R64)    :: t       ! normalized position ξ ∈ [0, 1]

    call Say_Setup("absorber.linear.cosinus")

    allocate (maskFunction(Grid_nPoints))

    xiMin = onset
    xiMax = Grid_Linear_xCoord(Grid_nPoints)

    do iGrid = 1, Grid_nPoints

      xi = abs(Grid_Linear_xCoord(iGrid))

      t = (xi - xiMin) / (xiMax - xiMin)

      if (xi < xiMin) then
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
