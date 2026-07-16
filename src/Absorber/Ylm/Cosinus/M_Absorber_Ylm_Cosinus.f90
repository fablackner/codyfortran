! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Cosine-profile radial absorber for a Ylm grid.
!>
!> @details
!> Implements a smooth, reflection-minimizing mask function based on a cosine
!> profile in the radial coordinate of a Ylm grid. The mask is unity in the
!> interior region and decays smoothly to zero inside a boundary layer of
!> configurable thickness. Because the mask depends only on r, it acts
!> identically on every spherical-harmonic (l,m) channel.
!>
!> @par Mathematical Form
!> The mask function is defined as:
!>
!>   M(r) = cos^(1/n)(π/2 · ξ(r))   for  r ≥ onset
!>        = 1                        for  r < onset
!>
!> where:
!> - `onset` is the radius where absorption begins
!> - `n` is the "order" parameter controlling taper steepness (larger n = sharper)
!> - ξ(r) = (r - onset) / (r_max - onset) is the normalized position within
!>   the absorbing layer, mapping [onset, r_max] → [0, 1]
!>
!> At the outer boundary (ξ=1), M → 0. In the interior (r < onset), M = 1.
!>
!> @note The exponent 1/n (not n) is used, so **larger order values produce
!>       sharper transitions** (less gradual). Order 1 gives pure cos(πξ/2).
!>
!> @par JSON Configuration
!> @code{.json}
!> {
!>   "absorber": {
!>     "ylm": {
!>       "cosinus": {
!>         "onset": 16.0,
!>         "order": 6
!>       }
!>     }
!>   }
!> }
!> @endcode
!>
!> | Parameter | Type    | Default    | Description                                  |
!> |-----------|---------|------------|----------------------------------------------|
!> | onset     | real    | 0.8 · rmax | Radius where absorption begins (r ≥ onset)   |
!> | order     | integer | 6          | Exponent denominator (larger = sharper taper)|
!>
!> @see M_Absorber, M_Absorber_Ylm
module M_Absorber_Ylm_Cosinus
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> @brief Install the cosine-profile Ylm absorber.
    !>
    !> @details
    !> Binds `M_Absorber`'s procedure pointers to this implementation and reads
    !> configuration parameters (onset, order) from JSON.
    !>
    !> @pre  `Grid_Fabricate` must have run (default onset uses `Grid_Ylm_rmax`).
    !> @post `Absorber_Setup` → `Setup`, `Absorber_ApplyAbsorber` → `ApplyAbsorber`
    module subroutine Absorber_Ylm_Cosinus_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
