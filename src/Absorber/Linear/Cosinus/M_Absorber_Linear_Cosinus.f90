! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Cosine-profile absorber for a 1D linear grid.
!>
!> @details
!> Implements a smooth, reflection-minimizing mask function based on a cosine
!> profile for both ends of a 1D linear grid. The mask is unity in the interior
!> region and decays smoothly to zero inside boundary layers of configurable
!> thickness, thereby absorbing wavefunctions and reducing spurious reflections.
!>
!> @par Mathematical Form
!> The mask function is defined as:
!>
!>   M(x) = cos^(1/n)(π/2 · ξ(x))   for  |x| ≥ onset
!>        = 1                        for  |x| < onset
!>
!> where:
!> - `onset` is the coordinate where absorption begins
!> - `n` is the "order" parameter controlling taper steepness (larger n = sharper)
!> - ξ(x) = (|x| - onset) / (x_max - onset) is the normalized position within
!>   the absorbing layer, mapping [onset, x_max] → [0, 1]
!>
!> At the boundaries (ξ=1), M → 0. In the interior (|x| < onset), M = 1.
!>
!> @note The exponent 1/n (not n) is used, so **larger order values produce
!>       sharper transitions** (less gradual). Order 1 gives pure cos(πξ/2).
!>
!> @par JSON Configuration
!> @code{.json}
!> {
!>   "absorber": {
!>     "linear": {
!>       "cosinus": {
!>         "onset": 80.0,
!>         "order": 6
!>       }
!>     }
!>   }
!> }
!> @endcode
!>
!> | Parameter | Type    | Default | Description                                   |
!> |-----------|---------|---------|-----------------------------------------------|
!> | onset     | real    | 100.0   | Coordinate where absorption begins (|x| ≥ onset) |
!> | order     | integer | 6       | Exponent denominator (larger = sharper taper) |
!>
!> @see M_Absorber, M_Absorber_Linear
module M_Absorber_Linear_Cosinus
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> @brief Install the cosine-profile linear absorber.
    !>
    !> @details
    !> Binds `M_Absorber`'s procedure pointers to this implementation and reads
    !> configuration parameters (onset, order) from JSON.
    !>
    !> @post `Absorber_Setup` → `Setup`, `Absorber_ApplyAbsorber` → `ApplyAbsorber`
    module subroutine Absorber_Linear_Cosinus_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
