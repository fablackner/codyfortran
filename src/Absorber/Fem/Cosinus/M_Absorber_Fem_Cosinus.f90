! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Cosine-profile boundary-layer absorber for the FEM grid.
!>
!> @details
!> The mask is a product of one-dimensional layer profiles along the three
!> axes,
!>
!>   M(x, y, z) = m_x(x) m_y(y) m_z(z),
!>   m_a(s) = cos^(1/order)(pi/2 (width - d_a)/width)   for d_a < width,
!>          = 1                                          otherwise,
!>
!> where d_a = min(s - lower_a, upper_a - s) is the distance of the node to the
!> nearer face of axis a. The mask is 1 in the inner box, decays
!> monotonically inside the layers of thickness `width` along the faces and
!> would vanish on the faces themselves (which are Dirichlet nodes, not
!> unknowns); edges and corners of the box get the product of the profiles.
!> It is built once in Setup (Grid_nPoints values) and ApplyAbsorber
!> multiplies every orbital column by it.
!>
!> @par JSON Configuration
!> @code{.json}
!> {
!>   "absorber": {
!>     "fem": {
!>       "cosinus": {
!>         "width": 4.0,
!>         "order": 8
!>       }
!>     }
!>   }
!> }
!> @endcode
!>
!> | Parameter | Type    | Default                    | Description                                     |
!> |-----------|---------|----------------------------|-------------------------------------------------|
!> | width     | real    | 0.2 · smallest box side    | Layer thickness, 0 < width <= smallest side / 2 |
!> | order     | integer | 6                          | Exponent denominator (larger = sharper), >= 1   |
!>
!> @see M_Absorber, M_Absorber_Fem
module M_Absorber_Fem_Cosinus
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> @brief Install the cosine-profile FEM absorber.
    !>
    !> @details
    !> Reads and validates width and order (the default width uses the box of
    !> grid.fem, so `Grid_Fabricate` must have run) and binds
    !> `Absorber_Setup` → `Setup`, `Absorber_ApplyAbsorber` → `ApplyAbsorber`.
    module subroutine Absorber_Fem_Cosinus_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Thickness of the absorbing layer along every face of the box.
  real(R64) :: Absorber_Fem_Cosinus_width = 0.0_R64

  !> Exponent denominator of the profile cos^(1/order) (larger = sharper).
  integer(I32) :: Absorber_Fem_Cosinus_order = 6

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
