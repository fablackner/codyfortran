! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Slater-type-orbital (LCAO) initializer on grid.fem.
!>
!> @details
!> The input is a list of terms given as parallel arrays (JSON supports 1D
!> arrays only). Orbital ind is the sum over the terms t with
!> orbital(t) = ind,
!>
!>   phi_ind(r) = sum_t coeff_t rho^(n_t - l_t - 1) S_(l_t m_t)(r - R_t) exp(-zeta_t rho),
!>   rho = |r - R_t|,   R_t = Grid_Fem_centers(:, center_t),
!>
!> with the real regular solid harmonics S_lm in Racah normalization
!> (M_Utils_SolidHarmonics: S_00 = 1, S_10 = z, S_11 = x, S_1-1 = y,
!> S_20 = (3z^2 - r^2)/2, S_2-2 = sqrt(3) xy, ...), so rho^(n-l-1) S_lm is
!> r^(n-1) times an angular factor: a polynomial in the coordinates for
!> n - l - 1 = 0 and never a division by rho (regular on the nucleus nodes).
!> The terms are not normalized individually (coeff weights the unnormalized
!> primitives); the family normalizes the whole orbital on the grid. The same
!> orbitals are used for every body type.
!>
!> JSON Configuration (water: O 1s, 2s, 2p_x, 2p_z, 2p_y + H 1s mixing in 5)
!> ------------------
!>   "orbsInit": { "fem": { "sto": {
!>       "orbital": [1, 2, 3, 4, 5, 5, 5],
!>       "center":  [1, 1, 1, 1, 1, 2, 3],
!>       "n":       [1, 2, 2, 2, 2, 1, 1],
!>       "l":       [0, 0, 1, 1, 1, 0, 0],
!>       "m":       [0, 0, 1, 0, -1, 0, 0],
!>       "zeta":    [7.66, 2.25, 2.23, 2.23, 2.23, 1.0, 1.0],
!>       "coeff":   [1, 1, 1, 1, 1, 0.3, 0.3]
!>   } } }
!>
!> | Key       | Type       | Default | Description                                         |
!> |-----------|------------|---------|-----------------------------------------------------|
!> | `orbital` | integer(:) | —       | Orbital index of each term, >= 1; 1..max covered    |
!> | `center`  | integer(:) | —       | grid.fem center of each term, 1..nCenters           |
!> | `n`       | integer(:) | —       | Principal quantum number, n > l                     |
!> | `l`       | integer(:) | —       | Angular momentum, l >= 0                            |
!> | `m`       | integer(:) | —       | Real solid harmonic index, \|m\| <= l (m < 0: sin)  |
!> | `zeta`    | real(:)    | —       | Exponent, > 0                                       |
!> | `coeff`   | real(:)    | 1 each  | Coefficient of the term in its orbital              |
!>
!> All arrays have one entry per term. The model defines max(orbital)
!> orbitals; requesting more (method.mb.orbBased.nOrbs) stops.
module M_OrbsInit_Fem_Sto
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Read and validate the term arrays and bind the InitFunction.
    module subroutine OrbsInit_Fem_Sto_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Orbital index of every term (>= 1; every index up to the maximum occurs).
  integer(I32), allocatable :: OrbsInit_Fem_Sto_orbital(:)
  !> grid.fem center of every term (1..Grid_Fem_nCenters).
  integer(I32), allocatable :: OrbsInit_Fem_Sto_center(:)
  !> Principal quantum number n of every term (n > l).
  integer(I32), allocatable :: OrbsInit_Fem_Sto_n(:)
  !> Angular momentum l of every term (l >= 0).
  integer(I32), allocatable :: OrbsInit_Fem_Sto_l(:)
  !> Real solid harmonic index m of every term (|m| <= l).
  integer(I32), allocatable :: OrbsInit_Fem_Sto_m(:)
  !> Slater exponent zeta of every term (> 0).
  real(R64), allocatable :: OrbsInit_Fem_Sto_zeta(:)
  !> Coefficient of every term in its orbital.
  real(R64), allocatable :: OrbsInit_Fem_Sto_coeff(:)

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
