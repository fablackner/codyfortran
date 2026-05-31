! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Hydrogen-like radial orbital initializer for spherical grids.
!>
!> @details
!> Implements orbital initialization using the exact radial wavefunctions of
!> hydrogen-like (one-electron) atoms:
!>
!>   R_nl(r) ∝ ρ^l · L_{n-l-1}^{2l+1}(ρ) · exp(-ρ/2)
!>
!> where:
!> - ρ = 2Zr/(na₀)   is the dimensionless radial coordinate
!> - Z              is the effective nuclear charge
!> - a₀ = 1         is the Bohr radius (atomic units)
!> - L_k^α(ρ)       is the generalized Laguerre polynomial
!> - n, l, m        are the principal, angular, and magnetic quantum numbers
!>
!> The quantum numbers for each orbital are specified as arrays in JSON config.
!> The function returns zero for grid points where (l, m) doesn't match the
!> orbital's target quantum numbers.
!>
!> Physical Constraints
!> --------------------
!> - n ≥ 1          (principal quantum number)
!> - 0 ≤ l ≤ n-1    (angular momentum)
!> - -l ≤ m ≤ l     (magnetic quantum number)
!>
!> JSON Configuration
!> ------------------
!>   "orbsInit": {
!>     "ylm": {
!>       "hydrogenLike": {
!>         "charge": 1.0,     // effective Z (default: 1.0)
!>         "n": [1, 2, 2],    // principal QN per orbital
!>         "l": [0, 0, 1],    // angular QN per orbital
!>         "m": [0, 0, 0]     // magnetic QN per orbital
!>       }
!>     }
!>   }
!>
!> @note For multi-electron atoms, Z can be adjusted for screening effects.
module M_OrbsInit_Ylm_HydrogenLike
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Fabricate a hydrogen-like (Coulombic) Ylm initializer.
    module subroutine OrbsInit_Ylm_HydrogenLike_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Effective nuclear charge Z (in atomic units, default 1.0 for hydrogen).
  real(R64) :: OrbsInit_Ylm_HydrogenLike_charge

  !> Principal quantum number n for each orbital (array, size = nOrbitals).
  integer(I32), allocatable :: OrbsInit_Ylm_HydrogenLike_n(:)

  !> Orbital angular momentum l for each orbital (must satisfy 0 ≤ l ≤ n-1).
  integer(I32), allocatable :: OrbsInit_Ylm_HydrogenLike_l(:)

  !> Magnetic quantum number m for each orbital (must satisfy -l ≤ m ≤ l).
  integer(I32), allocatable :: OrbsInit_Ylm_HydrogenLike_m(:)

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
