! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief LCAO orbital initializer for the prolate-spheroidal grid.
!>
!> @details
!> Builds gerade/ungerade linear combinations of hydrogen-like atomic orbitals
!> centered on the two nuclei (foci):
!>
!>   psi = R_nl(r1) Theta_lm(cos th1) + s * R_nl(r2) Theta_lm(cos th2),
!>   s   = +(-1)^l for gerade, -(-1)^l for ungerade
!>
!> with r1 = a(xi+eta), cos th1 = (1+xi eta)/(xi+eta) (focus at z = -a) and
!> r2 = a(xi-eta), cos th2 = (xi eta-1)/(xi-eta) (focus at z = +a). The result
!> is placed in the azimuthal channel m; subsequent orthonormalization is
!> handled by the Orbs module.
!>
!> JSON Configuration
!> ------------------
!>   "orbsInit": {
!>     "prolate": {
!>       "lcao": {
!>         "charge": 1.0,           // effective Z per center
!>         "n": [1, 1],             // principal QN per orbital
!>         "l": [0, 0],             // angular QN of the atomic orbital
!>         "m": [0, 0],             // azimuthal QN (channel) per orbital
!>         "geradeQ": [true, false] // inversion parity per orbital
!>       }
!>     }
!>   }
module M_OrbsInit_Prolate_Lcao
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Fabricate the LCAO initializer.
    module subroutine OrbsInit_Prolate_Lcao_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Effective nuclear charge Z of each center.
  real(R64) :: OrbsInit_Prolate_Lcao_charge

  !> Principal quantum number n for each orbital (array, size = nOrbitals).
  integer(I32), allocatable :: OrbsInit_Prolate_Lcao_n(:)

  !> Atomic angular momentum l for each orbital (0 <= l <= n-1).
  integer(I32), allocatable :: OrbsInit_Prolate_Lcao_l(:)

  !> Azimuthal quantum number m for each orbital (|m| <= min(l, grid mmax)).
  integer(I32), allocatable :: OrbsInit_Prolate_Lcao_m(:)

  !> Inversion parity per orbital (true: gerade, false: ungerade).
  logical, allocatable :: OrbsInit_Prolate_Lcao_geradeQ(:)

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
