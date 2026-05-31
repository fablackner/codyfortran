! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Coulomb interaction (1/r) for Ylm-basis atomic/molecular systems.
!>
!> @details Implements the electron-electron Coulomb repulsion using the
!> multipole expansion (Laplace expansion of 1/|r₁-r₂|).
!>
!> **Mathematical formulation:**
!> The radial Poisson equation for each (l,m) channel is:
!>    (d²/dr² - l(l+1)/r²) u(r) = -4πr ρₗₘ(r)
!> where V_lm(r) = u(r)/r satisfies the boundary conditions.
!>
!> **Available radial solvers:**
!>   - `StdImpl`: Direct O(N²) Green's function integration (reference)
!>   - `TwoScan`: O(N) forward/backward scan algorithm
!>   - `FullEq`: FEDVR-based differential equation with LU solve
!>   - `BlockEq`: Block-factorized FEDVR with precomputed Schur complement
!>
!> **JSON configuration:**
!> ```json
!> "sysInteraction": {
!>   "ylm": {
!>     "lmax": 4,
!>     "coulomb": {
!>       "strength": 1.0,
!>       "twoScan": {}
!>     }
!>   }
!> }
!> ```
module M_SysInteraction_Ylm_Coulomb
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Register a Coulomb solver variant (BlockEq, FullEq, TwoScan, ...)
    !> and propagate the strength and other settings to the back-end.
    module subroutine SysInteraction_Ylm_Coulomb_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Overall Coulomb coupling constant. Sign and units follow the global
  !> Hamiltonian conventions used in the simulation.
  real(R64) :: SysInteraction_Ylm_Coulomb_Strength

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
