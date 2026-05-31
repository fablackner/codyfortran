! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Full-equation FEDVR solver for radial Poisson equation.
!>
!> @details Solves the radial Poisson equation as a boundary-value problem
!> using the Finite Element Discrete Variable Representation (FEDVR):
!>    (d²/dr² - l(l+1)/r²) u(r) = -4πr ρ(r)
!> with Robin boundary condition at r_max for asymptotic decay.
!>
!> **Algorithm:**
!>   1. Assemble global FEDVR matrix (weak form)
!>   2. Apply boundary conditions (Dirichlet at r=0, Robin at r_max)
!>   3. LU factorize and solve
!>
!> **Complexity:** O(N³) for LU factorization (done per call, not cached).
!> **Use case:** High accuracy with FEDVR grids, moderate system sizes.
module M_SysInteraction_Ylm_Coulomb_FullEq
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Bind the full-equation Coulomb solver to the core Ylm interaction hooks.
    module subroutine SysInteraction_Ylm_Coulomb_FullEq_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
