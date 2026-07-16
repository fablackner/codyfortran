! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Full-equation FEDVR-ECS solver for the radial Poisson equation.
!>
!> @details Solves the radial Poisson equation as a boundary-value problem
!> along the Exterior Complex Scaling contour z(r):
!>    (d²/dz² - l(l+1)/z²) u(z) = -4πz ρ(z)
!> with Robin boundary condition at the contour end for asymptotic decay.
!> The analytic continuation of the Coulomb kernel makes the entire solve
!> complex: element sizes, quadrature weights, and coordinates are the
!> complex contour quantities provided by Grid_Ylm_FedvrEcs.
!>
!> **Algorithm:**
!>   1. Assemble global complex FEDVR matrix (weak form) on the contour
!>   2. Apply boundary conditions (Dirichlet at z=0, Robin at contour end)
!>   3. Complex LU factorize and solve
!>
!> **Complexity:** O(N³) for LU factorization (done per call, not cached).
!> **Use case:** Coulomb mean fields on FEDVR-ECS grids (requires
!> `grid.ylm.fedvrEcs`); the sources are expected to carry the complex
!> contour metric weights (c-product convention).
module M_SysInteraction_Ylm_Coulomb_FullEqEcs
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Bind the full-equation ECS Coulomb solver to the core Ylm interaction hooks.
    module subroutine SysInteraction_Ylm_Coulomb_FullEqEcs_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
