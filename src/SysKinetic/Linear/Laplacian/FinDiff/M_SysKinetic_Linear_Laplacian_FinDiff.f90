! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file M_SysKinetic_Linear_Laplacian_FinDiff.f90
!> @brief Finite-difference Laplacian for linear grids.
!>
!> @details
!> Discretizes ∇² on uniform 1D grids using a central finite-difference stencil.
!> The default is a 5-point stencil (O(h⁴) accuracy):
!>
!>     ∇²ψ ≈ (−ψ_{i-2} + 16ψ_{i-1} − 30ψ_i + 16ψ_{i+1} − ψ_{i+2}) / (12 h²)
!>
!> Boundary Conditions
!> -------------------
!> Uses Dirichlet (zero) boundary conditions: ψ(x_min) = ψ(x_max) = 0. The
!> derivative routine handles boundary points by assuming values outside the
!> domain are zero.
!>
!> Requirements
!> ------------
!> Requires `grid.linear.const` to be configured (uniform grid spacing).
!>
!> @see M_Utils_DerivativeFinDiff for the underlying stencil implementation
module M_SysKinetic_Linear_Laplacian_FinDiff
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> @brief Configure the FD Laplacian and bind the operator implementation.
    !>
    !> Validates that `grid.linear.const` is available, then assigns
    !> `SysKinetic_MultiplyWithKineticOp` and `SysKinetic_Setup`.
    module subroutine SysKinetic_Linear_Laplacian_FinDiff_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
