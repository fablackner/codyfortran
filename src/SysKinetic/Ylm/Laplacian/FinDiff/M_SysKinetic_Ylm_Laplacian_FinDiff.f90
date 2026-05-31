! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file M_SysKinetic_Ylm_Laplacian_FinDiff.f90
!> @brief Finite-difference radial Laplacian for Ylm channels.
!>
!> @details
!> Discretizes the radial part of the Laplacian on a uniform radial grid using
!> central finite differences. Uses the g(r) = r·f(r) transformation to avoid
!> the 2/r singularity.
!>
!> Boundary Conditions
!> -------------------
!> Assumes f(0) = 0 (regularity at origin) and f(r_max) → 0 (Dirichlet).
!> The grid starts at r = dr (first nonzero point).
!>
!> Requirements
!> ------------
!> Requires `grid.ylm.const` to be configured (uniform radial spacing).
!>
!> @see M_Utils_DerivativeFinDiff
module M_SysKinetic_Ylm_Laplacian_FinDiff
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> @brief Configure the FinDiff radial Laplacian and bind the operator.
    !>
    !> Validates that `grid.ylm.const` is available, then assigns
    !> `SysKinetic_Ylm_MultiplyWithRadialKineticOp` and `SysKinetic_Setup`.
    module subroutine SysKinetic_Ylm_Laplacian_FinDiff_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
