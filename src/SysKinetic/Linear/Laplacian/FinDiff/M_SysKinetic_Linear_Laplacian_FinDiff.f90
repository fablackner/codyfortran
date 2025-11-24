! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Finite-difference Laplacian for linear grids.
!>
!> Provides a fabrication hook for an FD-based ∇² on uniform linear grids. The
!> concrete operator is selected and assigned by the linear fabric; this module
!> focuses on method-specific configuration (stencil, spacing, boundaries).
module M_SysKinetic_Linear_Laplacian_FinDiff
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Configure the FD Laplacian and bind the operator implementation.
    !>
    !> Expected settings (by convention): grid spacing(s), FD stencil order,
    !> and boundary conditions (Dirichlet/Neumann/periodic). The selected
    !> operator will be registered with the linear kinetic facade.
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
