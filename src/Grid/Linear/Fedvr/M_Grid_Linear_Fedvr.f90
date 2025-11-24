! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> 1D linear grid using Finite-Element Discrete Variable Representation (FEDVR).
!>
!> FEDVR combines finite elements with a discrete variable representation to
!> achieve high accuracy and spectral-like properties on each element. This
!> module stores the configuration for a 1D FEDVR grid and wires the related
!> operators via its Fabricate routine.
module M_Grid_Linear_Fedvr
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Instantiate and configure a 1D FEDVR grid.
    !>
    !> Allocates element-local structures, initializes quadrature nodes and
    !> weights, and associates FEDVR-specific operators (mass/derivative, …).
    module subroutine Grid_Linear_Fedvr_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Number of finite elements in the 1D FEDVR grid (global partition).
  !> Increasing elements refines the mesh and local support.
  integer(I32) :: Grid_Linear_Fedvr_nElements

  !> Polynomial order within each element (number of local basis points - 1).
  !> Higher order increases accuracy for smooth functions at higher cost.
  integer(I32) :: Grid_Linear_Fedvr_order

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
