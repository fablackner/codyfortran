! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Prolate-spheroidal grid with FEDVR discretization of the xi coordinate.
!>
!> Uses Finite-Element Discrete Variable Representation (FEDVR) along xi on
!> [1, ximax]. The endpoint xi = 1 (the internuclear axis segment) is kept as
!> a grid point because sigma orbitals are nonzero there; ximax is excluded
!> (Dirichlet box boundary). The weak-form boundary term at xi = 1 vanishes
!> analytically since (xi^2 - 1) = 0 there.
module M_Grid_Prolate_Fedvr
  use M_Utils_Types
  use M_Utils_Fedvr, only: T_Fedvr_Ctx
  use M_Utils_DerivativeFedvr, only: T_DerivativeFedvr_Ctx

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Instantiate and configure the prolate xi FEDVR back-end.
    module subroutine Grid_Prolate_Fedvr_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Module-level FEDVR grid context (xi nodes, weights, element meta).
  type(T_Fedvr_Ctx) :: Grid_Prolate_Fedvr_fedvrCtx

  !> Module-level FEDVR derivative operator context (reference matrices).
  type(T_DerivativeFedvr_Ctx) :: Grid_Prolate_Fedvr_derivativeCtx

  !> Number of finite elements in the xi partition.
  integer(I32) :: Grid_Prolate_Fedvr_nElements

  !> Number of Gauss-Lobatto nodes per element (local basis size).
  integer(I32) :: Grid_Prolate_Fedvr_nLocals

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
