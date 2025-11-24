! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Ylm grid with FEDVR radial discretization.
!>
!> Uses Finite-Element Discrete Variable Representation (FEDVR) along the radial
!> coordinate to achieve high accuracy with localized basis support. Well-suited
!> for atomic/ionic systems requiring specialized treatment of r.
module M_Grid_Ylm_Fedvr
  use M_Utils_Types
  use M_Utils_Fedvr, only: T_Fedvr_Ctx
  use M_Utils_DerivativeFedvr, only: T_DerivativeFedvr_Ctx

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Instantiate and configure a Ylm radial FEDVR back-end.
    module subroutine Grid_Ylm_Fedvr_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Module-level FEDVR grid context (radial nodes, weights, element meta).
  type(T_Fedvr_Ctx) :: Grid_Ylm_Fedvr_fedvrCtx

  !> Module-level FEDVR derivative operator context (1st/2nd derivatives, ...).
  type(T_DerivativeFedvr_Ctx) :: Grid_Ylm_Fedvr_derivativeCtx

  !> Number of finite elements in the radial grid partition.
  integer(I32) :: Grid_Ylm_Fedvr_nElements

  !> Number of Gauss–Lobatto nodes per element (local basis size).
  integer(I32) :: Grid_Ylm_Fedvr_nLocals

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
