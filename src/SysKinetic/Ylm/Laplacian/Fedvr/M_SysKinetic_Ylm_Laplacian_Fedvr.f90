! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> FEDVR-based radial Laplacian for ylm channels.
!>
!> Provides a fabrication hook for the finite-element discrete variable
!> representation (FEDVR) discretization of the radial Laplacian. Suitable for
!> nonuniform radial grids with high-order local accuracy and sparse operators.
module M_SysKinetic_Ylm_Laplacian_Fedvr
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Configure the FEDVR radial Laplacian and bind the operator.
    !>
    !> Expects element partitions, local polynomial orders, and quadrature
    !> definitions. Precomputes derivative matrices and radial weights.
    module subroutine SysKinetic_Ylm_Laplacian_Fedvr_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
