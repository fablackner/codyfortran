! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file M_SysKinetic_Ylm_Laplacian_FedvrEcs.f90
!> @brief FEDVR-ECS radial Laplacian for Ylm channels.
!>
!> @details
!> Provides the radial Laplacian on a finite-element DVR grid whose radial
!> coordinate follows an Exterior Complex Scaling (ECS) contour. All radial
!> factors (the g(z) = z·f transformation, the centrifugal term, and the
!> derivative element sizes) are evaluated at the complex contour points, so
!> outgoing flux is absorbed beyond the ECS radius without an explicit mask.
!>
!> Requirements
!> ------------
!> Requires `grid.ylm.fedvrEcs` to be configured. The FEDVR-ECS context,
!> contour points, and derivative matrices are provided by the
!> Grid_Ylm_FedvrEcs module.
!>
!> @see M_Utils_DerivativeFedvrEcs, M_Grid_Ylm_FedvrEcs
module M_SysKinetic_Ylm_Laplacian_FedvrEcs
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> @brief Configure the FEDVR-ECS radial Laplacian and bind the operator.
    !>
    !> Validates that `grid.ylm.fedvrEcs` is available, then assigns
    !> `SysKinetic_Ylm_MultiplyWithRadialKineticOp`. Setup is handled by
    !> the Grid_Ylm_FedvrEcs module (no separate Setup needed here).
    module subroutine SysKinetic_Ylm_Laplacian_FedvrEcs_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
