! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Ylm grid with FEDVR radial discretization along a complex-scaled contour (ECS).
!>
!> Uses FEDVR in radius combined with Exterior Complex Scaling (ECS) beyond a
!> chosen radius to absorb outgoing flux. Configurable rotation onset, ramp,
!> and angle. Provides complex-valued contour points and weights.
module M_Grid_Ylm_FedvrEcs
  use M_Utils_Types
  use M_Utils_FedvrEcs, only: T_FedvrEcs_Ctx
  use M_Utils_DerivativeFedvrEcs, only: T_DerivativeFedvrEcs_Ctx

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Instantiate and configure a Ylm radial FEDVR-ECS back-end.
    module subroutine Grid_Ylm_FedvrEcs_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Module-level FEDVR-ECS grid context (complex contour, nodes, weights).
  type(T_FedvrEcs_Ctx) :: Grid_Ylm_FedvrEcs_FedvrEcsCtx

  !> Module-level derivative operator context adapted to ECS.
  type(T_DerivativeFedvrEcs_Ctx) :: Grid_Ylm_FedvrEcs_derivativeCtx

  !> Number of finite elements in the radial grid.
  integer(I32) :: Grid_Ylm_FedvrEcs_nElements

  !> Number of Gauss–Lobatto nodes per element.
  integer(I32) :: Grid_Ylm_FedvrEcs_nLocals

  !> First element index where ECS rotation starts (nElements+1 disables ECS).
  integer(I32) :: Grid_Ylm_FedvrEcs_firstEcsElement

  !> Number of elements used to ramp the ECS angle from 0 to theta.
  integer(I32) :: Grid_Ylm_FedvrEcs_transitionElements

  !> ECS rotation angle theta in radians.
  real(R64) :: Grid_Ylm_FedvrEcs_theta

  !> Radius r0 at which the ECS contour departs from the real axis.
  real(R64) :: Grid_Ylm_FedvrEcs_ecsRadius

  !> Complex-valued radial grid points following the ECS contour.
  complex(R64), allocatable :: Grid_Ylm_FedvrEcs_contourPoints(:)

  !> Complex-valued diagonal mass weights for the radial grid.
  complex(R64), allocatable :: Grid_Ylm_FedvrEcs_massWeights(:)

  !> Complex-valued volume weights (w*r^2) used in radial integrals.
  complex(R64), allocatable :: Grid_Ylm_FedvrEcs_volumeWeights(:)

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
