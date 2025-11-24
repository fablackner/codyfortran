! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Submodule implementing Finite-Element Discrete Variable Representation (FEDVR)
!> in one dimension. FEDVR combines finite element methods with discrete variable
!> representation for highly accurate spatial discretization.
submodule(M_Grid_Linear_Fedvr) S_Grid_Linear_Fedvr
  use M_Utils_Fedvr, only: T_Fedvr_Ctx
  use M_Utils_DerivativeFedvr, only: T_DerivativeFedvr_Ctx

  implicit none

  type(T_Fedvr_Ctx) :: fedvrCtx
  type(T_DerivativeFedvr_Ctx) :: derivativeCtx

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !> Sets up the FEDVR grid based on parameters from the JSON configuration
  module subroutine Grid_Linear_Fedvr_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Grid

    call Say_Fabricate("grid.linear.fedvr1d")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    ! Get configuration parameters from JSON
    Grid_Linear_Fedvr_nElements = Json_Get("grid.spatial.fedvr1d.nElements", 10)
    Grid_Linear_Fedvr_order = Json_Get("grid.spatial.fedvr1d.order", 8)

    ! Set total grid size based on elements and order
    ! Each element has (order+1) points, but adjacent elements share points
    Grid_nPoints = Grid_Linear_Fedvr_nElements * Grid_Linear_Fedvr_order + 1

    ! Set procedure pointers
    Grid_Setup => Setup

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !> Initializes the FEDVR grid by setting up points, weights, and operators
  subroutine Setup
    use M_Utils_Say
    use M_Utils_Fedvr
    use M_Utils_DerivativeFedvr
    use M_Grid
    use M_Grid_Linear

    integer(I32) :: i

    call Say_Setup("grid.linear.fedvr")

    ! Create the FEDVR grid using the utility module
    call Fedvr_CreateCtx(fedvrCtx, Grid_Linear_xmin, Grid_Linear_xmax, &
                         Grid_Linear_Fedvr_nElements, Grid_Linear_Fedvr_order, .false., .false.)

    call DerivativeFedvr_CreateCtx(derivativeCtx, fedvrCtx)

    ! Allocate grid arrays
    allocate (Grid_Linear_xCoord(Grid_nPoints))  ! Store 1D coordinates
    allocate (Grid_Linear_weights(Grid_nPoints))

    ! Set grid coordinates and weights
    do i = 1, Grid_nPoints
      Grid_Linear_xCoord(i) = fedvrCtx % points(i)
      Grid_Linear_weights(i) = fedvrCtx % weights(i)
    end do

  end subroutine

end submodule
