! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Finite-Element Discrete Variable Representation (FEDVR) grid utilities.
!>
!> Provides types and builders for 1D FEDVR grids, including element-local nodes
!> and weights as well as global concatenated nodes/weights suitable for
!> integration and operator representation on nonuniform element partitions.
module M_Utils_Fedvr
  use M_Utils_Types

  implicit none

  !> Type definition for a single finite element in the FEDVR grid
  type :: T_Fedvr_Element
    !> Left boundary of the element
    real(R64) :: leftBound
    !> Right boundary of the element
    real(R64) :: rightBound
    !> Size of the element
    real(R64) :: size
    !> Global index of the first point for this element
    integer(I32) :: iStart
    !> Global index of the last point for this element
    integer(I32) :: iEnd
    !> element-local physical weights on [leftBound,rightBound]
    real(R64), allocatable :: weights(:)
  end type

  !> Type definition for an FEDVR grid configuration
  type :: T_Fedvr_Ctx
    !> Number of finite elements
    integer(I32) :: nElements
    !> uniform number of local points for all elements
    integer(I32) :: nLocals
    !> Flag to exclude global lower endpoint (xmin)
    logical :: xminExcludedQ
    !> Flag to exclude global upper endpoint (xmax)
    logical :: xmaxExcludedQ
    !> Total number of grid points
    integer(I32) :: nPoints
    !> Lower bound of the domain
    real(R64) :: xmin
    !> Upper bound of the domain
    real(R64) :: xmax

    ! Per-element data
    !> Array of element data
    type(T_Fedvr_Element), allocatable :: elements(:)

    ! Global grid points and weights
    !> Global grid points
    real(R64), allocatable :: points(:)
    !> Integration weights for each point
    real(R64), allocatable :: weights(:)
  end type

  ! Public interface
  public :: T_Fedvr_Ctx
  public :: T_Fedvr_Element
  public :: Fedvr_CreateCtx
  public :: Fedvr_CreateCtxFromBounds
  public :: Fedvr_DestroyCtx

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Creates and initializes an FEDVR grid of nE equally sized elements
!> @param ctx The FEDVR grid structure to initialize
!> @param xmin Lower bound of the domain
!> @param xmax Upper bound of the domain
!> @param nE Number of finite elements
!> @param nLocals Number of local points within each element
!> @param xminExcludedQ Logical flag to determine if global lower endpoint (xmin) is excluded
!> @param xmaxExcludedQ Logical flag to determine if global upper endpoint (xmax) is excluded
  subroutine Fedvr_CreateCtx(ctx, xmin, xmax, nE, nLocals, xminExcludedQ, xmaxExcludedQ)

    type(T_Fedvr_Ctx), intent(out) :: ctx
    real(R64), intent(in) :: xmin, xmax
    integer(I32), intent(in) :: nE, nLocals
    logical, intent(in) :: xminExcludedQ
    logical, intent(in) :: xmaxExcludedQ
    integer(I32) :: iE
    real(R64) :: elementSize
    real(R64), allocatable :: bounds(:)

    ! Uniform element boundaries xmin + iE * elementSize
    elementSize = (xmax - xmin) / nE

    allocate (bounds(0:nE))
    do iE = 0, nE
      bounds(iE) = xmin + iE * elementSize
    end do

    call Fedvr_CreateCtxFromBounds(ctx, bounds, nLocals, xminExcludedQ, xmaxExcludedQ)

    ! Keep the uniform bookkeeping bitwise as before: the requested domain
    ! limits (bounds(nE) may differ from xmax in the last bit) and the common
    ! element size instead of the individual bound differences
    ctx % xmin = xmin
    ctx % xmax = xmax
    ctx % elements(:) % size = elementSize

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Creates and initializes an FEDVR grid on explicit, possibly non-uniform
!> element boundaries (element iE spans [bounds(iE-1), bounds(iE)])
!> @param ctx The FEDVR grid structure to initialize
!> @param bounds Element boundaries bounds(0:nE), strictly increasing; nE = size(bounds) - 1
!> @param nLocals Number of local points within each element (>= 2)
!> @param xminExcludedQ Logical flag to determine if global lower endpoint (bounds(0)) is excluded
!> @param xmaxExcludedQ Logical flag to determine if global upper endpoint (bounds(nE)) is excluded
  subroutine Fedvr_CreateCtxFromBounds(ctx, bounds, nLocals, xminExcludedQ, xmaxExcludedQ)
    use stdlib_quadrature, only: gauss_legendre_lobatto

    type(T_Fedvr_Ctx), intent(out) :: ctx
    real(R64), intent(in) :: bounds(0:)
    integer(I32), intent(in) :: nLocals
    logical, intent(in) :: xminExcludedQ
    logical, intent(in) :: xmaxExcludedQ
    integer(I32) :: nE, iE, iLocal, iGrid
    real(R64), allocatable :: elementPoints(:)

    nE = size(bounds) - 1

    if (nE < 1) error stop "Fedvr_CreateCtxFromBounds: at least two element boundaries are required"
    if (nLocals < 2) error stop "Fedvr_CreateCtxFromBounds: nLocals must be >= 2"
    ! NaN-safe: a NaN bound fails the comparison
    if (.not. all(bounds(1:nE) > bounds(0:nE - 1))) error stop "Fedvr_CreateCtxFromBounds: bounds must be strictly increasing"

    ! Set grid parameters
    ctx % nElements = nE
    ctx % xmin = bounds(0)
    ctx % xmax = bounds(nE)
    ctx % nLocals = nLocals
    ctx % xminExcludedQ = xminExcludedQ
    ctx % xmaxExcludedQ = xmaxExcludedQ

    ! Allocate temporary arrays for element points and weights
    allocate (elementPoints(nLocals))

    ! Allocate and initialize elements
    allocate (ctx % elements(nE))
    do iE = 1, nE
      associate (element => ctx % elements(iE))
        element % leftBound = bounds(iE - 1)
        element % rightBound = bounds(iE)
        element % size = bounds(iE) - bounds(iE - 1)
      end associate
    end do

    ! Standard FEDVR has shared points at element boundaries
    ! Total is nE*nLocals - (nE-1) shared interior points
    ctx % nPoints = nE * (nLocals - 1) + 1
    if (ctx % xminExcludedQ) ctx % nPoints = ctx % nPoints - 1
    if (ctx % xmaxExcludedQ) ctx % nPoints = ctx % nPoints - 1

    allocate (ctx % points(ctx % nPoints))
    allocate (ctx % weights(ctx % nPoints))
    ctx % weights = 0.0_R64

    ! Compute ctx % points and ctx % weights
    iGrid = 0
    do iE = 1, nE
      associate (element => ctx % elements(iE))

        allocate (element % weights(nLocals))
        ! Compute quadrature points and weights for this element
        call gauss_legendre_lobatto(elementPoints, element % weights, [element % leftBound, element % rightBound])

        element % iStart = 0

        do iLocal = 1, ctx % nLocals
          if (iE .eq. 1 .and. iLocal .eq. 1 .and. ctx % xminExcludedQ) cycle
          if (iE .eq. nE .and. iLocal .eq. nLocals .and. ctx % xmaxExcludedQ) cycle

          if (iE .eq. 1 .or. iLocal .ne. 1) iGrid = iGrid + 1

          if (element % iStart .eq. 0) element % iStart = iGrid
          ctx % points(iGrid) = elementPoints(iLocal)
          ctx % weights(iGrid) = ctx % weights(iGrid) + element % weights(iLocal)
        end do

        ! Set end index for this element
        element % iEnd = iGrid
      end associate
    end do

    ! Deallocate temporary arrays
    deallocate (elementPoints)
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Destroys and deallocates an FEDVR grid context
!> @param ctx The FEDVR grid structure to destroy
  subroutine Fedvr_DestroyCtx(ctx)
    type(T_Fedvr_Ctx), intent(inout) :: ctx

    deallocate (ctx % elements, ctx % points, ctx % weights)

    ! Reset scalar members
    ctx % nElements = 0
    ctx % nLocals = 0
    ctx % nPoints = 0
    ctx % xmin = 0.0_R64
    ctx % xmax = 0.0_R64
    ctx % xminExcludedQ = .false.
    ctx % xmaxExcludedQ = .false.
  end subroutine
end module
