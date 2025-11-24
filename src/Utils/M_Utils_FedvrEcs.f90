! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> FEDVR with Exterior Complex Scaling (ECS) grid utilities.
!>
!> Provides complex-valued element bounds, weights, and global integration
!> arrays to support exterior complex scaling regions for absorptive boundaries.
module M_Utils_FedvrEcs
  use M_Utils_Types

  implicit none

  !> Type definition for a single finite element in the FedvrEcs grid
  type :: T_FedvrEcs_Element
    !> Left boundary of the element (real)
    real(R64)    :: leftBoundReal
    !> Right boundary of the element (real)
    real(R64)    :: rightBoundReal
    !> Left boundary of the element (complex)
    complex(R64) :: leftBound
    !> Right boundary of the element (complex)
    complex(R64) :: rightBound
    !> Size of the element (complex)
    complex(R64) :: size
    !> Scaling factor for the element
    complex(R64) :: scaleFactor
    !> Global index of the first point for this element
    integer(I32) :: iStart
    !> Global index of the last point for this element
    integer(I32) :: iEnd
    !> element-local physical weights on [leftBound,rightBound]
    complex(R64), allocatable :: weights(:)
  end type

  !> Type definition for an FedvrEcs grid configuration
  type :: T_FedvrEcs_Ctx
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

    ! Configuration for exterior complex scaling
    !> First element for which complex scaling is applied
    integer(I32) :: firstEcsElement
    !> Number of elements over which to transition scaling
    integer(I32) :: transitionElements
    !> Angle for complex scaling
    real(R64)    :: ecsAngle
    !> Radius for complex scaling
    real(R64)    :: ecsRadius
    !> Flag to enable/disable complex scaling
    logical      :: ecsEnabledQ

    ! Per-element data
    !> Array of element data
    type(T_FedvrEcs_Element), allocatable :: elements(:)

    ! Global grid points and weights
    !> Global grid points (real)
    real(R64), allocatable :: pointsReal(:)
    !> Global grid points (complex)
    complex(R64), allocatable :: points(:)
    !> Integration weights for each point
    complex(R64), allocatable :: weights(:)
  end type

  ! Public interface
  public :: T_FedvrEcs_Ctx
  public :: T_FedvrEcs_Element
  public :: FedvrEcs_CreateCtx
  public :: FedvrEcs_DestroyCtx

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Creates and initializes an FedvrEcs grid
!> @param ctx The FedvrEcs grid structure to initialize
!> @param xmin Lower bound of the domain
!> @param xmax Upper bound of the domain
!> @param nE Number of finite elements
!> @param nLocals Number of local points within each element
!> @param xminExcludedQ Logical flag to determine if global lower endpoint (xmin) is excluded
!> @param xmaxExcludedQ Logical flag to determine if global upper endpoint (xmax) is excluded
!> @param ecsStartElement First element for complex scaling (optional)
!> @param ecsAngle Angle for complex scaling (optional)
!> @param ecsTransitionElements Number of elements for scaling transition (optional)
  subroutine FedvrEcs_CreateCtx(ctx, xmin, xmax, nE, nLocals, xminExcludedQ, xmaxExcludedQ, &
                                ecsStartElement, ecsAngle, ecsTransitionElements)
    use stdlib_quadrature, only: gauss_legendre_lobatto

    type(T_FedvrEcs_Ctx), intent(out) :: ctx
    real(R64), intent(in) :: xmin, xmax
    integer(I32), intent(in) :: nE, nLocals
    logical, intent(in) :: xminExcludedQ
    logical, intent(in) :: xmaxExcludedQ
    integer(I32), intent(in), optional :: ecsStartElement
    real(R64), intent(in), optional :: ecsAngle
    integer(I32), intent(in), optional :: ecsTransitionElements

    integer(I32) :: iE, iLocal, iGrid
    real(R64) :: elementSizeReal
    real(R64), allocatable :: elementPoints(:), elementWeights(:)
    real(R64) :: xi, localAngle
    integer(I32) :: firstEcsElement, transitionElements
    real(R64) :: targetAngle
    logical :: ecsEnabled
    complex(R64) :: contourLeft, weight, point
    complex(R64) :: scaleFactor, elementLength

    ! Set grid parameters
    ctx % nElements = nE
    ctx % xmin = xmin
    ctx % xmax = xmax
    ctx % nLocals = nLocals
    ctx % xminExcludedQ = xminExcludedQ
    ctx % xmaxExcludedQ = xmaxExcludedQ
    elementSizeReal = (xmax - xmin) / real(nE, R64)

    ! Determine complex scaling parameters
    firstEcsElement = nE + 1
    if (present(ecsStartElement)) firstEcsElement = ecsStartElement
    firstEcsElement = max(1, min(firstEcsElement, nE + 1))

    targetAngle = 0.0_R64
    if (present(ecsAngle)) targetAngle = ecsAngle

    transitionElements = 0
    if (present(ecsTransitionElements)) transitionElements = max(0, ecsTransitionElements)

    ecsEnabled = (abs(targetAngle) > 0.0_R64) .and. (firstEcsElement <= nE)
    ctx % ecsEnabledQ = ecsEnabled
    ctx % firstEcsElement = firstEcsElement
    ctx % transitionElements = transitionElements
    ctx % ecsAngle = targetAngle
    if (ecsEnabled) then
      ctx % ecsRadius = xmin + real(firstEcsElement - 1, R64) * elementSizeReal
    else
      ctx % ecsRadius = xmax
    end if

    ! Allocate temporary arrays for element points and weights
    allocate (elementPoints(nLocals), elementWeights(nLocals))
    ! Allocate and initialize elements
    allocate (ctx % elements(nE))
    contourLeft = cmplx(xmin, 0.0_R64, kind=R64)
    do iE = 1, nE
      associate (element => ctx % elements(iE))
        element % leftBoundReal = xmin + real(iE - 1, R64) * elementSizeReal
        element % rightBoundReal = xmin + real(iE, R64) * elementSizeReal

        if (.not. ecsEnabled .or. iE < firstEcsElement) then
          localAngle = 0.0_R64
        else if (transitionElements > 0 .and. iE < firstEcsElement + transitionElements) then
          localAngle = targetAngle * real(iE - firstEcsElement + 1, R64) / real(transitionElements + 1, R64)
        else
          localAngle = targetAngle
        end if

        scaleFactor = cmplx(cos(localAngle), sin(localAngle), kind=R64)
        elementLength = scaleFactor * elementSizeReal

        element % scaleFactor = scaleFactor
        element % leftBound = contourLeft
        element % rightBound = contourLeft + elementLength
        element % size = elementLength
        contourLeft = element % rightBound

        allocate (element % weights(nLocals))
        element % weights = (0.0_R64, 0.0_R64)
        element % iStart = 0
        element % iEnd = 0
      end associate
    end do

    ! Standard FedvrEsc has shared points at element boundaries
    ! Total is nE*nLocals - (nE-1) shared interior points
    ctx % nPoints = nE * (nLocals - 1) + 1
    if (ctx % xminExcludedQ) ctx % nPoints = ctx % nPoints - 1
    if (ctx % xmaxExcludedQ) ctx % nPoints = ctx % nPoints - 1

    allocate (ctx % points(ctx % nPoints))
    allocate (ctx % pointsReal(ctx % nPoints))
    allocate (ctx % weights(ctx % nPoints))
    ctx % points = (0.0_R64, 0.0_R64)
    ctx % pointsReal = 0.0_R64
    ctx % weights = (0.0_R64, 0.0_R64)

    ! Compute ctx % points and ctx % weights
    iGrid = 0
    do iE = 1, nE
      associate (element => ctx % elements(iE))

        ! Compute quadrature points and weights for this element
        call gauss_legendre_lobatto(elementPoints, elementWeights, [element % leftBoundReal, element % rightBoundReal])
        do iLocal = 1, ctx % nLocals
          if (iE .eq. 1 .and. iLocal .eq. 1 .and. ctx % xminExcludedQ) cycle
          if (iE .eq. nE .and. iLocal .eq. ctx % nLocals .and. ctx % xmaxExcludedQ) cycle

          if (iE .eq. 1 .or. iLocal .ne. 1) iGrid = iGrid + 1

          xi = (elementPoints(iLocal) - element % leftBoundReal) / (element % rightBoundReal - element % leftBoundReal)
          point = element % leftBound + xi * element % size
          weight = elementWeights(iLocal) * element % scaleFactor

          if (element % iStart .eq. 0) element % iStart = iGrid
          element % weights(iLocal) = weight

          ctx % points(iGrid) = point
          ctx % pointsReal(iGrid) = elementPoints(iLocal)
          ctx % weights(iGrid) = ctx % weights(iGrid) + weight
        end do

        ! Set end index for this element
        element % iEnd = iGrid
      end associate
    end do

    ! Deallocate temporary arrays
    deallocate (elementPoints, elementWeights)
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Destroys and deallocates an FedvrEcs grid context
!> @param ctx The FedvrEcs grid structure to destroy
  subroutine FedvrEcs_DestroyCtx(ctx)
    type(T_FedvrEcs_Ctx), intent(inout) :: ctx
    integer(I32) :: iE

    if (allocated(ctx % elements)) then
      do iE = 1, size(ctx % elements)
        if (allocated(ctx % elements(iE) % weights)) deallocate (ctx % elements(iE) % weights)
      end do
      deallocate (ctx % elements)
    end if

    if (allocated(ctx % points)) deallocate (ctx % points)
    if (allocated(ctx % pointsReal)) deallocate (ctx % pointsReal)
    if (allocated(ctx % weights)) deallocate (ctx % weights)

    ! Reset scalar members
    ctx % nElements = 0
    ctx % nLocals = 0
    ctx % nPoints = 0
    ctx % xmin = 0.0_R64
    ctx % xmax = 0.0_R64
    ctx % xminExcludedQ = .false.
    ctx % xmaxExcludedQ = .false.
    ctx % firstEcsElement = 0
    ctx % transitionElements = 0
    ctx % ecsAngle = 0.0_R64
    ctx % ecsRadius = 0.0_R64
    ctx % ecsEnabledQ = .false.
  end subroutine

end module
