! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> FEDVR-ECS derivative operators with complex weights.
!>
!> Builds first/second derivative matrices in the reference element and applies
!> them within an Exterior Complex Scaling grid context.
module M_Utils_DerivativeFedvrEcs
  use M_Utils_Types
  implicit none

  type :: T_DerivativeFedvrEcs_Ctx
    ! Reference operators on [-0.5, 0.5], size (nLocals, nLocals)
    real(R64), allocatable :: firstOrderMatrix(:, :)
    real(R64), allocatable :: secondOrderMatrix(:, :)
  end type

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine DerivativeFedvrEcs_CreateCtx(ctx, FedvrEcsCtx)
    use M_Utils_FedvrEcs, only: T_FedvrEcs_Ctx
    use stdlib_quadrature, only: gauss_legendre_lobatto

    type(T_DerivativeFedvrEcs_Ctx), intent(out) :: ctx
    type(T_FedvrEcs_Ctx), intent(in) :: FedvrEcsCtx

    integer(I32) :: i, j, k, nLoc
    real(R64), allocatable :: xRef(:), wRef(:)

    ! Allocate operator storage based on uniform nLocals
    nLoc = FedvrEcsCtx % nLocals
    allocate (ctx % firstOrderMatrix(nLoc, nLoc))
    allocate (ctx % secondOrderMatrix(nLoc, nLoc))
    ctx % firstOrderMatrix = 0.0_R64
    ctx % secondOrderMatrix = 0.0_R64

    ! Build reference nodes/weights on [-0.5, 0.5]
    allocate (xRef(nLoc), wRef(nLoc))
    call gauss_legendre_lobatto(xRef, wRef, [-0.5_R64, 0.5_R64])

    ! DRef(i_node, j_basis) = dL_j/dxi evaluated at node i on reference element
    do i = 1, nLoc
      do j = 1, nLoc
        ctx % firstOrderMatrix(i, j) = LagrangeDerivative(xRef, j, i, nLoc)
      end do
    end do

    ! KRef(i_basis, j_basis) = - sum_k wRef(k) * DRef(k, i) * DRef(k, j)
    do i = 1, nLoc
      do j = 1, nLoc
        do k = 1, nLoc
          ctx % secondOrderMatrix(i, j) = ctx % secondOrderMatrix(i, j) - &
                                          wRef(k) * ctx % firstOrderMatrix(k, i) * ctx % firstOrderMatrix(k, j)
        end do
      end do
    end do

    deallocate (xRef, wRef)
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine DerivativeFedvrEcs_DestroyCtx(ctx)
    type(T_DerivativeFedvrEcs_Ctx), intent(inout) :: ctx
    if (allocated(ctx % firstOrderMatrix)) deallocate (ctx % firstOrderMatrix)
    if (allocated(ctx % secondOrderMatrix)) deallocate (ctx % secondOrderMatrix)
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine DerivativeFedvrEcs_Do1stDerivative(dfDx, f, ctx, FedvrEcsCtx)
    use M_Utils_FedvrEcs, only: T_FedvrEcs_Ctx

    complex(R64), intent(out) :: dfDx(:)
    complex(R64), intent(in)  :: f(:)
    type(T_DerivativeFedvrEcs_Ctx), intent(in) :: ctx
    type(T_FedvrEcs_Ctx), intent(in) :: FedvrEcsCtx

    call ApplyDerivativeOperator(dfDx, f, ctx % firstOrderMatrix, FedvrEcsCtx)
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine DerivativeFedvrEcs_Do2ndDerivative(d2fDx2, f, ctx, FedvrEcsCtx)
    use M_Utils_FedvrEcs, only: T_FedvrEcs_Ctx

    complex(R64), intent(out) :: d2fDx2(:)
    complex(R64), intent(in) :: f(:)
    type(T_DerivativeFedvrEcs_Ctx), intent(in) :: ctx
    type(T_FedvrEcs_Ctx), intent(in) :: FedvrEcsCtx

    call ApplyDerivativeOperator(d2fDx2, f, ctx % secondOrderMatrix, FedvrEcsCtx)
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ApplyDerivativeOperator(res, f, operatorMatrix, FedvrEcsCtx)
    use M_Utils_FedvrEcs, only: T_FedvrEcs_Ctx

    complex(R64), intent(out) :: res(:)
    complex(R64), intent(in) :: f(:)
    real(R64), intent(in) :: operatorMatrix(:, :)
    type(T_FedvrEcs_Ctx), intent(in) :: FedvrEcsCtx

    integer(I32) :: iE, iStart, iEnd, nLoc, nE
    integer(I32) :: iLocalStart, iLocalEnd
    complex(R64) :: dx

    nLoc = FedvrEcsCtx % nLocals
    nE = FedvrEcsCtx % nElements
    res = (0.0_R64, 0.0_R64)

    do iE = 1, nE
      associate (element => FedvrEcsCtx % elements(iE))
        iStart = element % iStart
        iEnd = element % iEnd
        dx = element % size

        ! Set local index range based on element position
        iLocalStart = 1
        iLocalEnd = nLoc
        if (FedvrEcsCtx % xminExcludedQ .and. iE .eq. 1) iLocalStart = 2
        if (FedvrEcsCtx % xmaxExcludedQ .and. iE .eq. nE) iLocalEnd = nLoc - 1

        ! Apply operator matrix to f for this element
        res(iStart:iEnd) = res(iStart:iEnd) + &
                           matmul(operatorMatrix(iLocalStart:iLocalEnd, iLocalStart:iLocalEnd), &
                                  f(iStart:iEnd)) / dx
      end associate
    end do

    res(:) = res(:) / FedvrEcsCtx % weights(:)
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Function to compute the derivative of the Lagrange polynomial at a point
  !> @param x Array of grid points
  !> @param iLocal Index of the basis function (Lagrange polynomial)
  !> @param iLocalNode Index of the evaluation point
  !> @param nLocals Total number of local points
  pure function LagrangeDerivative(x, iLocal, iLocalNode, nLocals) result(res)
    real(R64), intent(in) :: x(:)
    integer(I32), intent(in) :: iLocal, iLocalNode, nLocals
    integer(I32) :: iNode
    real(R64) :: res

    if (iLocal .eq. iLocalNode) then
      res = 0.0_R64
      do iNode = 1, nLocals
        if (iNode .ne. iLocal) then
          res = res + 1.0_R64 / (x(iLocalNode) - x(iNode))
        end if
      end do
    else
      res = 1.0_R64 / (x(iLocal) - x(iLocalNode))
      do iNode = 1, nLocals
        if (iNode .ne. iLocal .and. iNode .ne. iLocalNode) then
          res = res * (x(iLocalNode) - x(iNode)) / (x(iLocal) - x(iNode))
        end if
      end do
    end if
  end function

end module
