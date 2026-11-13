! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Grid_Prolate_Fedvr) S_Grid_Prolate_Fedvr

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Grid_Prolate_Fedvr_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Grid
    use M_Grid_Prolate

    call Say_Fabricate("grid.prolate.fedvr")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    Grid_Prolate_Fedvr_nElements = Json_Get("grid.prolate.fedvr.nElements", 10)
    Grid_Prolate_Fedvr_nLocals = Json_Get("grid.prolate.fedvr.nLocals", 7)
    Grid_Prolate_Fedvr_grading = Json_Get("grid.prolate.fedvr.grading", 1.0_R64)
    if (.not. (Grid_Prolate_Fedvr_grading > 0.0_R64)) error stop "grid.prolate.fedvr.grading must be positive"

    ! Total xi points: nE*(nLocals-1) + 1 minus the excluded ximax endpoint
    Grid_Prolate_nXi = Grid_Prolate_Fedvr_nElements * (Grid_Prolate_Fedvr_nLocals - 1)

    ! Update total grid size to include the azimuthal channels
    Grid_nPoints = Grid_Prolate_nXi * Grid_Prolate_nEta * (2 * Grid_Prolate_mmax + 1)

    ! Set the setup procedure
    Grid_Prolate_Setup => Setup

    !------------------------------------
    ! branch
    !------------------------------------

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Setup
    use stdlib_quadrature, only: gauss_legendre_lobatto
    use M_Utils_Say
    use M_Utils_Fedvr
    use M_Utils_DerivativeFedvr
    use M_Grid
    use M_Grid_Prolate

    integer(I32) :: iE, kLocal, iLocal1, iLocal2, iLocalStart, iLocalEnd, iGlobal1, iGlobal2
    real(R64) :: tNode, mVal, mDer
    real(R64), allocatable :: elementPoints(:), elementWeights(:), gOdd(:, :), vOdd(:, :)

    call Say_Setup("grid.prolate.fedvr")

    allocate (Grid_Prolate_xiPoints(Grid_Prolate_nXi))
    allocate (Grid_Prolate_xiWeights(Grid_Prolate_nXi))

    ! xi = 1 is kept (sigma orbitals are nonzero on the internuclear axis),
    ! ximax is excluded (Dirichlet box boundary)
    call Fedvr_CreateCtxGraded(Grid_Prolate_Fedvr_fedvrCtx, &
                               1.0_R64, &
                               Grid_Prolate_ximax, &
                               Grid_Prolate_Fedvr_nElements, &
                               Grid_Prolate_Fedvr_nLocals, &
                               Grid_Prolate_Fedvr_grading, &
                               .false., &
                               .true.)

    call DerivativeFedvr_CreateCtx(Grid_Prolate_Fedvr_derivativeCtx, Grid_Prolate_Fedvr_fedvrCtx)

    Grid_Prolate_xiPoints(:) = Grid_Prolate_Fedvr_fedvrCtx % points(1:Grid_Prolate_nXi)
    Grid_Prolate_xiWeights(:) = Grid_Prolate_Fedvr_fedvrCtx % weights(1:Grid_Prolate_nXi)

    !------------------------------------
    ! xi matrices of the m channel, K_m = S_xi(:,:,p) + m^2 C_xi(:,:,p) (p = mod(|m|, 2)):
    !   S_xi(i,j,p) = sum_k w_k (xi_k^2 - 1) chi_i'(xi_k) chi_j'(xi_k)
    !   C_xi(i,j,p) = sum_k w_k chi_i(xi_k) chi_j(xi_k) / (xi_k^2 - 1)
    !------------------------------------

    ! Accumulated element by element in weak form; the quadrature sum runs over
    ! all local nodes (including the excluded ximax endpoint), while the basis
    ! indices skip excluded endpoints. The reference derivative matrix lives on
    ! [-0.5, 0.5], so physical derivatives carry a factor 1/elementSize.
    !
    ! p = 0 (even |m|): chi_i = L_i, the Lagrange polynomials of the element;
    !   C_xi = diag(w_i / (xi_i^2 - 1)) without the axis point (excluded for m /= 0).
    ! p = 1 (odd |m|): chi_i = t M_i / t_i with t = sqrt(xi^2 - 1), the axis factor
    !   of odd-m fields. In the first element M_i are the Lagrange polynomials of
    !   the nodes xi > 1 only, so that g = f / t stays free at the axis (with the
    !   L_i of all nodes, g(1) = 0 would be forced); elsewhere M_i = L_i. At the
    !   nodes t chi_i' = (xi M_i + (xi^2 - 1) M_i') / t_i. The axis node has no
    !   basis function; it contributes the rank-one term w_1 M_i(1) M_j(1) / (t_i t_j)
    !   to C_xi, while mass and potential terms keep their diagonal DVR form (they
    !   carry t^2, which vanishes there).
    allocate (Grid_Prolate_xiKinMatrix(Grid_Prolate_nXi, Grid_Prolate_nXi, 0:1))
    allocate (Grid_Prolate_xiCentrifugalMatrix(Grid_Prolate_nXi, Grid_Prolate_nXi, 0:1))
    Grid_Prolate_xiKinMatrix = 0.0_R64
    Grid_Prolate_xiCentrifugalMatrix = 0.0_R64

    allocate (elementPoints(Grid_Prolate_Fedvr_nLocals))
    allocate (elementWeights(Grid_Prolate_Fedvr_nLocals))
    allocate (gOdd(Grid_Prolate_Fedvr_nLocals, Grid_Prolate_Fedvr_nLocals))
    allocate (vOdd(Grid_Prolate_Fedvr_nLocals, Grid_Prolate_Fedvr_nLocals))

    do iE = 1, Grid_Prolate_Fedvr_nElements
      associate (element => Grid_Prolate_Fedvr_fedvrCtx % elements(iE), &
                 dRef => Grid_Prolate_Fedvr_derivativeCtx % firstOrderMatrix, &
                 nLoc => Grid_Prolate_Fedvr_nLocals)

        call gauss_legendre_lobatto(elementPoints, elementWeights, [element % leftBound, element % rightBound])

        iLocalStart = 1
        iLocalEnd = nLoc
        if (iE .eq. Grid_Prolate_Fedvr_nElements) iLocalEnd = nLoc - 1  ! ximax excluded

        ! odd basis at the nodes: gOdd(k,i) = t chi_i'(xi_k), vOdd(k,i) = chi_i(xi_k) / t(xi_k)
        gOdd = 0.0_R64
        vOdd = 0.0_R64
        do iLocal1 = iLocalStart, iLocalEnd
          tNode = sqrt(max(elementPoints(iLocal1)**2 - 1.0_R64, 0.0_R64))
          if (.not. (tNode > 0.0_R64)) cycle ! axis node
          do kLocal = 1, nLoc
            if (iE .eq. 1) then
              mVal = LagrangeValueSkipFirst(elementPoints, iLocal1, elementPoints(kLocal))
              mDer = LagrangeDerivativeSkipFirst(elementPoints, iLocal1, elementPoints(kLocal))
            else
              mVal = 0.0_R64
              if (kLocal .eq. iLocal1) mVal = 1.0_R64
              mDer = dRef(kLocal, iLocal1) / element % size
            end if
            gOdd(kLocal, iLocal1) = (elementPoints(kLocal) * mVal + (elementPoints(kLocal)**2 - 1.0_R64) * mDer) / tNode
            vOdd(kLocal, iLocal1) = mVal / tNode
          end do
        end do

        do iLocal2 = iLocalStart, iLocalEnd
          iGlobal2 = element % iStart + (iLocal2 - iLocalStart)
          do iLocal1 = iLocalStart, iLocalEnd
            iGlobal1 = element % iStart + (iLocal1 - iLocalStart)

            do kLocal = 1, nLoc
              Grid_Prolate_xiKinMatrix(iGlobal1, iGlobal2, 0) = Grid_Prolate_xiKinMatrix(iGlobal1, iGlobal2, 0) + &
                                                                elementWeights(kLocal) * (elementPoints(kLocal)**2 - 1.0_R64) * &
                                                                dRef(kLocal, iLocal1) * dRef(kLocal, iLocal2) / element % size**2
              Grid_Prolate_xiKinMatrix(iGlobal1, iGlobal2, 1) = Grid_Prolate_xiKinMatrix(iGlobal1, iGlobal2, 1) + &
                                                                elementWeights(kLocal) * gOdd(kLocal, iLocal1) * gOdd(kLocal, iLocal2)
              Grid_Prolate_xiCentrifugalMatrix(iGlobal1, iGlobal2, 1) = &
                Grid_Prolate_xiCentrifugalMatrix(iGlobal1, iGlobal2, 1) + &
                elementWeights(kLocal) * vOdd(kLocal, iLocal1) * vOdd(kLocal, iLocal2)
            end do

          end do
        end do

      end associate
    end do

    do iGlobal1 = 2, Grid_Prolate_nXi
      Grid_Prolate_xiCentrifugalMatrix(iGlobal1, iGlobal1, 0) = &
        Grid_Prolate_xiWeights(iGlobal1) / (Grid_Prolate_xiPoints(iGlobal1)**2 - 1.0_R64)
    end do

    deallocate (elementPoints, elementWeights, gOdd, vOdd)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Lagrange polynomial M_i of the nodes x(2:n) (node x(1) left out), evaluated at y.
  pure function LagrangeValueSkipFirst(x, i, y) result(res)
    real(R64), intent(in) :: x(:)
    integer(I32), intent(in) :: i
    real(R64), intent(in) :: y
    real(R64) :: res
    integer(I32) :: j

    res = 1.0_R64
    do j = 2, size(x)
      if (j .ne. i) res = res * (y - x(j)) / (x(i) - x(j))
    end do
  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Derivative of LagrangeValueSkipFirst with respect to y.
  pure function LagrangeDerivativeSkipFirst(x, i, y) result(res)
    real(R64), intent(in) :: x(:)
    integer(I32), intent(in) :: i
    real(R64), intent(in) :: y
    real(R64) :: res, term
    integer(I32) :: j, l

    res = 0.0_R64
    do l = 2, size(x)
      if (l .eq. i) cycle
      term = 1.0_R64 / (x(i) - x(l))
      do j = 2, size(x)
        if (j .ne. i .and. j .ne. l) term = term * (y - x(j)) / (x(i) - x(j))
      end do
      res = res + term
    end do
  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end submodule
