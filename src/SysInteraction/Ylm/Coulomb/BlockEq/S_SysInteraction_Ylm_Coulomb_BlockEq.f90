! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysInteraction_Ylm_Coulomb_BlockEq) S_SysInteraction_Ylm_Coulomb_BlockEq

  implicit none

  !> Type for complex element LU factorization data
  type :: T_SubBlock
    !> Complex element matrix
    complex(R64), allocatable :: matrix(:, :)
    !> Pivot indices from LU factorization
    integer(I32), allocatable :: ipiv(:)
  end type

  type(T_SubBlock), allocatable :: subBlocks(:, :)            ! (nE, lmax+1)
  complex(R64), allocatable :: homogeneousResponse(:, :)
  complex(R64), allocatable :: unitResponseLeft(:, :, :)      ! (nLoc, nE, lmax+1)
  complex(R64), allocatable :: unitResponseRight(:, :, :)     ! (nLoc, nE, lmax+1)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysInteraction_Ylm_Coulomb_BlockEq_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysInteraction
    use M_SysInteraction_Ylm

    call Say_Fabricate("sysInteraction.ylm.coulomb.BlockEq")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysInteraction_Ylm_FillInteractionPotentialRadial => FillInteractionPotentialRadial
    SysInteraction_Setup => Setup

  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Setup
    use M_Utils_Fedvr
    use M_SysInteraction_Ylm
    use M_Grid_Ylm
    use M_Grid_Ylm_Fedvr, fedvrCtx => Grid_Ylm_Fedvr_fedvrCtx

    integer(I32) :: iE, iL
    integer(I32) :: nE
    integer(I32) :: lVal, lmaxPot
    integer(I32) :: nLoc

    nE = Grid_Ylm_Fedvr_nElements
    lmaxPot = SysInteraction_Ylm_lmax
    nLoc = Grid_Ylm_Fedvr_nLocals

    allocate (subBlocks(nE, lmaxPot + 1))
    allocate (homogeneousResponse(nLoc, nE))
    allocate (unitResponseLeft(nLoc, nE, lmaxPot + 1))
    allocate (unitResponseRight(nLoc, nE, lmaxPot + 1))

    do iL = 1, lmaxPot + 1
      lVal = iL - 1
      do iE = 1, nE
        ! Build and factorize element matrix for (l, element) using derivativeCtx
        call BuildSubBlock(subBlocks(iE, iL), fedvrCtx % elements(iE), lVal, (iE .eq. 1), (iE .eq. nE))

        ! Left boundary unit response
        unitResponseLeft(:, iE, iL) = cmplx(0.0_R64, 0.0_R64, kind=R64)
        unitResponseLeft(1, iE, iL) = cmplx(1.0_R64, 0.0_R64, kind=R64)
        call SolveElementSystem(subBlocks(iE, iL) % matrix, unitResponseLeft(:, iE, iL), subBlocks(iE, iL) % ipiv)

        ! Right boundary unit response
        unitResponseRight(:, iE, iL) = cmplx(0.0_R64, 0.0_R64, kind=R64)
        unitResponseRight(nLoc, iE, iL) = cmplx(1.0_R64, 0.0_R64, kind=R64)
        call SolveElementSystem(subBlocks(iE, iL) % matrix, unitResponseRight(:, iE, iL), subBlocks(iE, iL) % ipiv)
      end do
    end do
  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Solves the radial multipole Poisson equation using block LU decomposition
  !> (d²/dr² - ℓ(ℓ+1)/r²) u(r) = -4πrg(r)
  !> @param potLm Complex output potential for the (l,m) component
  !> @param srcLm Complex input density for the (l,m) component
  !> @param l Angular momentum quantum number
  subroutine FillInteractionPotentialRadial(potLm, srcLm, l, m, time, bt1_, bt2_)
    use M_Utils_UnusedVariables
    use M_Utils_Fedvr
    use M_Grid_Ylm
    use M_Grid_Ylm_Fedvr, fedvrCtx => Grid_Ylm_Fedvr_fedvrCtx, derivativeCtx => Grid_Ylm_Fedvr_derivativeCtx
    use M_SysInteraction_Ylm_Coulomb

    complex(R64), intent(out), contiguous :: potLm(:)
    complex(R64), intent(in), contiguous  :: srcLm(:)
    integer(I32), intent(in)  :: l
    integer(I32), intent(in)  :: m
    real(R64), intent(in)     :: time
    integer(I32), intent(in), optional  :: bt1_
    integer(I32), intent(in), optional  :: bt2_

    ! Local variables
    complex(R64), allocatable :: schurMatrix(:, :)
    complex(R64), allocatable :: schurWeightRight(:)
    complex(R64), allocatable :: rhsRight(:)
    integer(I32) :: nRad, nE, nLoc
    real(R64)    :: strength

    if (.false.) call UnusedVariables_Mark(m, time, bt1_, bt2_)

    ! Get dimensions
    nE = Grid_Ylm_Fedvr_nElements
    nRad = Grid_Ylm_nRadial
    nLoc = Grid_Ylm_Fedvr_nLocals
    strength = SysInteraction_Ylm_Coulomb_Strength

    ! Initialize homogeneousResponse with source terms
    homogeneousResponse = cmplx(0.0_R64, 0.0_R64, kind=R64)

    ! Allocate temporary arrays
    allocate (rhsRight(nE - 1))
    allocate (schurMatrix(nE - 1, nE - 1))
    allocate (schurWeightRight(nE - 1))

    ! Initialize arrays
    rhsRight = cmplx(0.0_R64, 0.0_R64, kind=R64)
    schurMatrix = cmplx(0.0_R64, 0.0_R64, kind=R64)
    schurWeightRight = cmplx(0.0_R64, 0.0_R64, kind=R64)

    ! Map source term to element-local right-hand sides and solve element systems
    call FillHomogeneousResponse(homogeneousResponse, rhsRight, srcLm, l)

    ! Build Schur complement system
    call BuildInterfaceSystem(schurMatrix, schurWeightRight, homogeneousResponse, &
                              unitResponseLeft(:, :, l + 1), unitResponseRight(:, :, l + 1), &
                              rhsRight)

    ! Solve the interface system
    if (nE >= 2) then
      call SolveInterfaceSystem(schurMatrix, schurWeightRight)
    end if

    ! Recover full solution from interface values
    call RecoverFullSolution(potLm, schurWeightRight, &
                             homogeneousResponse, &
                             unitResponseLeft(:, :, l + 1), unitResponseRight(:, :, l + 1), &
                             rhsRight)

    potLm(:) = strength * potLm(:)

    ! Clean up
    deallocate (rhsRight, schurMatrix, schurWeightRight)
  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Builds the element matrix for the Poisson operator (d²/dr² - ℓ(ℓ+1)/r²)
  !> and performs LU factorization
  !> @param subBlock Output structure containing matrix and pivots
  !> @param element FEDVR element data
  !> @param l Angular momentum quantum number
  !> @param firstQ Flag indicating if this is the first element
  !> @param lastQ Flag indicating if this is the last element
  subroutine BuildSubBlock(subBlock, element, l, firstQ, lastQ)
    use M_Utils_LapackLib
    use M_Utils_Fedvr
    use M_Grid_Ylm_Fedvr, fedvrCtx => Grid_Ylm_Fedvr_fedvrCtx, derivativeCtx => Grid_Ylm_Fedvr_derivativeCtx

    type(T_SubBlock), intent(inout) :: subBlock
    type(T_Fedvr_Element), intent(in) :: element
    integer(I32), intent(in) :: l
    logical :: firstQ, lastQ

    integer(I32) :: iLoc, nLoc, nE, iRad
    real(R64) :: r, lTerm, dx

    nE = Grid_Ylm_Fedvr_nElements
    nLoc = Grid_Ylm_Fedvr_nLocals
    dx = element % size

    allocate (subBlock % matrix(nLoc, nLoc), subBlock % ipiv(nLoc))

    ! Assemble weak-form second-order operator from reference matrix:
    ! K_e = (1/dx) * secondOrderMatrix_ref
    subBlock % matrix(:, :) = derivativeCtx % secondOrderMatrix(:, :) / dx

    do iRad = element % iStart, element % iEnd
      r = fedvrCtx % points(iRad)
      iLoc = iRad - element % iStart + 1
      if (firstQ) iLoc = iLoc + 1
      subBlock % matrix(iLoc, :) = subBlock % matrix(iLoc, :) * r
      subBlock % matrix(:, iLoc) = subBlock % matrix(:, iLoc) * r
    end do

    ! Add centrifugal term on interior rows only (skip physical boundary rows)
    if (l > 0) then
      lTerm = real(-l * (l + 1), kind=R64) ! the devision by / (r * r) is removed by multiplication with r
      do iRad = element % iStart, element % iEnd
        iLoc = iRad - element % iStart + 1
        if (firstQ) iLoc = iLoc + 1
        subBlock % matrix(iLoc, iLoc) = subBlock % matrix(iLoc, iLoc) + &
                                        lTerm * element % weights(iLoc)
      end do
    end if

    ! since the decomposition of the global matrix into block components is not
    ! unique, we add a shift to the first and last diagonal element to break the singularity
    ! of the second order derivative matrix (sum of each row is zero) so [1..1] is a null vector.
    ! The natural boundary condition of the second order derivative matrix seems to be Neumann
    ! which is degenerate due to the constant null vector.
    subBlock % matrix(1, 1) = subBlock % matrix(1, 1) + 3.0_R64
    subBlock % matrix(nLoc, nLoc) = subBlock % matrix(nLoc, nLoc) - 3.0_R64

    if (firstQ) then
      subBlock % matrix(:, 1) = 0.0_R64
      subBlock % matrix(1, :) = 0.0_R64
      subBlock % matrix(1, 1) = 1.0_R64
    end if

    if (lastQ) then
      subBlock % matrix(nLoc, :) = derivativeCtx % firstOrderMatrix(nLoc, :) / dx
      subBlock % matrix(nLoc, nLoc) = subBlock % matrix(nLoc, nLoc) + (l + 1) / &
                                      fedvrCtx % points(element % iEnd)
    end if
    call LapackLib_FactorizeLU(subBlock % matrix, subBlock % ipiv)
  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !> Solves a linear system using the precomputed LU factorization for an element
  !> @param matrix LU-factorized element matrix
  !> @param rhsInSolOut Input right-hand side vector, overwritten with solution
  !> @param ipiv Pivot indices from LU factorization
  subroutine SolveElementSystem(matrix, rhsInSolOut, ipiv)
    use M_Utils_LapackLib

    complex(R64), intent(inout), contiguous :: matrix(:, :)
    complex(R64), intent(inout), contiguous :: rhsInSolOut(:)
    integer(I32), intent(in), contiguous :: ipiv(:)

    ! Solve system using generic LAPACK wrapper
    call LapackLib_SolveFactorized(matrix, ipiv, rhsInSolOut, 'N')

  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Maps source term to element-local right-hand sides and solves element systems
  !> @param homogeneousResponse Output array for homogeneous solutions
  !> @param rhsRight Interface values of right-hand side
  !> @param srcLm Input density for the (l,m) component
  !> @param l Angular momentum quantum number
  subroutine FillHomogeneousResponse(homogeneousResponse, rhsRight, srcLm, l)
    use M_Utils_Constants, only: PI
    use M_Grid_Ylm_Fedvr, fedvrCtx => Grid_Ylm_Fedvr_fedvrCtx

    complex(R64), intent(out), contiguous :: homogeneousResponse(:, :)
    complex(R64), intent(out), contiguous :: rhsRight(:)
    complex(R64), intent(in), contiguous :: srcLm(:)
    integer(I32), intent(in) :: l

    integer(I32) :: nE, nLoc
    integer(I32) :: iE, iStart, iEnd
    integer(I32) :: iLocalStart, iLocalEnd

    nE = Grid_Ylm_Fedvr_nElements
    nLoc = Grid_Ylm_Fedvr_nLocals

    do iE = 1, nE
      associate (element => fedvrCtx % elements(iE))

        iStart = element % iStart
        iEnd = element % iEnd

        iLocalStart = 1
        iLocalEnd = nLoc
        if (iE .eq. 1) iLocalStart = 2

        homogeneousResponse(iLocalStart:iLocalEnd, iE) = -4.0_R64 * PI * srcLm(iStart:iEnd)

        if (iE < nE) rhsRight(iE) = homogeneousResponse(nLoc, iE)

        homogeneousResponse(1, iE) = cmplx(0.0_R64, 0.0_R64, kind=R64)
        homogeneousResponse(nLoc, iE) = cmplx(0.0_R64, 0.0_R64, kind=R64)

        call SolveElementSystem(subBlocks(iE, l + 1) % matrix, homogeneousResponse(:, iE), &
                                subBlocks(iE, l + 1) % ipiv)
      end associate
    end do
  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Builds the Schur complement system for interface coupling between elements
  !> @param schurMatrix Output Schur complement matrix
  !> @param schurWeightRight Output right-hand side vector for Schur system
  !> @param homogeneousResponse Homogeneous response from element solves
  !> @param unitResponseLeft Unit response vectors for left boundary conditions
  !> @param unitResponseRight Unit response vectors for right boundary conditions
  !> @param rhsRight Interface values of right-hand side
  !>
  !> Builds the Schur complement system for interface coupling based on derivation
  !  elemHomoSol(nLoc, iE)
  ! + unitResponseRight(nLoc, iE) * r(iE)
  ! + unitResponseLeft(nLoc, iE) * l(iE)
  ! ==
  !  elemHomoSol(1, rE)
  ! + unitResponseLeft(1, rE) * l(rE)
  ! + unitResponseRight(1, rE) * r(rE)

  ! l(rE) .eq. elemRhs(nLoc, iE) - r(iE)
  ! l(iE) .eq. elemRhs(1, iE) - r(lE)

  ! !---

  !  elemHomoSol(nLoc, iE)
  ! + unitResponseRight(nLoc, iE) * r(iE)
  ! + unitResponseLeft(nLoc, iE) * (elemRhs(1, iE) - r(lE))
  ! ==
  !  elemHomoSol(1, rE)
  ! + unitResponseLeft(1, rE) * (elemRhs(nLoc, iE) - r(iE))
  ! + unitResponseRight(1, rE) * r(rE)

  ! !---

  !  elemHomoSol(nLoc, iE)
  ! - elemHomoSol(1, rE)
  ! + unitResponseLeft(nLoc, iE) * elemRhs(1, iE)
  ! - unitResponseLeft(1, rE) * elemRhs(nLoc, iE)
  ! ==
  ! - unitResponseRight(nLoc, iE) * r(iE)
  ! + unitResponseLeft(nLoc, iE) * r(lE)
  ! - unitResponseLeft(1, rE) * r(iE)
  ! + unitResponseRight(1, rE) * r(rE)
  subroutine BuildInterfaceSystem(schurMatrix, schurWeightRight, homogeneousResponse, unitResponseLeft, &
                                  unitResponseRight, rhsRight)
    use M_Grid_Ylm_Fedvr, fedvrCtx => Grid_Ylm_Fedvr_fedvrCtx

    complex(R64), intent(inout) :: schurMatrix(:, :)
    complex(R64), intent(inout) :: schurWeightRight(:)
    complex(R64), intent(in) :: homogeneousResponse(:, :)
    complex(R64), intent(in) :: unitResponseLeft(:, :)
    complex(R64), intent(in) :: unitResponseRight(:, :)
    complex(R64), intent(in) :: rhsRight(:)

    integer(I32) :: iE, lE, rE
    integer(I32) :: nE, nLoc

    nE = Grid_Ylm_Fedvr_nElements
    nLoc = Grid_Ylm_Fedvr_nLocals

    do iE = 1, nE - 1 ! iterate over interface nodes
      lE = iE - 1
      rE = iE + 1

      schurWeightRight(iE) = homogeneousResponse(nLoc, iE) &
                             - homogeneousResponse(1, rE) &
                             - unitResponseLeft(1, rE) * rhsRight(iE)

      if (iE > 1) then
        schurWeightRight(iE) = schurWeightRight(iE) &
                               + unitResponseLeft(nLoc, iE) * rhsRight(lE)
        schurMatrix(iE, lE) = unitResponseLeft(nLoc, iE)
      end if

      schurMatrix(iE, iE) = -(unitResponseRight(nLoc, iE) + &
                              unitResponseLeft(1, rE))

      if (iE < nE - 1) then
        schurMatrix(iE, rE) = unitResponseRight(1, rE)
      end if
    end do
  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Solves the interface system for boundary matching
  !> @param matrix The complex Schur complement matrix
  !> @param rhs Complex right-hand side vector
  !> @param solution Complex interface values solution
  !> @param n Size of the system
  subroutine SolveInterfaceSystem(matrix, rhsInSolOut)
    use M_Utils_LapackLib

    complex(R64), intent(inout), contiguous :: matrix(:, :)
    complex(R64), intent(inout), contiguous :: rhsInSolOut(:)

    integer(I32), allocatable :: ipiv(:)

    ! For small systems, direct solve is fine
    allocate (ipiv(size(rhsInSolOut)))

    ! LU factorization and solve using generic LAPACK wrappers
    call LapackLib_FactorizeLU(matrix, ipiv)
    call LapackLib_SolveFactorized(matrix, ipiv, rhsInSolOut, 'N')

    deallocate (ipiv)
  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Recovers full solution from interface values
  subroutine RecoverFullSolution(potLm, schurWeightRight, homogeneousResponse, unitResponseLeft, &
                                 unitResponseRight, rhsRight)
    use M_Grid_Ylm_Fedvr, fedvrCtx => Grid_Ylm_Fedvr_fedvrCtx

    complex(R64), intent(out) :: potLm(:)
    complex(R64), intent(in) :: schurWeightRight(:)
    complex(R64), intent(in) :: homogeneousResponse(:, :)
    complex(R64), intent(in) :: unitResponseLeft(:, :)
    complex(R64), intent(in) :: unitResponseRight(:, :)
    complex(R64), intent(in) :: rhsRight(:)

    integer(I32) :: nE, nLoc
    integer(I32) :: iE, lE
    integer(I32) :: iStart, iEnd
    complex(R64) :: schurWeightLeft
    integer(I32) :: iLocalStart, iLocalEnd

    nE = Grid_Ylm_Fedvr_nElements
    nLoc = Grid_Ylm_Fedvr_nLocals

    potLm(:) = cmplx(0.0_R64, 0.0_R64, kind=R64)  ! Initialize solution to zero

    do iE = 1, nE
      associate (element => fedvrCtx % elements(iE))

        iStart = element % iStart
        iEnd = element % iEnd
        iLocalStart = 1
        iLocalEnd = nLoc
        if (iE .eq. 1) iLocalStart = 2

        ! Map local solution back to global vector
        potLm(iStart:iEnd) = homogeneousResponse(iLocalStart:iLocalEnd, iE)

        ! Add contribution from right interface value
        if (iE < nE) then
          potLm(iStart:iEnd) = potLm(iStart:iEnd) + unitResponseRight(iLocalStart:iLocalEnd, iE) * schurWeightRight(iE)
        end if

        ! Add contribution from left interface value
        if (iE > 1) then
          lE = iE - 1
          schurWeightLeft = rhsRight(lE) - schurWeightRight(lE)
          potLm(iStart:iEnd) = potLm(iStart:iEnd) + unitResponseLeft(iLocalStart:iLocalEnd, iE) * schurWeightLeft
        end if
      end associate
    end do
  end subroutine

end submodule
