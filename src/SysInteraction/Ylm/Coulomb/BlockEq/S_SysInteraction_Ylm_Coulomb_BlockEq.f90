! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Block-factorized FEDVR implementation with Schur complement.
!>
!> @details Precomputes element-level LU factorizations and interface responses
!> during setup. Each solve then requires only:
!>   1. Map source to element-local RHS and solve (using cached LU)
!>   2. Assemble and solve small Schur complement system at interfaces
!>   3. Back-substitute to recover full solution
submodule(M_SysInteraction_Ylm_Coulomb_BlockEq) S_SysInteraction_Ylm_Coulomb_BlockEq

  implicit none

  !> @brief Storage for precomputed LU factorization of element matrices.
  type :: T_SubBlock
    complex(R64), allocatable :: matrix(:, :)  !< LU-factorized element matrix
    integer(I32), allocatable :: ipiv(:)       !< Pivot indices
  end type

  !> Precomputed element factorizations for each (element, l) pair
  type(T_SubBlock), allocatable :: subBlocks(:, :)            ! (nE, lmax+1)

  !> Unit response to left boundary excitation
  complex(R64), allocatable :: unitResponseLeft(:, :, :)      ! (nLoc, nE, lmax+1)

  !> Unit response to right boundary excitation
  complex(R64), allocatable :: unitResponseRight(:, :, :)     ! (nLoc, nE, lmax+1)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Bind the block-equation Coulomb solver.
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
  !> @brief Precompute element factorizations and unit responses.
  !>
  !> For each (element, l) pair:
  !>   1. Build and LU-factorize the element Poisson matrix
  !>   2. Compute unit response to left boundary: solve with RHS = [1,0,...,0]
  !>   3. Compute unit response to right boundary: solve with RHS = [0,...,0,1]
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
    allocate (unitResponseLeft(nLoc, nE, lmaxPot + 1))
    allocate (unitResponseRight(nLoc, nE, lmaxPot + 1))

    do iL = 1, lmaxPot + 1
      lVal = iL - 1
      do iE = 1, nE
        call BuildSubBlock(subBlocks(iE, iL), fedvrCtx % elements(iE), lVal, (iE .eq. 1), (iE .eq. nE))

        ! Left boundary unit response: excite left DOF
        unitResponseLeft(:, iE, iL) = cmplx(0.0_R64, 0.0_R64, kind=R64)
        unitResponseLeft(1, iE, iL) = cmplx(1.0_R64, 0.0_R64, kind=R64)
        call SolveElementSystem(subBlocks(iE, iL) % matrix, unitResponseLeft(:, iE, iL), subBlocks(iE, iL) % ipiv)

        ! Right boundary unit response: excite right DOF
        unitResponseRight(:, iE, iL) = cmplx(0.0_R64, 0.0_R64, kind=R64)
        unitResponseRight(nLoc, iE, iL) = cmplx(1.0_R64, 0.0_R64, kind=R64)
        call SolveElementSystem(subBlocks(iE, iL) % matrix, unitResponseRight(:, iE, iL), subBlocks(iE, iL) % ipiv)
      end do
    end do
  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Solve radial Poisson equation using precomputed block factorizations.
  !>
  !> **Algorithm:**
  !>   1. For each element, solve with zero boundary conditions (homogeneous)
  !>   2. Build Schur complement system for interface values
  !>   3. Solve interface system
  !>   4. Recover full solution by superposition
  !>
  !> @param[out] potLm   Radial potential component Vₗₘ(r)
  !> @param[in]  srcLm   Radial density component ρₗₘ(r) (includes weights)
  !> @param[in]  l       Angular momentum quantum number
  !> @param[in]  m       Magnetic quantum number (unused)
  !> @param[in]  time    Physical time (unused)
  !> @param[in]  bt1_    Target body type (unused)
  !> @param[in]  bt2_    Source body type (unused)
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

    complex(R64), allocatable :: schurMatrix(:, :)
    complex(R64), allocatable :: schurWeightRight(:)
    complex(R64), allocatable :: rhsRight(:)
    ! Workspace for homogeneous element solutions; local (not module state) so
    ! that concurrent solves for different (l,m) channels stay thread-safe
    complex(R64), allocatable :: homogeneousResponse(:, :)
    integer(I32) :: nRad, nE, nLoc
    real(R64)    :: strength

    if (.false.) call UnusedVariables_Mark(m, time, bt1_, bt2_)

    nE = Grid_Ylm_Fedvr_nElements
    nRad = Grid_Ylm_nRadial
    nLoc = Grid_Ylm_Fedvr_nLocals
    strength = SysInteraction_Ylm_Coulomb_Strength

    allocate (homogeneousResponse(nLoc, nE))
    homogeneousResponse = cmplx(0.0_R64, 0.0_R64, kind=R64)

    allocate (rhsRight(nE - 1))
    allocate (schurMatrix(nE - 1, nE - 1))
    allocate (schurWeightRight(nE - 1))

    rhsRight = cmplx(0.0_R64, 0.0_R64, kind=R64)
    schurMatrix = cmplx(0.0_R64, 0.0_R64, kind=R64)
    schurWeightRight = cmplx(0.0_R64, 0.0_R64, kind=R64)

    ! Step 1: Solve element systems with homogeneous BCs
    call FillHomogeneousResponse(homogeneousResponse, rhsRight, srcLm, l)

    ! Step 2: Build Schur complement system
    call BuildInterfaceSystem(schurMatrix, schurWeightRight, homogeneousResponse, &
                              unitResponseLeft(:, :, l + 1), unitResponseRight(:, :, l + 1), &
                              rhsRight)

    ! Step 3: Solve interface system
    if (nE >= 2) then
      call SolveInterfaceSystem(schurMatrix, schurWeightRight)
    end if

    ! Step 4: Recover full solution
    call RecoverFullSolution(potLm, schurWeightRight, &
                             homogeneousResponse, &
                             unitResponseLeft(:, :, l + 1), unitResponseRight(:, :, l + 1), &
                             rhsRight)

    potLm(:) = strength * potLm(:)

    deallocate (rhsRight, schurMatrix, schurWeightRight)
  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Build and factorize element Poisson matrix.
  !>
  !> Constructs the FEDVR weak-form matrix for (d²/dr² - l(l+1)/r²) on one
  !> element, applies boundary modifications, and LU-factorizes for fast solves.
  !>
  !> @param[inout] subBlock  Output structure with matrix and pivots
  !> @param[in]    element   FEDVR element data
  !> @param[in]    l         Angular momentum quantum number
  !> @param[in]    firstQ    True if this is the first element (r=0 boundary)
  !> @param[in]    lastQ     True if this is the last element (r_max boundary)
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

    ! Start with second-derivative matrix scaled by element size
    subBlock % matrix(:, :) = derivativeCtx % secondOrderMatrix(:, :) / dx

    ! Apply spherical r-weighting
    do iRad = element % iStart, element % iEnd
      r = fedvrCtx % points(iRad)
      iLoc = iRad - element % iStart + 1
      if (firstQ) iLoc = iLoc + 1
      subBlock % matrix(iLoc, :) = subBlock % matrix(iLoc, :) * r
      subBlock % matrix(:, iLoc) = subBlock % matrix(:, iLoc) * r
    end do

    ! Add centrifugal term
    if (l > 0) then
      lTerm = real(-l * (l + 1), kind=R64)
      do iRad = element % iStart, element % iEnd
        iLoc = iRad - element % iStart + 1
        if (firstQ) iLoc = iLoc + 1
        subBlock % matrix(iLoc, iLoc) = subBlock % matrix(iLoc, iLoc) + &
                                        lTerm * element % weights(iLoc)
      end do
    end if

    ! Regularization shift to break degeneracy
    subBlock % matrix(1, 1) = subBlock % matrix(1, 1) + 3.0_R64
    subBlock % matrix(nLoc, nLoc) = subBlock % matrix(nLoc, nLoc) - 3.0_R64

    ! Dirichlet BC at r=0
    if (firstQ) then
      subBlock % matrix(:, 1) = 0.0_R64
      subBlock % matrix(1, :) = 0.0_R64
      subBlock % matrix(1, 1) = 1.0_R64
    end if

    ! Robin BC at r_max
    if (lastQ) then
      subBlock % matrix(nLoc, :) = derivativeCtx % firstOrderMatrix(nLoc, :) / dx
      subBlock % matrix(nLoc, nLoc) = subBlock % matrix(nLoc, nLoc) + (l + 1) / &
                                      fedvrCtx % points(element % iEnd)
    end if

    call LapackLib_FactorizeLU(subBlock % matrix, subBlock % ipiv)
  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Solve element system using precomputed LU factorization.
  subroutine SolveElementSystem(matrix, rhsInSolOut, ipiv)
    use M_Utils_LapackLib

    complex(R64), intent(inout), contiguous :: matrix(:, :)
    complex(R64), intent(inout), contiguous :: rhsInSolOut(:)
    integer(I32), intent(in), contiguous :: ipiv(:)

    call LapackLib_SolveFactorized(matrix, ipiv, rhsInSolOut, 'N')

  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Solve element systems with homogeneous boundary conditions.
  !>
  !> For each element, sets up the RHS from the source term, imposes zero
  !> boundary values, and solves using the precomputed LU factorization.
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

        ! Impose homogeneous BCs at element boundaries
        homogeneousResponse(1, iE) = cmplx(0.0_R64, 0.0_R64, kind=R64)
        homogeneousResponse(nLoc, iE) = cmplx(0.0_R64, 0.0_R64, kind=R64)

        call SolveElementSystem(subBlocks(iE, l + 1) % matrix, homogeneousResponse(:, iE), &
                                subBlocks(iE, l + 1) % ipiv)
      end associate
    end do
  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Build Schur complement system for interface coupling.
  !>
  !> The interface values must satisfy continuity constraints between elements.
  !> This builds a tridiagonal system relating adjacent interface DOFs.
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

    do iE = 1, nE - 1
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
  !> @brief Solve the interface system via LU factorization.
  subroutine SolveInterfaceSystem(matrix, rhsInSolOut)
    use M_Utils_LapackLib

    complex(R64), intent(inout), contiguous :: matrix(:, :)
    complex(R64), intent(inout), contiguous :: rhsInSolOut(:)

    integer(I32), allocatable :: ipiv(:)

    allocate (ipiv(size(rhsInSolOut)))

    call LapackLib_FactorizeLU(matrix, ipiv)
    call LapackLib_SolveFactorized(matrix, ipiv, rhsInSolOut, 'N')

    deallocate (ipiv)
  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Recover full solution from interface values via superposition.
  !>
  !> Full solution = homogeneous + left_response × left_value + right_response × right_value
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

    potLm(:) = cmplx(0.0_R64, 0.0_R64, kind=R64)

    do iE = 1, nE
      associate (element => fedvrCtx % elements(iE))

        iStart = element % iStart
        iEnd = element % iEnd
        iLocalStart = 1
        iLocalEnd = nLoc
        if (iE .eq. 1) iLocalStart = 2

        potLm(iStart:iEnd) = homogeneousResponse(iLocalStart:iLocalEnd, iE)

        if (iE < nE) then
          potLm(iStart:iEnd) = potLm(iStart:iEnd) + unitResponseRight(iLocalStart:iLocalEnd, iE) * schurWeightRight(iE)
        end if

        if (iE > 1) then
          lE = iE - 1
          schurWeightLeft = rhsRight(lE) - schurWeightRight(lE)
          potLm(iStart:iEnd) = potLm(iStart:iEnd) + unitResponseLeft(iLocalStart:iLocalEnd, iE) * schurWeightLeft
        end if
      end associate
    end do
  end subroutine

end submodule
