! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Full FEDVR-ECS implementation of the radial Poisson solver.
submodule(M_SysInteraction_Ylm_Coulomb_FullEqEcs) S_SysInteraction_Ylm_Coulomb_FullEqEcs

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Bind the full-equation ECS Coulomb solver, validate grid requirements.
  module subroutine SysInteraction_Ylm_Coulomb_FullEqEcs_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysInteraction
    use M_SysInteraction_Ylm

    call Say_Fabricate("sysInteraction.ylm.coulomb.fullEqEcs")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysInteraction_Ylm_FillInteractionPotentialRadial => FillInteractionPotentialRadial

    !------------------------------------
    ! validate grid requirements
    !------------------------------------

    if (.not. Json_GetExistence("grid.ylm.fedvrEcs")) then
      error stop "grid.ylm.fedvrEcs is required for sysInteraction.ylm.coulomb.fullEqEcs"
    end if

  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Solve the radial Poisson equation on the ECS contour via FEDVR weak form.
  !>
  !> Assembles the full complex FEDVR matrix for the operator
  !> (d²/dz² - l(l+1)/z²) in weak form along the contour, applies boundary
  !> conditions, and solves via complex LU decomposition.
  !>
  !> **Boundary conditions:**
  !>   - z = 0: Dirichlet u(0) = 0 (regularity)
  !>   - contour end: Robin dV/dz + (l+1)/z × V = 0 (asymptotic decay ~ z^-(l+1))
  !>
  !> @param[out] potLm   Radial potential component Vₗₘ(z) on the contour
  !> @param[in]  srcLm   Radial source component ρₗₘ(z) (includes complex metric weights)
  !> @param[in]  l       Angular momentum quantum number
  !> @param[in]  m       Magnetic quantum number (unused)
  !> @param[in]  time    Physical time (unused)
  !> @param[in]  bt1_    Target body type (unused)
  !> @param[in]  bt2_    Source body type (unused)
  subroutine FillInteractionPotentialRadial(potLm, srcLm, l, m, time, bt1_, bt2_)
    use M_Utils_Constants, only: PI
    use M_Utils_UnusedVariables
    use M_Utils_LapackLib
    use M_Utils_FedvrEcs
    use M_Grid_Ylm
    use M_Grid_Ylm_FedvrEcs, fedvrEcsCtx => Grid_Ylm_FedvrEcs_fedvrEcsCtx, &
      derivativeCtx => Grid_Ylm_FedvrEcs_derivativeCtx, &
      contourPoints => Grid_Ylm_FedvrEcs_contourPoints
    use M_SysInteraction_Ylm_Coulomb

    complex(R64), intent(out), contiguous :: potLm(:)
    complex(R64), intent(in), contiguous  :: srcLm(:)
    integer(I32), intent(in)  :: l
    integer(I32), intent(in)  :: m
    real(R64), intent(in)     :: time
    integer(I32), intent(in), optional  :: bt1_
    integer(I32), intent(in), optional  :: bt2_

    integer(I32) :: iRad
    complex(R64), allocatable :: poissonMatrix(:, :), rhsInSolOut(:)
    integer(I32), allocatable :: ipiv(:)
    integer(I32) :: nRad, iE
    integer(I32) :: nLoc, nE, iStart, iEnd
    integer(I32) :: iLocalStart, iLocalEnd
    complex(R64) :: z, zEnd, lTerm
    complex(R64) :: dx
    real(R64) :: strength

    if (.false.) call UnusedVariables_Mark(m, time, bt1_, bt2_)

    nE = Grid_Ylm_FedvrEcs_nElements
    nRad = Grid_Ylm_nRadial
    nLoc = Grid_Ylm_FedvrEcs_nLocals
    strength = SysInteraction_Ylm_Coulomb_Strength

    allocate (poissonMatrix(nRad, nRad))
    allocate (rhsInSolOut(nRad))
    allocate (ipiv(nRad))

    poissonMatrix = (0.0_R64, 0.0_R64)

    ! Assemble weak-form Laplacian from element matrices along the contour
    do iE = 1, nE
      associate (element => fedvrEcsCtx % elements(iE))
        dx = element % size
        iStart = element % iStart
        iEnd = element % iEnd

        iLocalStart = 1
        iLocalEnd = nLoc
        if (iE .eq. 1) iLocalStart = 2

        poissonMatrix(iStart:iEnd, iStart:iEnd) = poissonMatrix(iStart:iEnd, iStart:iEnd) + &
                                                  derivativeCtx % secondOrderMatrix(iLocalStart:iLocalEnd, &
                                                                                    iLocalStart:iLocalEnd) / dx
      end associate
    end do

    ! Apply z-weighting for spherical coordinates on the contour
    do iRad = 1, nRad
      z = contourPoints(iRad)
      poissonMatrix(:, iRad) = poissonMatrix(:, iRad) * z
      poissonMatrix(iRad, :) = poissonMatrix(iRad, :) * z
    end do

    ! Add centrifugal term -l(l+1)/z² in weak form (complex contour weights)
    if (l > 0) then
      lTerm = cmplx(-l * (l + 1), 0.0_R64, kind=R64)
      do iRad = 1, nRad
        poissonMatrix(iRad, iRad) = poissonMatrix(iRad, iRad) + lTerm * fedvrEcsCtx % weights(iRad)
      end do
    end if

    ! Set up RHS: -4π ρ (fully complex solve)
    rhsInSolOut(:) = -4.0_R64 * PI * srcLm(:)

    ! Apply Robin BC at the contour end: dV/dz + (l+1)/z × V = 0
    associate (element => fedvrEcsCtx % elements(nE))
      dx = element % size
      iStart = element % iStart
      iEnd = element % iEnd

      iLocalStart = 1
      iLocalEnd = nLoc
      if (nE .eq. 1) iLocalStart = 2

      zEnd = contourPoints(iEnd)

      poissonMatrix(iEnd, :) = (0.0_R64, 0.0_R64)
      rhsInSolOut(iEnd) = (0.0_R64, 0.0_R64)

      poissonMatrix(iEnd, iStart:iEnd) = derivativeCtx % firstOrderMatrix(iLocalEnd, &
                                                                          iLocalStart:iLocalEnd) / dx
      poissonMatrix(iEnd, iEnd) = poissonMatrix(iEnd, iEnd) + (l + 1) / zEnd
    end associate

    ! Solve via complex LU factorization
    call LapackLib_FactorizeLU(poissonMatrix, ipiv)
    call LapackLib_SolveFactorized(poissonMatrix, ipiv, rhsInSolOut, 'N')

    potLm(:) = strength * rhsInSolOut(:)

    deallocate (poissonMatrix, rhsInSolOut, ipiv)
  end subroutine
end submodule
