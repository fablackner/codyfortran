! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Full FEDVR implementation of radial Poisson solver.
submodule(M_SysInteraction_Ylm_Coulomb_FullEq) S_SysInteraction_Ylm_Coulomb_FullEq

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Bind the full-equation FEDVR Coulomb solver.
  module subroutine SysInteraction_Ylm_Coulomb_FullEq_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysInteraction
    use M_SysInteraction_Ylm

    call Say_Fabricate("sysInteraction.ylm.coulomb.fullEq")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysInteraction_Ylm_FillInteractionPotentialRadial => FillInteractionPotentialRadial

  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Solve radial Poisson equation via FEDVR weak-form discretization.
  !>
  !> Assembles the full FEDVR matrix for the operator (d²/dr² - l(l+1)/r²) in
  !> weak form, applies boundary conditions, and solves via LU decomposition.
  !>
  !> **Boundary conditions:**
  !>   - r = 0: Dirichlet u(0) = 0 (regularity)
  !>   - r = r_max: Robin du/dr + (l+1)/r × u = 0 (asymptotic decay ~ r^-(l+1))
  !>
  !> @param[out] potLm   Radial potential component Vₗₘ(r)
  !> @param[in]  srcLm   Radial density component ρₗₘ(r) (includes weights)
  !> @param[in]  l       Angular momentum quantum number
  !> @param[in]  m       Magnetic quantum number (unused)
  !> @param[in]  time    Physical time (unused)
  !> @param[in]  bt1_    Target body type (unused)
  !> @param[in]  bt2_    Source body type (unused)
  subroutine FillInteractionPotentialRadial(potLm, srcLm, l, m, time, bt1_, bt2_)
    use M_Utils_Constants, only: PI
    use M_Utils_UnusedVariables
    use M_Utils_LapackLib
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

    integer(I32) :: iRad
    real(R64), allocatable :: poissonMatrix(:, :), rhsInSolOut(:, :)
    integer(I32), allocatable :: ipiv(:)
    integer(I32) :: nRad, iE
    integer(I32) :: nLoc, nE, iStart, iEnd
    integer(I32) :: iLocalStart, iLocalEnd
    real(R64) :: r, lTerm, rEnd
    real(R64) :: dx
    real(R64) :: strength

    if (.false.) call UnusedVariables_Mark(m, time, bt1_, bt2_)

    nE = Grid_Ylm_Fedvr_nElements
    nRad = Grid_Ylm_nRadial
    nLoc = Grid_Ylm_Fedvr_nLocals
    strength = SysInteraction_Ylm_Coulomb_Strength

    allocate (poissonMatrix(nRad, nRad))
    allocate (rhsInSolOut(nRad, 2))
    allocate (ipiv(nRad))

    poissonMatrix = 0.0_R64

    ! Assemble weak-form Laplacian from element matrices
    do iE = 1, nE
      associate (element => fedvrCtx % elements(iE))
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

    ! Apply r-weighting for spherical coordinates
    do iRad = 1, nRad
      r = Grid_Ylm_radialPoints(iRad)
      poissonMatrix(:, iRad) = poissonMatrix(:, iRad) * r
      poissonMatrix(iRad, :) = poissonMatrix(iRad, :) * r
    end do

    ! Add centrifugal term -l(l+1)/r² in weak form
    if (l > 0) then
      lTerm = real(-l * (l + 1), kind=R64)
      do iRad = 1, nRad
        poissonMatrix(iRad, iRad) = poissonMatrix(iRad, iRad) + lTerm * fedvrCtx % weights(iRad)
      end do
    end if

    ! Set up RHS: -4π ρ (real and imaginary parts separately)
    rhsInSolOut(:, 1) = real(-4.0_R64 * PI * srcLm(:))
    rhsInSolOut(:, 2) = aimag(-4.0_R64 * PI * srcLm(:))

    ! Apply Robin BC at r_max: du/dr + (l+1)/r × u = 0
    associate (element => fedvrCtx % elements(nE))
      dx = element % size
      iStart = element % iStart
      iEnd = element % iEnd

      iLocalStart = 1
      iLocalEnd = nLoc
      if (nE .eq. 1) iLocalStart = 2

      rEnd = Grid_Ylm_radialPoints(iEnd)

      poissonMatrix(iEnd, :) = 0.0_R64
      rhsInSolOut(iEnd, :) = 0.0_R64

      poissonMatrix(iEnd, iStart:iEnd) = derivativeCtx % firstOrderMatrix(iLocalEnd, &
                                                                          iLocalStart:iLocalEnd) / dx
      poissonMatrix(iEnd, iEnd) = poissonMatrix(iEnd, iEnd) + (l + 1) / rEnd
    end associate

    ! Solve via LU factorization
    call LapackLib_FactorizeLU(poissonMatrix, ipiv)
    call LapackLib_SolveFactorized(poissonMatrix, ipiv, rhsInSolOut, 'N')

    potLm(:) = strength * cmplx(rhsInSolOut(:, 1), rhsInSolOut(:, 2), kind=R64)

    deallocate (poissonMatrix, rhsInSolOut, ipiv)
  end subroutine
end submodule
