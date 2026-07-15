! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Implementation of the standard prolate Coulomb xi solver.
submodule(M_SysInteraction_Prolate_Coulomb_StdImpl) S_SysInteraction_Prolate_Coulomb_StdImpl

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Bind the channel solver and the setup hook, validate the grid.
  module subroutine SysInteraction_Prolate_Coulomb_StdImpl_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysInteraction
    use M_SysInteraction_Prolate

    call Say_Fabricate("sysInteraction.prolate.coulomb.stdImpl")

    !------------------------------------
    ! bind procedure pointers
    !------------------------------------

    SysInteraction_Prolate_FillInteractionPotentialChannel => FillInteractionPotentialChannel
    SysInteraction_Setup => Setup

    !------------------------------------
    ! validate grid requirements
    !------------------------------------

    if (.not. Json_GetExistence("grid.prolate.fedvr")) then
      error stop "grid.prolate.fedvr is required for sysInteraction.prolate.coulomb.stdImpl"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Precompute Legendre tables, tail coefficients, and LU factors.
  subroutine Setup
    use M_Utils_Say
    use M_Utils_Legendre
    use M_Utils_LapackLib
    use M_Grid_Prolate
    use M_SysInteraction_Prolate
    use M_SysInteraction_Prolate_Coulomb

    integer(I32) :: lmax, mmaxAbs, m, tau, iXi, jEta, nSolve, iStart, i1, i2
    real(R64) :: kappa, qDiag
    real(R64), allocatable :: pTmp(:), qTmp(:)

    call Say_Setup("sysInteraction.prolate.coulomb.stdImpl")

    lmax = SysInteraction_Prolate_Coulomb_lmax
    mmaxAbs = min(SysInteraction_Prolate_mmax, lmax)

    allocate (SysInteraction_Prolate_Coulomb_StdImpl_pBarEta(0:lmax, Grid_Prolate_nEta, 0:mmaxAbs))
    allocate (SysInteraction_Prolate_Coulomb_StdImpl_pXi(0:lmax, Grid_Prolate_nXi, 0:mmaxAbs))
    allocate (SysInteraction_Prolate_Coulomb_StdImpl_corrCoef(0:lmax, 0:mmaxAbs))
    allocate (SysInteraction_Prolate_Coulomb_StdImpl_solvers(0:lmax, 0:mmaxAbs))

    SysInteraction_Prolate_Coulomb_StdImpl_pBarEta = 0.0_R64
    SysInteraction_Prolate_Coulomb_StdImpl_pXi = 0.0_R64
    SysInteraction_Prolate_Coulomb_StdImpl_corrCoef = 0.0_R64

    allocate (pTmp(0:lmax), qTmp(0:lmax))

    do m = 0, mmaxAbs

      !------------------------------------
      ! Legendre tables
      !------------------------------------

      ! Normalized on-cut Pbar_tau^m at the eta nodes
      do jEta = 1, Grid_Prolate_nEta
        call Legendre_FillP(pTmp, lmax, m, Grid_Prolate_etaPoints(jEta))
        do tau = m, lmax
          SysInteraction_Prolate_Coulomb_StdImpl_pBarEta(tau, jEta, m) = Legendre_NormFactorP(tau, m) * pTmp(tau)
        end do
      end do

      ! Off-cut P_tau^m at the xi nodes (P = 0 at xi = 1 for m /= 0)
      do iXi = 1, Grid_Prolate_nXi
        call Legendre_FillP(pTmp, lmax, m, Grid_Prolate_xiPoints(iXi))
        SysInteraction_Prolate_Coulomb_StdImpl_pXi(m:lmax, iXi, m) = pTmp(m:lmax)
      end do

      ! Exterior-tail coefficients at the box boundary ximax
      call Legendre_FillP(pTmp, lmax, m, Grid_Prolate_ximax)
      call Legendre_FillQ(qTmp, lmax, m, Grid_Prolate_ximax)
      do tau = m, lmax
        kappa = (-1.0_R64)**(m + 1) * Legendre_FactorialRatio(tau, m)
        SysInteraction_Prolate_Coulomb_StdImpl_corrCoef(tau, m) = qTmp(tau) / (kappa * pTmp(tau))
      end do

      !------------------------------------
      ! LU factorizations of A = S_xi + diag(w_xi (tau(tau+1) + m^2/(xi^2-1)))
      !------------------------------------

      ! For m /= 0 the xi = 1 point (index 1) is excluded from the system
      iStart = 1
      if (m > 0) iStart = 2
      nSolve = Grid_Prolate_nXi - iStart + 1

      do tau = m, lmax
        associate (solver => SysInteraction_Prolate_Coulomb_StdImpl_solvers(tau, m))

          allocate (solver % lu(nSolve, nSolve))
          allocate (solver % ipiv(nSolve))

          do i2 = 1, nSolve
            do i1 = 1, nSolve
              solver % lu(i1, i2) = Grid_Prolate_xiKinMatrix(iStart - 1 + i1, iStart - 1 + i2)
            end do
          end do

          do i1 = 1, nSolve
            associate (xi => Grid_Prolate_xiPoints(iStart - 1 + i1))
              qDiag = tau * (tau + 1.0_R64)
              if (m > 0) qDiag = qDiag + m * m / (xi**2 - 1.0_R64)
              solver % lu(i1, i1) = solver % lu(i1, i1) + Grid_Prolate_xiWeights(iStart - 1 + i1) * qDiag
            end associate
          end do

          call LapackLib_FactorizeLU(solver % lu, solver % ipiv)

        end associate
      end do

    end do

    deallocate (pTmp, qTmp)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Solve one azimuthal channel of the interaction potential.
  !>
  !> @details
  !> For each degree tau = |m|..lmax:
  !>   1. Project the weighted source on Pbar_tau^m(eta):  g(xi)
  !>   2. Solve the weak-form system A v = (4 pi / a) g (Dirichlet at ximax)
  !>   3. Add the exterior tail: v += moment * corrCoef * P_tau^m(xi) with
  !>      moment = -(4 pi / a) sum_i P_tau^m(xi_i) g(xi_i)
  !>   4. Accumulate V_m(xi, eta) += v(xi) Pbar_tau^m(eta)
  subroutine FillInteractionPotentialChannel(potM, srcM, m, time, bt1_, bt2_)
    use M_Utils_Constants
    use M_Utils_UnusedVariables
    use M_Utils_LapackLib
    use M_Grid_Prolate
    use M_SysInteraction_Prolate_Coulomb

    complex(R64), intent(out), contiguous :: potM(:)
    complex(R64), intent(in), contiguous :: srcM(:)
    integer(I32), intent(in) :: m
    real(R64), intent(in) :: time
    integer(I32), intent(in), optional  :: bt1_
    integer(I32), intent(in), optional  :: bt2_

    integer(I32) :: mAbs, tau, iStart, nSolve, iXi, jEta, iGrid
    complex(R64) :: moment
    real(R64) :: sourceFactor
    complex(R64), allocatable :: g(:), v(:)
    real(R64), allocatable :: rhs(:, :)

    if (.false.) call UnusedVariables_Mark(time, bt1_, bt2_)

    potM = (0.0_R64, 0.0_R64)

    mAbs = abs(m)
    if (mAbs > SysInteraction_Prolate_Coulomb_lmax) return

    iStart = 1
    if (mAbs > 0) iStart = 2
    nSolve = Grid_Prolate_nXi - iStart + 1

    sourceFactor = 4.0_R64 * PI / Grid_Prolate_a * SysInteraction_Prolate_Coulomb_strength

    allocate (g(Grid_Prolate_nXi), v(Grid_Prolate_nXi))
    allocate (rhs(nSolve, 2))

    do tau = mAbs, SysInteraction_Prolate_Coulomb_lmax
      associate (pBarEta => SysInteraction_Prolate_Coulomb_StdImpl_pBarEta(tau, :, mAbs), &
                 pXi => SysInteraction_Prolate_Coulomb_StdImpl_pXi(tau, :, mAbs), &
                 solver => SysInteraction_Prolate_Coulomb_StdImpl_solvers(tau, mAbs))

        ! 1. eta projection of the weighted source
        g = (0.0_R64, 0.0_R64)
        iGrid = 0
        do jEta = 1, Grid_Prolate_nEta
          do iXi = 1, Grid_Prolate_nXi
            iGrid = iGrid + 1
            g(iXi) = g(iXi) + pBarEta(jEta) * srcM(iGrid)
          end do
        end do

        if (all(abs(g) < 1.0e-16_R64)) cycle

        ! 2. weak-form solve with homogeneous Dirichlet boundary
        rhs(:, 1) = sourceFactor * real(g(iStart:), kind=R64)
        rhs(:, 2) = sourceFactor * aimag(g(iStart:))

        call LapackLib_SolveFactorized(solver % lu, solver % ipiv, rhs, 'N')

        v = (0.0_R64, 0.0_R64)
        v(iStart:) = cmplx(rhs(:, 1), rhs(:, 2), kind=R64)

        ! 3. exterior tail correction (restores the Q_tau^m(xi) behavior at ximax)
        moment = -sourceFactor * sum(pXi(:) * g(:))
        v(:) = v(:) + moment * SysInteraction_Prolate_Coulomb_StdImpl_corrCoef(tau, mAbs) * pXi(:)

        ! 4. accumulate the channel potential
        iGrid = 0
        do jEta = 1, Grid_Prolate_nEta
          do iXi = 1, Grid_Prolate_nXi
            iGrid = iGrid + 1
            potM(iGrid) = potM(iGrid) + v(iXi) * pBarEta(jEta)
          end do
        end do

      end associate
    end do

    deallocate (g, v, rhs)

  end subroutine

end submodule
