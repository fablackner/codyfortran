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
  !> @brief Precompute Legendre tables, tail coefficients, and the generalized
  !> eigen-decomposition K_m U = W U Lambda of every channel |m|.
  subroutine Setup
    use M_Utils_Say
    use M_Utils_Legendre
    use M_Utils_LapackLib
    use M_Grid_Prolate
    use M_SysInteraction_Prolate
    use M_SysInteraction_Prolate_Coulomb

    integer(I32) :: lmax, mmaxAbs, m, tau, iTau, nTau, iXi, jEta, nSolve, iStart, i1, i2
    real(R64) :: kappa
    real(R64), allocatable :: pTmp(:), qTmp(:), kMat(:, :), wInvSqrt(:), evecs(:, :)

    call Say_Setup("sysInteraction.prolate.coulomb.stdImpl")

    lmax = SysInteraction_Prolate_Coulomb_lmax
    mmaxAbs = min(SysInteraction_Prolate_mmax, lmax)

    if (allocated(SysInteraction_Prolate_Coulomb_StdImpl_channels)) deallocate (SysInteraction_Prolate_Coulomb_StdImpl_channels)
    allocate (SysInteraction_Prolate_Coulomb_StdImpl_channels(0:mmaxAbs))

    allocate (pTmp(0:lmax), qTmp(0:lmax))

    do m = 0, mmaxAbs
      associate (channel => SysInteraction_Prolate_Coulomb_StdImpl_channels(m))

        nTau = lmax - m + 1

        !------------------------------------
        ! Legendre tables and tail coefficients
        !------------------------------------

        allocate (channel % tauFactor(nTau), channel % corrCoef(nTau))
        allocate (channel % pBarEta(Grid_Prolate_nEta, nTau), channel % pXi(Grid_Prolate_nXi, nTau))

        do iTau = 1, nTau
          tau = m + iTau - 1
          channel % tauFactor(iTau) = tau * (tau + 1.0_R64)
        end do

        ! Normalized on-cut Pbar_tau^m at the eta nodes
        do jEta = 1, Grid_Prolate_nEta
          call Legendre_FillP(pTmp, lmax, m, Grid_Prolate_etaPoints(jEta))
          do iTau = 1, nTau
            tau = m + iTau - 1
            channel % pBarEta(jEta, iTau) = Legendre_NormFactorP(tau, m) * pTmp(tau)
          end do
        end do

        ! Off-cut P_tau^m at the xi nodes (P = 0 at xi = 1 for m /= 0)
        do iXi = 1, Grid_Prolate_nXi
          call Legendre_FillP(pTmp, lmax, m, Grid_Prolate_xiPoints(iXi))
          channel % pXi(iXi, :) = pTmp(m:lmax)
        end do

        ! Exterior-tail coefficients at the box boundary ximax
        call Legendre_FillP(pTmp, lmax, m, Grid_Prolate_ximax)
        call Legendre_FillQ(qTmp, lmax, m, Grid_Prolate_ximax)
        do iTau = 1, nTau
          tau = m + iTau - 1
          kappa = (-1.0_R64)**(m + 1) * Legendre_FactorialRatio(tau, m)
          channel % corrCoef(iTau) = qTmp(tau) / (kappa * pTmp(tau))
        end do

        !------------------------------------
        ! K_m = S_xi + m^2 C_xi on the solved xi points (matrices of the parity
        ! of m, axis-adapted for odd m); the
        ! symmetric form W^-1/2 K_m W^-1/2 = Q Lambda Q^T gives U = W^-1/2 Q
        !------------------------------------

        ! For m /= 0 the xi = 1 point (index 1) is excluded from the system
        iStart = 1
        if (m > 0) iStart = 2
        nSolve = Grid_Prolate_nXi - iStart + 1

        allocate (kMat(nSolve, nSolve), wInvSqrt(nSolve))

        do i1 = 1, nSolve
          wInvSqrt(i1) = 1.0_R64 / sqrt(Grid_Prolate_xiWeights(iStart - 1 + i1))
        end do

        do i2 = 1, nSolve
          do i1 = 1, nSolve
            kMat(i1, i2) = Grid_Prolate_xiKinMatrix(iStart - 1 + i1, iStart - 1 + i2, mod(m, 2)) + &
                           m * m * Grid_Prolate_xiCentrifugalMatrix(iStart - 1 + i1, iStart - 1 + i2, mod(m, 2))
          end do
        end do

        do i2 = 1, nSolve
          do i1 = 1, nSolve
            kMat(i1, i2) = wInvSqrt(i1) * kMat(i1, i2) * wInvSqrt(i2)
          end do
        end do

        call LapackLib_DiagonalizeGeneric(channel % lambda, evecs, kMat, .true.)

        allocate (channel % u(nSolve, nSolve))
        do i2 = 1, nSolve
          channel % u(:, i2) = wInvSqrt(:) * evecs(:, i2)
        end do

        deallocate (kMat, wInvSqrt, evecs)

      end associate
    end do

    deallocate (pTmp, qTmp)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Solve one azimuthal channel of the interaction potential.
  !>
  !> @details
  !> For all degrees tau = |m|..lmax at once (real and imaginary parts as
  !> separate columns):
  !>   1. Project the weighted source on Pbar_tau^m(eta):  G = srcM Pbar
  !>   2. Solve A_tau v_tau = (4 pi / a) g_tau with homogeneous Dirichlet
  !>      boundary: v_tau = (4 pi / a) U (Lambda + tau(tau+1))^-1 U^T g_tau
  !>   3. Add the exterior tail: v_tau += moment * corrCoef * P_tau^m(xi) with
  !>      moment = -(4 pi / a) sum_i P_tau^m(xi_i) g_tau(xi_i)
  !>   4. Synthesize V_m(xi, eta) = sum_tau v_tau(xi) Pbar_tau^m(eta) = V Pbar^T
  subroutine FillInteractionPotentialChannel(potM, srcM, m, time, bt1_, bt2_)
    use M_Utils_Constants
    use M_Utils_UnusedVariables
    use M_Utils_BlasLib
    use M_Grid_Prolate
    use M_SysInteraction_Prolate_Coulomb

    complex(R64), intent(out), contiguous :: potM(:)
    complex(R64), intent(in), contiguous :: srcM(:)
    integer(I32), intent(in) :: m
    real(R64), intent(in) :: time
    integer(I32), intent(in), optional  :: bt1_
    integer(I32), intent(in), optional  :: bt2_

    integer(I32) :: mAbs, iStart, nSolve, nTau, nXi, nEta, iTau, iPart
    real(R64) :: sourceFactor, moment
    real(R64), allocatable :: srcPart(:, :), g(:, :), h(:, :), v(:, :), potPart(:, :)

    if (.false.) call UnusedVariables_Mark(time, bt1_, bt2_)

    potM = (0.0_R64, 0.0_R64)

    mAbs = abs(m)
    if (mAbs > SysInteraction_Prolate_Coulomb_lmax) return

    associate (channel => SysInteraction_Prolate_Coulomb_StdImpl_channels(mAbs))

      nXi = Grid_Prolate_nXi
      nEta = Grid_Prolate_nEta
      nTau = size(channel % tauFactor)

      iStart = 1
      if (mAbs > 0) iStart = 2
      nSolve = nXi - iStart + 1

      sourceFactor = 4.0_R64 * PI / Grid_Prolate_a * SysInteraction_Prolate_Coulomb_strength

      ! columns 1:nTau real part, nTau+1:2 nTau imaginary part
      allocate (srcPart(nXi, nEta), g(nXi, 2 * nTau), h(nSolve, 2 * nTau), v(nXi, 2 * nTau), potPart(nXi, nEta))

      ! 1. eta projection of the weighted source
      srcPart = reshape(real(srcM, kind=R64), [nXi, nEta])
      call BlasLib_CalcMatOnMat(g(:, 1:nTau), srcPart, channel % pBarEta, 'N', 'N')
      srcPart = reshape(aimag(srcM), [nXi, nEta])
      call BlasLib_CalcMatOnMat(g(:, nTau + 1:2 * nTau), srcPart, channel % pBarEta, 'N', 'N')

      ! 2. Dirichlet solve of every degree in the eigenbasis of K_m
      call BlasLib_CalcMatOnMat(h, channel % u, g(iStart:nXi, :), 'T', 'N')
      do iPart = 0, 1
        do iTau = 1, nTau
          h(:, iPart * nTau + iTau) = h(:, iPart * nTau + iTau) / (channel % lambda(:) + channel % tauFactor(iTau))
        end do
      end do
      v = 0.0_R64
      call BlasLib_CalcMatOnMat(v(iStart:nXi, :), channel % u, h, 'N', 'N')
      v = sourceFactor * v

      ! 3. exterior tail correction (restores the Q_tau^m(xi) behavior at ximax)
      do iPart = 0, 1
        do iTau = 1, nTau
          moment = -sourceFactor * sum(channel % pXi(:, iTau) * g(:, iPart * nTau + iTau))
          v(:, iPart * nTau + iTau) = v(:, iPart * nTau + iTau) + &
                                      moment * channel % corrCoef(iTau) * channel % pXi(:, iTau)
        end do
      end do

      ! 4. synthesize the channel potential
      call BlasLib_CalcMatOnMat(potPart, v(:, 1:nTau), channel % pBarEta, 'N', 'T')
      potM(:) = reshape(potPart, [nXi * nEta])
      call BlasLib_CalcMatOnMat(potPart, v(:, nTau + 1:2 * nTau), channel % pBarEta, 'N', 'T')
      potM(:) = cmplx(real(potM, kind=R64), reshape(potPart, [nXi * nEta]), kind=R64)

      deallocate (srcPart, g, h, v, potPart)

    end associate

  end subroutine

end submodule
