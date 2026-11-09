! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file S_SysKinetic_Prolate_Laplacian_Fedvr.f90
!> @brief Implementation of the FEDVR kinetic operator for prolate channels.
submodule(M_SysKinetic_Prolate_Laplacian_Fedvr) S_SysKinetic_Prolate_Laplacian_Fedvr

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Bind FEDVR channel operator, validate grid requirements.
  module subroutine SysKinetic_Prolate_Laplacian_Fedvr_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysKinetic
    use M_SysKinetic_Prolate

    call Say_Fabricate("sysKinetic.prolate.laplacian.fedvr")

    !------------------------------------
    ! bind procedure pointers
    !------------------------------------

    SysKinetic_Prolate_MultiplyWithChannelKineticOp => MultiplyWithChannelKineticOp
    SysKinetic_SolveShifted => SolveShifted

    !------------------------------------
    ! validate grid requirements
    !------------------------------------

    if (.not. Json_GetExistence("grid.prolate.fedvr")) then
      error stop "grid.prolate.fedvr is required for sysKinetic.prolate.laplacian.fedvr"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Apply the kinetic operator to one azimuthal channel.
  !>
  !> @details
  !> Algorithm on the (nXi, nEta) channel matrix f:
  !>   1. Sturm-Liouville terms via the grid's symmetric DVR matrices:
  !>      [d/dxi (xi^2-1) d/dxi f]  = -(1/w_xi)  S_xi  f    (rows)
  !>      [d/deta (1-eta^2) d/deta f] = -(1/w_eta) f S_eta  (columns)
  !>   2. Centrifugal term m^2 (1/(xi^2-1) + 1/(1-eta^2)) f for m /= 0
  !>   3. Common prefactor -1/(2 M a^2 (xi^2 - eta^2))
  !>
  !> For m /= 0 the channel physically vanishes on the internuclear axis like
  !> (xi^2-1)^(|m|/2); the xi = 1 grid point (where the centrifugal term
  !> diverges) is therefore treated as a per-channel Dirichlet point: its
  !> input and output values are zeroed.
  subroutine MultiplyWithChannelKineticOp(dOrbM, orbM, m, time, bt_)
    use M_Utils_UnusedVariables
    use M_Grid
    use M_Grid_Prolate
    use M_SysKinetic_Prolate_Laplacian

    complex(R64), intent(out) :: dOrbM(:)
    complex(R64), intent(in)  :: orbM(:)
    integer(I32), intent(in)  :: m
    real(R64), intent(in)     :: time
    integer(I32), intent(in), optional :: bt_

    complex(R64), allocatable :: f(:, :), res(:, :)
    integer(I32) :: bt, iXi, jEta, parity
    real(R64) :: centrifugal

    if (.false.) call UnusedVariables_Mark(time)

    bt = 1
    if (present(bt_)) bt = bt_

    allocate (f(Grid_Prolate_nXi, Grid_Prolate_nEta))
    allocate (res(Grid_Prolate_nXi, Grid_Prolate_nEta))

    f = reshape(orbM, [Grid_Prolate_nXi, Grid_Prolate_nEta])

    ! Dirichlet at xi = 1 for m /= 0 (channel vanishes on the axis)
    if (m .ne. 0) f(1, :) = (0.0_R64, 0.0_R64)

    ! Sturm-Liouville derivative terms (positive semi-definite convention:
    ! -d/dxi (xi^2-1) d/dxi corresponds to +S_xi/w_xi)
    ! (the odd-|m| matrices use the axis-adapted basis, see M_Grid_Prolate)
    parity = mod(abs(m), 2)
    res = matmul(Grid_Prolate_xiKinMatrix(:, :, parity) + m * m * Grid_Prolate_xiCentrifugalMatrix(:, :, parity), f)
    do iXi = 1, Grid_Prolate_nXi
      res(iXi, :) = res(iXi, :) / Grid_Prolate_xiWeights(iXi)
    end do

    f = matmul(f, Grid_Prolate_etaKinMatrix(:, :, parity)) ! S_eta is symmetric
    do jEta = 1, Grid_Prolate_nEta
      res(:, jEta) = res(:, jEta) + f(:, jEta) / Grid_Prolate_etaWeights(jEta)
    end do

    ! Restore f for the centrifugal term
    f = reshape(orbM, [Grid_Prolate_nXi, Grid_Prolate_nEta])

    do jEta = 1, Grid_Prolate_nEta
      do iXi = 1, Grid_Prolate_nXi

        ! eta centrifugal term (the xi one is part of the matrix product above)
        if (m .ne. 0 .and. iXi > 1) then
          centrifugal = m * m / (1.0_R64 - Grid_Prolate_etaPoints(jEta)**2)
          res(iXi, jEta) = res(iXi, jEta) + centrifugal * f(iXi, jEta)
        end if

        ! Common prefactor: T = +1/(2 M a^2 (xi^2-eta^2)) [S terms + centrifugal]
        res(iXi, jEta) = res(iXi, jEta) / &
                         (2.0_R64 * SysKinetic_Prolate_Laplacian_bodyMass(bt) * Grid_Prolate_a**2 * &
                          (Grid_Prolate_xiPoints(iXi)**2 - Grid_Prolate_etaPoints(jEta)**2))
      end do
    end do

    if (m .ne. 0) res(1, :) = (0.0_R64, 0.0_R64)

    dOrbM(:) = reshape(res, [Grid_Prolate_nSpatial])

    deallocate (f, res)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Exact shifted solve dOrb = (T + shift)^-1 orb for a real shift.
  !>
  !> @details
  !> Per azimuthal channel, multiplying (T + sigma) F = R by the prefactor
  !> 2 M a^2 (xi^2 - eta^2) separates the operator (the shift term
  !> sigma 2 M a^2 (xi^2 - eta^2) splits into a xi and an eta part):
  !>
  !>   A_xi F + F B_eta = R .* 2 M a^2 (xi^2 - eta^2),
  !>   A_xi  = W_xi^-1 (S_xi + m^2 C_xi) + diag(2 M sigma a^2 xi^2),
  !>   B_eta = S_eta W_eta^-1 + diag(m^2/(1-eta^2) - 2 M sigma a^2 eta^2).
  !>
  !> With the symmetric forms W^-1/2 S W^-1/2 + diag = Q Lambda Q^T of both
  !> directions (fast diagonalization), X = Q_xi^T W_xi^1/2 F W_eta^1/2 Q_eta
  !> solves Lambda X + X M = (transformed right-hand side) elementwise. For
  !> m /= 0 the axis point xi = 1 is excluded (channel Dirichlet point, as in
  !> the operator). The 1D decompositions depend on the shift and are built
  !> per call (nXi and nEta are small).
  !>
  !> Real shifts only (the preconditioner of the eigensolvers); a complex
  !> shift (Crank-Nicolson) would need complex-symmetric decompositions.
  subroutine SolveShifted(dOrb, orb, shift, bt_)
    use M_Utils_LapackLib
    use M_Utils_BlasLib
    use M_Grid_Prolate
    use M_SysKinetic_Prolate_Laplacian

    complex(R64), intent(out), contiguous :: dOrb(:)
    complex(R64), intent(in), contiguous  :: orb(:)
    complex(R64), intent(in)              :: shift
    integer(I32), intent(in), optional    :: bt_

    integer(I32) :: bt, mChannel, mAbs, parity, iStart, nXiS, nEta, nXi, iXi, jEta, i1, i2, iPart, base
    real(R64) :: sigma, twoMa2
    real(R64), allocatable :: aXi(:, :), bEta(:, :), qXi(:, :), qEta(:, :), lambdaXi(:), muEta(:)
    real(R64), allocatable :: wXiSqrt(:), wEtaSqrt(:), rhs(:, :), tmp(:, :), x(:, :)

    bt = 1
    if (present(bt_)) bt = bt_

    if (abs(aimag(shift)) > 0.0_R64) &
      error stop "SysKinetic_Prolate_Laplacian_Fedvr: SolveShifted supports real shifts only"
    sigma = real(shift, kind=R64)

    nXi = Grid_Prolate_nXi
    nEta = Grid_Prolate_nEta
    twoMa2 = 2.0_R64 * SysKinetic_Prolate_Laplacian_bodyMass(bt) * Grid_Prolate_a**2

    allocate (wEtaSqrt(nEta), bEta(nEta, nEta))
    wEtaSqrt(:) = sqrt(Grid_Prolate_etaWeights(:))

    dOrb(:) = (0.0_R64, 0.0_R64)

    do mChannel = 0, 2 * Grid_Prolate_mmax
      mAbs = (mChannel + 1) / 2
      parity = mod(mAbs, 2)
      base = mChannel * Grid_Prolate_nSpatial

      iStart = 1
      if (mAbs > 0) iStart = 2
      nXiS = nXi - iStart + 1

      ! symmetric xi operator W^-1/2 S W^-1/2 + diag on the solved points
      allocate (wXiSqrt(nXiS), aXi(nXiS, nXiS))
      wXiSqrt(:) = sqrt(Grid_Prolate_xiWeights(iStart:nXi))
      do i2 = 1, nXiS
        do i1 = 1, nXiS
          aXi(i1, i2) = (Grid_Prolate_xiKinMatrix(iStart - 1 + i1, iStart - 1 + i2, parity) + &
                         mAbs**2 * Grid_Prolate_xiCentrifugalMatrix(iStart - 1 + i1, iStart - 1 + i2, parity)) / &
                        (wXiSqrt(i1) * wXiSqrt(i2))
        end do
      end do
      do i1 = 1, nXiS
        associate (xi => Grid_Prolate_xiPoints(iStart - 1 + i1))
          aXi(i1, i1) = aXi(i1, i1) + twoMa2 * sigma * xi**2
        end associate
      end do
      call LapackLib_DiagonalizeGeneric(lambdaXi, qXi, aXi, .true.)

      ! symmetric eta operator
      do i2 = 1, nEta
        do i1 = 1, nEta
          bEta(i1, i2) = Grid_Prolate_etaKinMatrix(i1, i2, parity) / (wEtaSqrt(i1) * wEtaSqrt(i2))
        end do
      end do
      do jEta = 1, nEta
        associate (eta => Grid_Prolate_etaPoints(jEta))
          bEta(jEta, jEta) = bEta(jEta, jEta) - twoMa2 * sigma * eta**2
          if (mAbs > 0) bEta(jEta, jEta) = bEta(jEta, jEta) + mAbs**2 / (1.0_R64 - eta**2)
        end associate
      end do
      call LapackLib_DiagonalizeGeneric(muEta, qEta, bEta, .true.)

      allocate (rhs(nXiS, nEta), tmp(nXiS, nEta), x(nXiS, nEta))

      do iPart = 1, 2
        ! right-hand side R .* 2 M a^2 (xi^2 - eta^2), scaled by W_xi^1/2 and W_eta^1/2
        do jEta = 1, nEta
          do i1 = 1, nXiS
            iXi = iStart - 1 + i1
            associate (r => orb(base + (jEta - 1) * nXi + iXi))
              if (iPart .eq. 1) then
                rhs(i1, jEta) = real(r, kind=R64)
              else
                rhs(i1, jEta) = aimag(r)
              end if
            end associate
            rhs(i1, jEta) = rhs(i1, jEta) * twoMa2 * (Grid_Prolate_xiPoints(iXi)**2 - Grid_Prolate_etaPoints(jEta)**2) * &
                            wXiSqrt(i1) * wEtaSqrt(jEta)
          end do
        end do

        ! X = Q_xi^T rhs Q_eta, divided by lambda_i + mu_j
        call BlasLib_CalcMatOnMat(tmp, qXi, rhs, 'T', 'N')
        call BlasLib_CalcMatOnMat(x, tmp, qEta, 'N', 'N')
        do jEta = 1, nEta
          x(:, jEta) = x(:, jEta) / (lambdaXi(:) + muEta(jEta))
        end do

        ! F = W_xi^-1/2 Q_xi X Q_eta^T W_eta^-1/2
        call BlasLib_CalcMatOnMat(tmp, qXi, x, 'N', 'N')
        call BlasLib_CalcMatOnMat(rhs, tmp, qEta, 'N', 'T')

        do jEta = 1, nEta
          do i1 = 1, nXiS
            iXi = iStart - 1 + i1
            associate (f => dOrb(base + (jEta - 1) * nXi + iXi))
              if (iPart .eq. 1) then
                f = cmplx(rhs(i1, jEta) / (wXiSqrt(i1) * wEtaSqrt(jEta)), 0.0_R64, kind=R64)
              else
                f = f + cmplx(0.0_R64, rhs(i1, jEta) / (wXiSqrt(i1) * wEtaSqrt(jEta)), kind=R64)
              end if
            end associate
          end do
        end do
      end do

      deallocate (wXiSqrt, aXi, qXi, lambdaXi, qEta, muEta, rhs, tmp, x)
    end do

  end subroutine

end submodule
