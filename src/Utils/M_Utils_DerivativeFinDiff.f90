! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> High-order finite-difference derivatives on uniform 1D grids.
!>
!> Provides precomputed stencils for 1st and 2nd derivatives with selectable
!> accuracy orders near boundaries and in the bulk. Useful as a non-FFT
!> alternative and for validation of spectral derivatives.
module M_Utils_DerivativeFinDiff

  use M_Utils_Types

  implicit none

  type :: T_DerivativeFinDiff_Ctx
    integer(I32) :: nG = 0
    real(R64)    :: dx = 0.0_R64
    real(R64)    :: dxSq = 0.0_R64

    ! First derivative coefficients (7-point one–sided, 8th-order centered interior)
    real(R64) :: cd1_1, cd1_2, cd1_3, cd1_4               ! centered antisymmetric
    real(R64) :: fd1_0, fd1_1, fd1_2, fd1_3, fd1_4, fd1_5, fd1_6  ! forward (left), backward use sign

    ! Second derivative coefficients (8th-order centered interior)
    real(R64) :: cd2_0, cd2_1, cd2_2, cd2_3, cd2_4
    real(R64) :: fd2_0, fd2_1, fd2_2, fd2_3, fd2_4, fd2_5, fd2_6, fd2_7
  end type

  public :: T_DerivativeFinDiff_Ctx

  public :: DerivativeFinDiff_CreateCtx, DerivativeFinDiff_DestroyCtx

  public :: DerivativeFinDiff_Do1stDerivative

  public :: DerivativeFinDiff_Do2ndDerivative
  public :: DerivativeFinDiff_Do2ndDerivativePeriodic
  public :: DerivativeFinDiff_Do2ndDerivativeZeroEnds

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine DerivativeFinDiff_CreateCtx(ctx, nG, dx)
    type(T_DerivativeFinDiff_Ctx), intent(out) :: ctx
    integer(I32), intent(in) :: nG
    real(R64), intent(in) :: dx

    ctx % nG = nG
    ctx % dx = dx
    ctx % dxSq = dx * dx

    ! First derivative coefficients
    ctx % cd1_1 = 4.0_R64 / 5.0_R64
    ctx % cd1_2 = -1.0_R64 / 5.0_R64
    ctx % cd1_3 = 4.0_R64 / 105.0_R64
    ctx % cd1_4 = -1.0_R64 / 280.0_R64

    ctx % fd1_0 = -49.0_R64 / 20.0_R64
    ctx % fd1_1 = 6.0_R64
    ctx % fd1_2 = -15.0_R64 / 2.0_R64
    ctx % fd1_3 = 20.0_R64 / 3.0_R64
    ctx % fd1_4 = -15.0_R64 / 4.0_R64
    ctx % fd1_5 = 6.0_R64 / 5.0_R64
    ctx % fd1_6 = -1.0_R64 / 6.0_R64

    ! Second derivative coefficients
    ctx % cd2_0 = -205.0_R64 / 72.0_R64
    ctx % cd2_1 = 8.0_R64 / 5.0_R64
    ctx % cd2_2 = -1.0_R64 / 5.0_R64
    ctx % cd2_3 = 8.0_R64 / 315.0_R64
    ctx % cd2_4 = -1.0_R64 / 560.0_R64

    ctx % fd2_0 = 469.0_R64 / 90.0_R64
    ctx % fd2_1 = -223.0_R64 / 10.0_R64
    ctx % fd2_2 = 879.0_R64 / 20.0_R64
    ctx % fd2_3 = -949.0_R64 / 18.0_R64
    ctx % fd2_4 = 41.0_R64
    ctx % fd2_5 = -201.0_R64 / 10.0_R64
    ctx % fd2_6 = 1019.0_R64 / 180.0_R64
    ctx % fd2_7 = -7.0_R64 / 10.0_R64
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine DerivativeFinDiff_DestroyCtx(ctx)
    use M_Utils_UnusedVariables

    type(T_DerivativeFinDiff_Ctx), intent(inout) :: ctx

    ! no deallocate needed
    if (.false.) call UnusedVariables_Mark(ctx)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine DerivativeFinDiff_Do1stDerivative(dfDx, f, ctx)
    complex(R64), intent(out) :: dfDx(:)
    complex(R64), intent(in)  :: f(:)
    type(T_DerivativeFinDiff_Ctx), intent(in) :: ctx

    integer(I32) :: n, i

    n = size(f)
    if (n .ne. ctx % nG) stop "DerivativeFinDiff_Do1stDerivative: size mismatch."
    if (n < 10) stop "DerivativeFinDiff_Do1stDerivative: need at least 10 points."

    dfDx = (0.0_R64, 0.0_R64)

    dfDx = (ctx % cd1_4 * (eoshift(f, 4) - eoshift(f, -4)) + &
            ctx % cd1_3 * (eoshift(f, 3) - eoshift(f, -3)) + &
            ctx % cd1_2 * (eoshift(f, 2) - eoshift(f, -2)) + &
            ctx % cd1_1 * (eoshift(f, 1) - eoshift(f, -1))) / ctx % dx

    do i = 1, 4
      dfDx(i) = (ctx % fd1_0 * f(i) + ctx % fd1_1 * f(i + 1) + ctx % fd1_2 * f(i + 2) + &
                 ctx % fd1_3 * f(i + 3) + ctx % fd1_4 * f(i + 4) + ctx % fd1_5 * f(i + 5) + &
                 ctx % fd1_6 * f(i + 6)) / ctx % dx
    end do

    do i = n - 3, n
      dfDx(i) = -(ctx % fd1_0 * f(i) + ctx % fd1_1 * f(i - 1) + ctx % fd1_2 * f(i - 2) + &
                  ctx % fd1_3 * f(i - 3) + ctx % fd1_4 * f(i - 4) + ctx % fd1_5 * f(i - 5) + &
                  ctx % fd1_6 * f(i - 6)) / ctx % dx
    end do
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine DerivativeFinDiff_Do2ndDerivative(d2fDx2, f, ctx)
    complex(R64), intent(out) :: d2fDx2(:)
    complex(R64), intent(in) :: f(:)
    type(T_DerivativeFinDiff_Ctx), intent(in) :: ctx

    integer(I32) :: n, i

    n = size(f)
    if (n .ne. ctx % nG) stop "DerivativeFinDiff_Do2ndDerivative: size mismatch."
    if (n < 11) stop "DerivativeFinDiff_Do2ndDerivative: need at least 11 points."

    d2fDx2 = (0.0_R64, 0.0_R64)

    d2fDx2 = (ctx % cd2_0 * f + &
              ctx % cd2_1 * (eoshift(f, 1) + eoshift(f, -1)) + &
              ctx % cd2_2 * (eoshift(f, 2) + eoshift(f, -2)) + &
              ctx % cd2_3 * (eoshift(f, 3) + eoshift(f, -3)) + &
              ctx % cd2_4 * (eoshift(f, 4) + eoshift(f, -4))) / ctx % dxSq

    do i = 1, 4
      d2fDx2(i) = (ctx % fd2_0 * f(i) + ctx % fd2_1 * f(i + 1) + ctx % fd2_2 * f(i + 2) + &
                   ctx % fd2_3 * f(i + 3) + ctx % fd2_4 * f(i + 4) + ctx % fd2_5 * f(i + 5) + &
                   ctx % fd2_6 * f(i + 6) + ctx % fd2_7 * f(i + 7)) / ctx % dxSq
    end do

    do i = n - 3, n
      d2fDx2(i) = (ctx % fd2_0 * f(i) + ctx % fd2_1 * f(i - 1) + ctx % fd2_2 * f(i - 2) + &
                   ctx % fd2_3 * f(i - 3) + ctx % fd2_4 * f(i - 4) + ctx % fd2_5 * f(i - 5) + &
                   ctx % fd2_6 * f(i - 6) + ctx % fd2_7 * f(i - 7)) / ctx % dxSq
    end do
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine DerivativeFinDiff_Do2ndDerivativePeriodic(d2fDx2, f, ctx)
    complex(R64), intent(out) :: d2fDx2(:)
    complex(R64), intent(in) :: f(:)
    type(T_DerivativeFinDiff_Ctx), intent(in) :: ctx

    integer(I32) :: n

    n = size(f)
    if (n .ne. ctx % nG) stop "DerivativeFinDiff_Do2ndDerivativePeriodic: size mismatch."
    if (n < 9) stop "DerivativeFinDiff_Do2ndDerivativePeriodic: need at least 9 points."

    d2fDx2 = (ctx % cd2_0 * f + &
              ctx % cd2_1 * (cshift(f, 1) + cshift(f, -1)) + &
              ctx % cd2_2 * (cshift(f, 2) + cshift(f, -2)) + &
              ctx % cd2_3 * (cshift(f, 3) + cshift(f, -3)) + &
              ctx % cd2_4 * (cshift(f, 4) + cshift(f, -4))) / ctx % dxSq
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine DerivativeFinDiff_Do2ndDerivativeZeroEnds(d2fDx2, f, ctx)
    complex(R64), intent(out) :: d2fDx2(:)
    complex(R64), intent(in) :: f(:)
    type(T_DerivativeFinDiff_Ctx), intent(in) :: ctx

    integer(I32) :: n

    n = size(f)
    if (n .ne. ctx % nG) stop "DerivativeFinDiff_Do2ndDerivativeZeroEnds: size mismatch."

    d2fDx2 = (0.0_R64, 0.0_R64)

    d2fDx2 = ctx % cd2_0 * f + &
             ctx % cd2_1 * (eoshift(f, 1, boundary=(0.0_R64, 0.0_R64)) + eoshift(f, -1, boundary=(0.0_R64, 0.0_R64))) + &
             ctx % cd2_2 * (eoshift(f, 2, boundary=(0.0_R64, 0.0_R64)) + eoshift(f, -2, boundary=(0.0_R64, 0.0_R64))) + &
             ctx % cd2_3 * (eoshift(f, 3, boundary=(0.0_R64, 0.0_R64)) + eoshift(f, -3, boundary=(0.0_R64, 0.0_R64))) + &
             ctx % cd2_4 * (eoshift(f, 4, boundary=(0.0_R64, 0.0_R64)) + eoshift(f, -4, boundary=(0.0_R64, 0.0_R64)))

    d2fDx2 = d2fDx2 / ctx % dxSq
  end subroutine

end module
