! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> FFT-based spectral derivatives on a uniform 1D grid.
!>
!> This module computes first and second spatial derivatives using FFTW:
!> - Create a context with grid size N and spacing dx.
!> - Forward FFT to k-space, multiply by (ik) or (-k^2), and inverse FFT.
!> The k-grid follows raw FFTW ordering as provided by `M_Utils_FftwLib`.
module M_Utils_DerivativeFftw

  use M_Utils_Types
  use M_Utils_FftwLib

  implicit none

  type :: T_DerivativeFftw_Ctx
    integer(I32) :: N = 0
    real(R64)    :: dx = 0.0_R64
    type(T_FftwLib_Ctx) :: fftCtx
    real(R64), allocatable :: k(:)
    complex(R64), allocatable :: work(:)
  end type

  public :: T_DerivativeFftw_Ctx

  public :: DerivativeFftw_CreateCtx, DerivativeFftw_DestroyCtx
  public :: DerivativeFftw_Do1stDerivative, DerivativeFftw_Do2ndDerivative

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Initialize FFT plans and allocate work buffers for spectral derivatives.
  subroutine DerivativeFftw_CreateCtx(ctx, N, dx)
    type(T_DerivativeFftw_Ctx), intent(out) :: ctx
    integer(I32), intent(in) :: N
    real(R64), intent(in) :: dx

    ctx % N = N
    ctx % dx = dx

    call FftwLib_CreateCtx(ctx % fftCtx, N)

    allocate (ctx % k(N), ctx % work(N))

    ! build k
    ! Removed dummy FFTW call; build k directly below:
    ! Simple rebuild (duplicate of internal function):
    ! NOTE: direct formula avoids exposing private routine.

    block
      integer(I32):: j, m
      real(R64):: dk

      dk = 2.0_R64 * 3.1415926535897932384626433832795_R64 / (N * dx)

      do j = 1, N
        m = j - 1
        if (m <= N / 2) then
          ctx % k(j) = m * dk
        else
          ctx % k(j) = (m - N) * dk
        end if
      end do
    end block

  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Destroy FFT plans and release buffers associated with the derivative context.
  subroutine DerivativeFftw_DestroyCtx(ctx)
    type(T_DerivativeFftw_Ctx), intent(inout) :: ctx

    if (allocated(ctx % k)) deallocate (ctx % k)
    if (allocated(ctx % work)) deallocate (ctx % work)

    call FftwLib_DestroyCtx(ctx % fftCtx)

    ctx % N = 0
    ctx % dx = 0.0_R64
  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Compute the first derivative df/dx via spectral multiplication by (i k).
  subroutine DerivativeFftw_Do1stDerivative(dfDx, f, ctx)
    complex(R64), intent(out) :: dfDx(:)
    complex(R64), intent(in)  :: f(:)
    type(T_DerivativeFftw_Ctx), intent(in) :: ctx

    complex(R64), allocatable :: ft(:), fTmp(:)
    integer(I32):: i, N

    N = ctx % N

    if (size(f) .ne. N .or. size(dfDx) .ne. N) stop 'DerivativeFftw_Do1stDerivative: size mismatch'

    allocate (ft(N), fTmp(N))

    call FftwLib_DoFt(ft, f, ctx % dx, ctx % fftCtx)

    do i = 1, N
      ft(i) = (0.0_R64, 1.0_R64) * ctx % k(i) * ft(i)
    end do

    call FftwLib_DoInverseFt(fTmp, ft, ctx % dx, ctx % fftCtx)

    dfDx = fTmp

    deallocate (ft, fTmp)
  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Compute the second derivative d²f/dx² via spectral multiplication by (−k²).
  subroutine DerivativeFftw_Do2ndDerivative(d2fDx2, f, ctx)
    complex(R64), intent(out) :: d2fDx2(:)
    complex(R64), intent(in)  :: f(:)
    type(T_DerivativeFftw_Ctx), intent(in) :: ctx

    complex(R64), allocatable :: ft(:), fTmp(:)
    integer(I32):: i, N

    N = ctx % N

    if (size(f) .ne. N .or. size(d2fDx2) .ne. N) stop 'DerivativeFftw_Do2ndDerivative: size mismatch'

    allocate (ft(N), fTmp(N))

    call FftwLib_DoFt(ft, f, ctx % dx, ctx % fftCtx)

    do i = 1, N
      ft(i) = -(ctx % k(i)**2) * ft(i)
    end do

    call FftwLib_DoInverseFt(fTmp, ft, ctx % dx, ctx % fftCtx)

    d2fDx2 = fTmp

    deallocate (ft, fTmp)
  end subroutine

end module
