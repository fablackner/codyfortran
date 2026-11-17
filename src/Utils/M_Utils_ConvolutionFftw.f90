! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> FFT-based 1D convolution utilities built on top of FFTW3.
!>
!> Provides a small context type that plans forward/inverse FFTs and allocates
!> work buffers to apply linear convolutions efficiently in O(N log N), using
!> raw FFTW ordering for spectral data.
module M_Utils_ConvolutionFftw
  use M_Utils_Types
  use M_Utils_FftwLib

  implicit none

  private
  public :: T_ConvolutionFftw_Ctx
  public :: ConvolutionFftw_CreateCtx
  public :: ConvolutionFftw_DestroyCtx
  public :: ConvolutionFftw_Apply

  type :: T_ConvolutionFftw_Ctx
    integer(I32) :: nG = 0
    integer(I32) :: nPad = 0
    real(R64)    :: dx = 0.0_R64
    type(T_FftwLib_Ctx) :: fftCtx
    complex(R64), allocatable :: kernelK(:)
    complex(R64), allocatable :: srcPadded(:), srcK(:), convK(:), convPadded(:)
  end type

  abstract interface
    pure function I_KernelFun(distance) result(res)
      import :: R64
      real(R64) :: res
      real(R64), intent(in) :: distance
    end function
  end interface

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ConvolutionFftw_CreateCtx(ctx, nG, dx, KernelFun)
    type(T_ConvolutionFftw_Ctx), intent(inout) :: ctx
    integer(I32), intent(in) :: nG
    real(R64), intent(in) :: dx
    procedure(I_KernelFun) :: KernelFun
    integer(I32) :: i
    real(R64) :: distance, kernelVal
    complex(R64), allocatable :: kernelPadded(:)

    ctx % nG = nG
    ctx % dx = dx
    ctx % nPad = 2 * nG
    call FftwLib_CreateCtx(ctx % fftCtx, ctx % nPad)

    allocate (kernelPadded(ctx % nPad))
    kernelPadded = (0.0_R64, 0.0_R64)

    do i = 0, nG - 1
      distance = real(i, R64) * dx
      kernelVal = KernelFun(distance)
      if (i .eq. 0) then
        kernelPadded(1) = cmplx(kernelVal, 0.0_R64, kind=R64)
      else
        kernelPadded(i + 1) = cmplx(kernelVal, 0.0_R64, kind=R64)
        kernelPadded(ctx % nPad - i + 1) = cmplx(kernelVal, 0.0_R64, kind=R64)
      end if
    end do

    allocate (ctx % kernelK(ctx % nPad))
    call FftwLib_DoFt(ctx % kernelK, kernelPadded, dx, ctx % fftCtx)
    deallocate (kernelPadded)

    allocate (ctx % srcPadded(ctx % nPad), ctx % srcK(ctx % nPad), ctx % convK(ctx % nPad), ctx % convPadded(ctx % nPad))
    ctx % srcPadded = (0.0_R64, 0.0_R64)
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ConvolutionFftw_DestroyCtx(ctx)
    type(T_ConvolutionFftw_Ctx), intent(inout) :: ctx
    if (allocated(ctx % kernelK)) deallocate (ctx % kernelK)
    if (allocated(ctx % srcPadded)) deallocate (ctx % srcPadded)
    if (allocated(ctx % srcK)) deallocate (ctx % srcK)
    if (allocated(ctx % convK)) deallocate (ctx % convK)
    if (allocated(ctx % convPadded)) deallocate (ctx % convPadded)
    call FftwLib_DestroyCtx(ctx % fftCtx)
    ctx % nG = 0; ctx % nPad = 0; ctx % dx = 0.0_R64
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ConvolutionFftw_Apply(dest, src, ctx)
    complex(R64), intent(out) :: dest(:)
    complex(R64), intent(in)  :: src(:)
    type(T_ConvolutionFftw_Ctx), intent(inout) :: ctx

    ctx % srcPadded = (0.0_R64, 0.0_R64)
    ctx % srcPadded(1:ctx % nG) = src

    call FftwLib_DoFt(ctx % srcK, ctx % srcPadded, ctx % dx, ctx % fftCtx)
    ctx % convK = (ctx % srcK * ctx % kernelK) / ctx % dx
    call FftwLib_DoInverseFt(ctx % convPadded, ctx % convK, ctx % dx, ctx % fftCtx)
    dest = ctx % convPadded(1:ctx % nG)
  end subroutine

end module
