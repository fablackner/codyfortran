! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Direct (real-space) convolution/integral utilities.
!>
!> Implements a kernel-application context to evaluate discrete convolutions
!> or integral operators K*f via precomputed kernel samples on a uniform grid
!> without FFTs (useful for short/support-limited kernels or validation).
module M_Utils_ConvolutionIntegral
  use M_Utils_Types
  implicit none

  private
  public :: T_ConvolutionIntegral_Ctx
  public :: ConvolutionIntegral_CreateCtx
  public :: ConvolutionIntegral_DestroyCtx
  public :: ConvolutionIntegral_Apply

  type :: T_ConvolutionIntegral_Ctx
    integer(I32) :: nG = 0
    real(R64)    :: dx = 0.0_R64
    real(R64)    :: srcThreshold = 0.0_R64
    real(R64), allocatable :: kernelVals(:)  ! kernelVals(k) = K(k*dx), k = 0..nG-1
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
  subroutine ConvolutionIntegral_CreateCtx(ctx, nG, dx, KernelFun, srcThreshold_)

    type(T_ConvolutionIntegral_Ctx), intent(inout) :: ctx
    integer(I32), intent(in) :: nG
    real(R64), intent(in) :: dx
    procedure(I_KernelFun) :: KernelFun
    real(R64), intent(in), optional :: srcThreshold_

    integer(I32) :: k

    ctx % nG = nG
    ctx % dx = dx
    if (present(srcThreshold_)) ctx % srcThreshold = srcThreshold_

    allocate (ctx % kernelVals(0:nG - 1))
    do k = 0, nG - 1
      ctx % kernelVals(k) = KernelFun(real(k, R64) * dx)
    end do
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ConvolutionIntegral_DestroyCtx(ctx)
    type(T_ConvolutionIntegral_Ctx), intent(inout) :: ctx

    if (allocated(ctx % kernelVals)) deallocate (ctx % kernelVals)

    ctx % nG = 0
    ctx % dx = 0.0_R64
    ctx % srcThreshold = 0.0_R64

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ConvolutionIntegral_Apply(dest, src, ctx)
    complex(R64), intent(out), contiguous :: dest(:)
    complex(R64), intent(in), contiguous :: src(:)
    type(T_ConvolutionIntegral_Ctx), intent(in) :: ctx

    integer(I32) :: i, j, nG, k
    real(R64) :: thr

    nG = ctx % nG
    thr = ctx % srcThreshold
    dest = (0.0_R64, 0.0_R64)

    do j = 1, nG
      if (abs(src(j)) > thr) then
        do i = 1, nG
          k = abs(i - j)
          dest(i) = dest(i) + ctx % kernelVals(k) * src(j)
        end do
      end if
    end do
  end subroutine

end module
