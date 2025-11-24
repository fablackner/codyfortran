! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Thin wrappers around FFTW3 for 1D transforms and wavenumber handling.
!>
!> Exposes plans and helpers that adhere strictly to raw FFTW ordering for
!> spectral arrays and the transform convention documented below. Use this when
!> full control over normalization and k-ordering is required.
!
! This version outputs k_(:) ALWAYS in raw FFTW order:
!   Fortran index j -> integer mode m = j-1
!   dk = 2π / (N * dx)
!   if m <= N/2 then k_ =  m * dk
!   if m  > N/2 then k_ = (m - N) * dk   (negative frequencies wrapped)
! Even N includes a single Nyquist mode at m = N/2 (k_ = +N/2*dk ≡ -N/2*dk).
!
! Transform convention implemented (discrete Riemann sum approximation):
!   F(k_j) = Σ_{n=0}^{N-1} f(x_n) * exp(-i k_j x_n) * dx / sqrt(2π)
!   f(x_n) = Σ_{j=0}^{N-1} F(k_j) * exp(+i k_j x_n) * dk / sqrt(2π)
!
! where dx = x(2)-x(1) and dk = 2π / (N * dx).
!
! NOTE:
! - If you need a centered spectrum for visualization, do it externally.
!
module M_Utils_FftwLib
  use M_Utils_Types
  use, intrinsic :: iso_c_binding
  implicit none

  include 'fftw3.f03'

  type :: T_FftwLib_Ctx
    integer(I32) :: N = 0
    type(c_ptr) :: plan_forward = c_null_ptr
    type(c_ptr) :: plan_backward = c_null_ptr
  end type

  public :: T_FftwLib_Ctx
  public :: FftwLib_CreateCtx
  public :: FftwLib_DestroyCtx

  interface
    subroutine fftw_execute_dft_mod(p, in, out) bind(C, name='fftw_execute_dft')
      import
      type(C_PTR), value :: p
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(in) :: in
      complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
    end subroutine

  end interface

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pure subroutine BuildRawWaveNumbers(kArray, dx, N)
    use M_Utils_Constants

    ! Build raw FFTW-order wavenumber array for a uniform grid.
    real(R64), intent(out) :: kArray(:)
    real(R64), intent(in)  :: dx
    integer(I32), intent(in) :: N

    integer(I32) :: j, m
    real(R64) :: dk

    if (N <= 0) return
    dk = 2.0_R64 * PI / (N * dx)

    do j = 1, N
      m = j - 1
      if (m <= N / 2) then
        kArray(j) = m * dk
      else
        kArray(j) = (m - N) * dk
      end if
    end do
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FftwLib_CreateCtx(ctx, N)
    type(T_FftwLib_Ctx), intent(out) :: ctx
    integer(I32), intent(in) :: N

    complex(R64), allocatable :: work_in(:), work_out(:)

    allocate (work_in(N), work_out(N))

    ctx % N = N
    ctx % plan_forward = fftw_plan_dft_1d(N, work_in, work_out, FFTW_FORWARD, FFTW_ESTIMATE)
    ctx % plan_backward = fftw_plan_dft_1d(N, work_in, work_out, FFTW_BACKWARD, FFTW_ESTIMATE)

    deallocate (work_in, work_out)
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FftwLib_DestroyCtx(ctx)
    type(T_FftwLib_Ctx), intent(inout) :: ctx

    call fftw_destroy_plan(ctx % plan_forward)
    call fftw_destroy_plan(ctx % plan_backward)

    ctx % plan_forward = c_null_ptr
    ctx % plan_backward = c_null_ptr

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FftwLib_DoFt(outputSpectrum, inputFunction, dx, ctx, k_)
    use M_Utils_Constants
    complex(R64), intent(out)              :: outputSpectrum(:)
    complex(R64), intent(in)               :: inputFunction(:)
    real(R64), intent(in)                  :: dx
    type(T_FftwLib_Ctx), intent(in)       :: ctx
    real(R64), intent(out), allocatable, optional :: k_(:)

    ! Execute directly on the provided arrays (no copy into ctx%work_in/out)
    call fftw_execute_dft_mod(ctx % plan_forward, inputFunction, outputSpectrum)
    outputSpectrum = outputSpectrum * dx

    if (present(k_)) call BuildRawWaveNumbers(k_, dx, ctx % N)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FftwLib_DoInverseFt(outputFunction, inputSpectrum, dx, ctx)
    use M_Utils_Constants
    complex(R64), intent(out)              :: outputFunction(:)
    complex(R64), intent(in)               :: inputSpectrum(:)
    real(R64), intent(in)                  :: dx
    type(T_FftwLib_Ctx), intent(in)       :: ctx

    ! Execute directly on the provided arrays (no copy into ctx%work_in/out)
    call fftw_execute_dft_mod(ctx % plan_backward, inputSpectrum, outputFunction)

    outputFunction = outputFunction / dx / ctx % N

  end subroutine

end module
