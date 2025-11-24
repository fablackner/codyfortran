! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Manual Fourier-transform utilities that mirror the conventions used by FFTW.
!>
!> Provides helpers to build raw FFTW-order wavenumber arrays and to perform
!> forward/inverse discrete Fourier transforms consistent with the convention:
!>   F(k_j) = Σ f(x_n) exp(-i k_j x_n) dx / sqrt(2π)
!>   f(x_n) = Σ F(k_j) exp(+i k_j x_n) dk / sqrt(2π)
!> where `k_(:)` is in raw FFTW order (wrapped negative frequencies).
module M_Utils_FtManual
  use M_Utils_Types

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Build raw FFTW-order wavenumbers for a uniform real-space grid.
!>
!> Parameters
!> - kArray [out]: real array of size N receiving k_j with raw FFTW ordering.
!> - dx     [in]:  grid spacing in real space; sets dk = 2π/(N*dx).
  pure subroutine BuildRawWaveNumbers(kArray, dx)
    use M_Utils_Constants

    ! Build raw FFTW-order wavenumber array for a uniform grid.
    real(R64), intent(out) :: kArray(:)
    real(R64), intent(in)  :: dx
    integer(I32) :: N, j, m
    real(R64) :: dk

    N = size(kArray)
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
!> Forward discrete Fourier transform consistent with the module’s convention.
!>
!> Parameters
!> - outputSpectrum [out]: complex spectrum F(k_j) in raw FFTW k-order.
!> - dx             [in]:  real-space grid spacing used in normalization.
!> - inputFunction  [in]:  real-space samples f(x_n).
!> - k_             [in,opt]: precomputed k-array (raw FFTW order); built internally if absent.
  subroutine FtManual_DoFt(outputSpectrum, dx, inputFunction, k_)
    use M_Utils_Constants
    ! Signature now matches FtManual_DoFt
    complex(R64), intent(out), allocatable :: outputSpectrum(:)
    real(R64), intent(in)               :: dx
    complex(R64), intent(in)               :: inputFunction(:)
    real(R64), intent(out), allocatable, optional :: k_(:)

    integer(I32) :: N, jK, jX
    complex(R64) :: sumVal
    real(R64), allocatable :: xLocal(:), kLocal(:)

    N = size(inputFunction)
    if (.not. allocated(outputSpectrum)) allocate (outputSpectrum(N))
    if (size(outputSpectrum) .ne. N) stop 'FtManual_DoFtManual: size(outputSpectrum) mismatch'

    if (present(k_)) then
      if (.not. allocated(k_)) allocate (k_(N))
      if (size(k_) .ne. N) stop 'FtManual_DoFtManual: size(k_) mismatch'
      call BuildRawWaveNumbers(k_, dx)
    else
      allocate (kLocal(N))
      call BuildRawWaveNumbers(kLocal, dx)
    end if

    allocate (xLocal(N))
    xLocal = [((jX - 1) * dx, jX=1, N)]

    do jK = 1, N
      sumVal = (0.0_R64, 0.0_R64)
      do jX = 1, N
        if (present(k_)) then
          sumVal = sumVal + inputFunction(jX) * exp(-IU * k_(jK) * xLocal(jX))
        else
          sumVal = sumVal + inputFunction(jX) * exp(-IU * kLocal(jK) * xLocal(jX))
        end if
      end do
      outputSpectrum(jK) = sumVal * dx
    end do

    deallocate (xLocal)
    if (.not. present(k_)) deallocate (kLocal)
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Inverse discrete Fourier transform consistent with the module’s convention.
!>
!> Parameters
!> - outputFunction [out]: real-space reconstruction f(x_n).
!> - dx             [in]:  real-space spacing (used to compute dk and scaling).
!> - inputSpectrum  [in]:  complex spectrum F(k_j) in raw FFTW k-order.
  subroutine FtManual_DoInverseFt(outputFunction, dx, inputSpectrum)
    use M_Utils_Constants
    ! Signature now matches FtManual_DoInverseFt
    complex(R64), intent(out), allocatable :: outputFunction(:)
    real(R64), intent(in)               :: dx
    complex(R64), intent(in)               :: inputSpectrum(:)  ! raw ordering

    integer(I32) :: N, jX, jK
    real(R64) :: dk
    complex(R64) :: sumVal
    real(R64), allocatable :: xLocal(:), kLocal(:)

    N = size(inputSpectrum)
    if (.not. allocated(outputFunction)) allocate (outputFunction(N))
    if (size(outputFunction) .ne. N) stop 'FtManual_DoManualInverseFT: size(outputFunction) mismatch'

    if (N > 0) dk = 2.0_R64 * PI / (N * dx)

    allocate (xLocal(N), kLocal(N))
    xLocal = [((jX - 1) * dx, jX=1, N)]
    call BuildRawWaveNumbers(kLocal, dx)

    do jX = 1, N
      sumVal = (0.0_R64, 0.0_R64)
      do jK = 1, N
        sumVal = sumVal + inputSpectrum(jK) * exp(IU * kLocal(jK) * xLocal(jX))
      end do
      outputFunction(jX) = sumVal / (dx * N)
    end do

    deallocate (xLocal, kLocal)
  end subroutine

end module
