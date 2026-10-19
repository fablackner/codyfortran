! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Minimal GPU linear algebra through cuBLAS (double precision).
!>
!> Device memory is addressed by `type(c_ptr)` handles plus element offsets
!> (real(R64) elements), matrices are column-major as in Fortran/BLAS. The
!> implementation submodule is chosen at build time: S_Utils_CublasLib_Cuda
!> (CMake option CODY_USE_CUDA=ON, links cuBLAS and the CUDA runtime through
!> bind(C), no CUDA compiler needed) or S_Utils_CublasLib_Stub (default; no GPU,
!> every call except `CublasLib_AvailableQ` stops).
module M_Utils_CublasLib
  use, intrinsic :: iso_c_binding, only: c_ptr
  use M_Utils_Types

  implicit none

  interface
    !> True if a CUDA device is usable; initializes cuBLAS on the first call.
    module function CublasLib_AvailableQ() result(res)
      logical :: res
    end function

    !> Allocates `nElements` real(R64) elements of device memory.
    module subroutine CublasLib_Allocate(dPtr, nElements)
      !> Device memory handle.
      type(c_ptr), intent(out) :: dPtr
      !> Number of real(R64) elements.
      integer(I64), intent(in) :: nElements
    end subroutine

    !> Releases device memory allocated by `CublasLib_Allocate`.
    module subroutine CublasLib_Deallocate(dPtr)
      !> Device memory handle; null on exit.
      type(c_ptr), intent(inout) :: dPtr
    end subroutine

    !> Copies host data to device memory starting at element `offset`.
    module subroutine CublasLib_CopyToDevice(dPtr, hostData, offset)
      !> Device memory handle.
      type(c_ptr), intent(in) :: dPtr
      !> Host data (all elements are copied).
      real(R64), intent(in), contiguous, target :: hostData(:)
      !> Element offset into the device memory.
      integer(I64), intent(in) :: offset
    end subroutine

    !> Copies device memory starting at element `offset` to host data.
    module subroutine CublasLib_CopyToHost(hostData, dPtr, offset)
      !> Host data (size(hostData) elements are copied).
      real(R64), intent(out), contiguous, target :: hostData(:)
      !> Device memory handle.
      type(c_ptr), intent(in) :: dPtr
      !> Element offset into the device memory.
      integer(I64), intent(in) :: offset
    end subroutine

    !> C = alpha op(A) op(B) + beta C on the device (cublasDgemm), with A, B, C
    !> starting at the given element offsets of their device memory.
    module subroutine CublasLib_Dgemm(transA, transB, m, n, k, alpha, dA, offsetA, lda, &
                                      dB, offsetB, ldb, beta, dC, offsetC, ldc)
      character(len=1), intent(in) :: transA
      character(len=1), intent(in) :: transB
      integer(I32), intent(in) :: m
      integer(I32), intent(in) :: n
      integer(I32), intent(in) :: k
      real(R64), intent(in) :: alpha
      type(c_ptr), intent(in) :: dA
      integer(I64), intent(in) :: offsetA
      integer(I32), intent(in) :: lda
      type(c_ptr), intent(in) :: dB
      integer(I64), intent(in) :: offsetB
      integer(I32), intent(in) :: ldb
      real(R64), intent(in) :: beta
      type(c_ptr), intent(in) :: dC
      integer(I64), intent(in) :: offsetC
      integer(I32), intent(in) :: ldc
    end subroutine
  end interface

end module
