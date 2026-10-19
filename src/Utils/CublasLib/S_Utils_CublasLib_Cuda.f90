! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> cuBLAS implementation of M_Utils_CublasLib (CMake option CODY_USE_CUDA=ON):
!> the CUDA runtime and cuBLAS C APIs are called through bind(C) interfaces.
submodule(M_Utils_CublasLib) S_Utils_CublasLib_Cuda
  use, intrinsic :: iso_c_binding, only: c_int, c_size_t, c_intptr_t, c_double, c_null_ptr, &
                                         c_associated, c_loc

  implicit none

  !=============================================================================
  ! local data
  !=============================================================================

  !> cuBLAS handle, created on the first `CublasLib_AvailableQ`.
  type(c_ptr) :: handle = c_null_ptr

  !> cudaMemcpyKind and cublasOperation_t values.
  integer(c_int), parameter :: memcpyHostToDevice = 1
  integer(c_int), parameter :: memcpyDeviceToHost = 2
  integer(c_int), parameter :: opN = 0
  integer(c_int), parameter :: opT = 1

  !> Bytes per real(R64) element.
  integer(I64), parameter :: elementBytes = 8

  interface
    function cudaGetDeviceCount(count) result(res) bind(C, name="cudaGetDeviceCount")
      import :: c_int
      integer(c_int), intent(out) :: count
      integer(c_int) :: res
    end function

    function cudaMalloc(ptr, size) result(res) bind(C, name="cudaMalloc")
      import :: c_ptr, c_size_t, c_int
      type(c_ptr), intent(out) :: ptr
      integer(c_size_t), value :: size
      integer(c_int) :: res
    end function

    function cudaFree(ptr) result(res) bind(C, name="cudaFree")
      import :: c_ptr, c_int
      type(c_ptr), value :: ptr
      integer(c_int) :: res
    end function

    function cudaMemcpy(dst, src, count, kind) result(res) bind(C, name="cudaMemcpy")
      import :: c_ptr, c_size_t, c_int
      type(c_ptr), value :: dst
      type(c_ptr), value :: src
      integer(c_size_t), value :: count
      integer(c_int), value :: kind
      integer(c_int) :: res
    end function

    function cublasCreate(handle) result(res) bind(C, name="cublasCreate_v2")
      import :: c_ptr, c_int
      type(c_ptr), intent(out) :: handle
      integer(c_int) :: res
    end function

    function cublasDgemm(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc) &
      result(res) bind(C, name="cublasDgemm_v2")
      import :: c_ptr, c_int, c_double
      type(c_ptr), value :: handle
      integer(c_int), value :: transa
      integer(c_int), value :: transb
      integer(c_int), value :: m
      integer(c_int), value :: n
      integer(c_int), value :: k
      real(c_double), intent(in) :: alpha
      type(c_ptr), value :: A
      integer(c_int), value :: lda
      type(c_ptr), value :: B
      integer(c_int), value :: ldb
      real(c_double), intent(in) :: beta
      type(c_ptr), value :: C
      integer(c_int), value :: ldc
      integer(c_int) :: res
    end function
  end interface

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module function CublasLib_AvailableQ() result(res)
    logical :: res

    integer(c_int) :: nDevices

    res = .false.
    if (c_associated(handle)) then
      res = .true.
      return
    end if

    if (cudaGetDeviceCount(nDevices) .ne. 0) return
    if (nDevices < 1) return
    if (cublasCreate(handle) .ne. 0) return

    res = .true.

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine CublasLib_Allocate(dPtr, nElements)
    type(c_ptr), intent(out) :: dPtr
    integer(I64), intent(in) :: nElements

    if (cudaMalloc(dPtr, int(nElements * elementBytes, kind=c_size_t)) .ne. 0) &
      error stop "CublasLib_Allocate: cudaMalloc failed (device memory exhausted?)"

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine CublasLib_Deallocate(dPtr)
    type(c_ptr), intent(inout) :: dPtr

    if (.not. c_associated(dPtr)) return
    if (cudaFree(dPtr) .ne. 0) error stop "CublasLib_Deallocate: cudaFree failed"
    dPtr = c_null_ptr

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine CublasLib_CopyToDevice(dPtr, hostData, offset)
    type(c_ptr), intent(in) :: dPtr
    real(R64), intent(in), contiguous, target :: hostData(:)
    integer(I64), intent(in) :: offset

    if (cudaMemcpy(Shifted(dPtr, offset), c_loc(hostData), &
                   int(size(hostData, kind=I64) * elementBytes, kind=c_size_t), memcpyHostToDevice) .ne. 0) &
      error stop "CublasLib_CopyToDevice: cudaMemcpy failed"

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine CublasLib_CopyToHost(hostData, dPtr, offset)
    real(R64), intent(out), contiguous, target :: hostData(:)
    type(c_ptr), intent(in) :: dPtr
    integer(I64), intent(in) :: offset

    if (cudaMemcpy(c_loc(hostData), Shifted(dPtr, offset), &
                   int(size(hostData, kind=I64) * elementBytes, kind=c_size_t), memcpyDeviceToHost) .ne. 0) &
      error stop "CublasLib_CopyToHost: cudaMemcpy failed"

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

    if (cublasDgemm(handle, Operation(transA), Operation(transB), m, n, k, alpha, &
                    Shifted(dA, offsetA), lda, Shifted(dB, offsetB), ldb, beta, &
                    Shifted(dC, offsetC), ldc) .ne. 0) &
      error stop "CublasLib_Dgemm: cublasDgemm failed"

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Device address `offset` real(R64) elements behind `dPtr`.
  function Shifted(dPtr, offset) result(res)
    type(c_ptr), intent(in) :: dPtr
    integer(I64), intent(in) :: offset
    type(c_ptr) :: res

    integer(c_intptr_t) :: address

    address = transfer(dPtr, address)
    address = address + int(offset * elementBytes, kind=c_intptr_t)
    res = transfer(address, res)

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> cublasOperation_t of a BLAS transpose flag ('N' or 'T').
  function Operation(trans) result(res)
    character(len=1), intent(in) :: trans
    integer(c_int) :: res

    select case (trans)
    case ('N', 'n')
      res = opN
    case ('T', 't')
      res = opT
    case default
      error stop "CublasLib_Dgemm: transpose flag must be 'N' or 'T'"
    end select

  end function

end submodule
