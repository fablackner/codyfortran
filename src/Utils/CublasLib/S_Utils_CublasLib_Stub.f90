! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Stub implementation of M_Utils_CublasLib for builds without CUDA (default,
!> CMake option CODY_USE_CUDA=OFF): no device is available and every other
!> call stops.
submodule(M_Utils_CublasLib) S_Utils_CublasLib_Stub
  use, intrinsic :: iso_c_binding, only: c_null_ptr
  use M_Utils_UnusedVariables

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module function CublasLib_AvailableQ() result(res)
    logical :: res

    res = .false.

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine CublasLib_Allocate(dPtr, nElements)
    type(c_ptr), intent(out) :: dPtr
    integer(I64), intent(in) :: nElements

    if (.false.) call UnusedVariables_Mark(nElements)
    dPtr = c_null_ptr
    error stop "CublasLib: CodyFortran was built without CUDA (cmake -DCODY_USE_CUDA=ON)"

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine CublasLib_Deallocate(dPtr)
    type(c_ptr), intent(inout) :: dPtr

    dPtr = c_null_ptr

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine CublasLib_CopyToDevice(dPtr, hostData, offset)
    type(c_ptr), intent(in) :: dPtr
    real(R64), intent(in), contiguous, target :: hostData(:)
    integer(I64), intent(in) :: offset

    if (.false.) call UnusedVariables_Mark(dPtr, hostData, offset)
    error stop "CublasLib: CodyFortran was built without CUDA (cmake -DCODY_USE_CUDA=ON)"

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine CublasLib_CopyToHost(hostData, dPtr, offset)
    real(R64), intent(out), contiguous, target :: hostData(:)
    type(c_ptr), intent(in) :: dPtr
    integer(I64), intent(in) :: offset

    if (.false.) call UnusedVariables_Mark(dPtr, offset)
    hostData(:) = 0.0_R64
    error stop "CublasLib: CodyFortran was built without CUDA (cmake -DCODY_USE_CUDA=ON)"

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

    if (.false.) call UnusedVariables_Mark(transA, transB, m, n, k)
    if (.false.) call UnusedVariables_Mark(alpha, dA, offsetA, lda, beta)
    if (.false.) call UnusedVariables_Mark(dB, offsetB, ldb)
    if (.false.) call UnusedVariables_Mark(dC, offsetC, ldc)
    error stop "CublasLib: CodyFortran was built without CUDA (cmake -DCODY_USE_CUDA=ON)"

  end subroutine

end submodule
