program T_Utils_CublasLib
  use, intrinsic :: iso_c_binding, only: c_ptr
  use M_Utils_Types
  use M_Utils_CublasLib
  use testdrive, only: check, error_type
  implicit none

  type(error_type), allocatable :: error

  ! Device GEMM with offsets and both transpose flags against matmul; skipped
  ! without a CUDA device or in builds without -DCODY_USE_CUDA=ON
  if (.not. CublasLib_AvailableQ()) then
    print *, "No CUDA device available (or built without CUDA): test skipped"
    stop
  end if

  call test_dgemm(error)
  if (allocated(error)) error stop 'CublasLib GEMM test failed'

contains

  subroutine test_dgemm(error)
    type(error_type), allocatable, intent(out) :: error

    integer(I32), parameter :: n = 7, nCols = 3
    real(R64) :: a(n, n), b(n, nCols), c(n, nCols), packed(n * n + n * nCols)
    real(R64) :: maxErr
    type(c_ptr) :: dAB, dC
    integer(I32) :: i, j

    print *, "Testing CublasLib_Dgemm on a 7x7 times 7x3 product..."

    do j = 1, n
      do i = 1, n
        a(i, j) = real(i - 2 * j, kind=R64) / 5.0_R64
      end do
    end do
    do j = 1, nCols
      do i = 1, n
        b(i, j) = 1.0_R64 / real(i + j, kind=R64)
      end do
    end do

    ! A and B in one device array (B at an element offset)
    packed(1:n * n) = reshape(a, [n * n])
    packed(n * n + 1:) = reshape(b, [n * nCols])
    call CublasLib_Allocate(dAB, int(size(packed), kind=I64))
    call CublasLib_Allocate(dC, int(n * nCols, kind=I64))
    call CublasLib_CopyToDevice(dAB, packed, 0_I64)

    maxErr = 0.0_R64

    call CublasLib_Dgemm('N', 'N', n, nCols, n, 1.0_R64, dAB, 0_I64, n, dAB, int(n * n, kind=I64), n, &
                         0.0_R64, dC, 0_I64, n)
    call CopyBack(c, dC)
    maxErr = max(maxErr, maxval(abs(c - matmul(a, b))))

    call CublasLib_Dgemm('T', 'N', n, nCols, n, 1.0_R64, dAB, 0_I64, n, dAB, int(n * n, kind=I64), n, &
                         0.0_R64, dC, 0_I64, n)
    call CopyBack(c, dC)
    maxErr = max(maxErr, maxval(abs(c - matmul(transpose(a), b))))

    call CublasLib_Deallocate(dAB)
    call CublasLib_Deallocate(dC)

    print *, "Maximum deviation from matmul: ", maxErr

    call check(error, maxErr < 1.0e-13_R64)

  end subroutine test_dgemm

  !> Copies the device matrix dC back into cOut (same shape).
  subroutine CopyBack(cOut, dC)
    real(R64), intent(out) :: cOut(:, :)
    type(c_ptr), intent(in) :: dC

    real(R64), allocatable :: flat(:)

    allocate (flat(size(cOut)))
    call CublasLib_CopyToHost(flat, dC, 0_I64)
    cOut = reshape(flat, shape(cOut))

  end subroutine

end program T_Utils_CublasLib
