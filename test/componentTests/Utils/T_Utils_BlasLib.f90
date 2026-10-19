program T_Utils_BlasLib
  use M_Utils_Types
  use M_Utils_BlasLib
  use testdrive, only: check, error_type
  implicit none

  type(error_type), allocatable :: error

  ! Matrix-vector products with a rectangular matrix for every transpose flag
  call test_mat_on_vec(error)
  if (allocated(error)) error stop 'BlasLib matrix-vector test failed'

contains

  subroutine test_mat_on_vec(error)
    type(error_type), allocatable, intent(out) :: error

    integer(I32), parameter :: nRows = 3, nCols = 5
    complex(R64) :: zA(nRows, nCols), zX(nCols), zXt(nRows), zY(nRows), zYt(nCols)
    real(R64) :: dA(nRows, nCols), dX(nCols), dXt(nRows), dY(nRows), dYt(nCols)
    real(R64) :: maxErr
    integer(I32) :: i, j

    print *, "Testing BlasLib_CalcMatOnVec on a 3x5 matrix..."

    do j = 1, nCols
      do i = 1, nRows
        dA(i, j) = real(i + 2 * j, kind=R64) / 7.0_R64
        zA(i, j) = cmplx(dA(i, j), real(i - j, kind=R64), kind=R64)
      end do
    end do
    do j = 1, nCols
      dX(j) = real(j, kind=R64) - 2.5_R64
      zX(j) = cmplx(dX(j), 0.5_R64 * j, kind=R64)
    end do
    do i = 1, nRows
      dXt(i) = 1.0_R64 / i
      zXt(i) = cmplx(dXt(i), -1.0_R64 * i, kind=R64)
    end do

    maxErr = 0.0_R64

    call BlasLib_CalcMatOnVec(zY, zA, zX, 'N')
    maxErr = max(maxErr, maxval(abs(zY - matmul(zA, zX))))

    call BlasLib_CalcMatOnVec(zYt, zA, zXt, 'T')
    maxErr = max(maxErr, maxval(abs(zYt - matmul(transpose(zA), zXt))))

    call BlasLib_CalcMatOnVec(zYt, zA, zXt, 'C')
    maxErr = max(maxErr, maxval(abs(zYt - matmul(conjg(transpose(zA)), zXt))))

    call BlasLib_CalcMatOnVec(dY, dA, dX, 'N')
    maxErr = max(maxErr, maxval(abs(dY - matmul(dA, dX))))

    call BlasLib_CalcMatOnVec(dYt, dA, dXt, 'T')
    maxErr = max(maxErr, maxval(abs(dYt - matmul(transpose(dA), dXt))))

    print *, "Maximum deviation from matmul: ", maxErr

    call check(error, maxErr < 1.0e-13_R64)

  end subroutine test_mat_on_vec

end program T_Utils_BlasLib
