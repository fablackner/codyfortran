! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
module M_Utils_RdmCholeskyExtract
  use M_Utils_Types

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmCholeskyExtract_FillCholeskyVectorFromRdm3aab(choleskyVector, rdm3aab)
    !> Converts a 3-RDM with (aab,aab) spin signature to its Cholesky
    !> decomposition factors, returned as a linear vector of the lower
    !> triangular matrix.
    !> Assumes antisymmetry in the first two alpha indices: i1<i2, j1<j2.
    use M_Utils_LapackLib
    implicit none

    complex(R64), intent(out), allocatable :: choleskyVector(:)
    complex(R64), intent(in), contiguous :: rdm3aab(:, :, :, :, :, :)

    integer(I32) :: dim, i1, i2, i3, j1, j2, j3, i, j, k
    integer(I32) :: nOa, nOb
    complex(R64), allocatable :: mMat(:, :)

    nOa = size(rdm3aab, 1)
    nOb = size(rdm3aab, 3)

    dim = (nOa * (nOa - 1) / 2) * nOb
    allocate (mMat(dim, dim))

    ! Reshape 3-RDM into a 2D matrix mMat
    ! i = (i1,i2,i3) with i1<i2
    ! j = (j1,j2,j3) with j1<j2
    i = 0
    do i3 = 1, nOb
      do i2 = 2, nOa
        do i1 = 1, i2 - 1
          i = i + 1
          j = 0
          do j3 = 1, nOb
            do j2 = 2, nOa
              do j1 = 1, j2 - 1
                j = j + 1
                mMat(i, j) = rdm3aab(i1, i2, i3, j1, j2, j3)
              end do
            end do
          end do
        end do
      end do
    end do

    ! Cholesky decomposition mMat = L * L^H
    call LapackLib_FactorizeCholesky(mMat, 'L')

    ! Extract lower triangular elements into a vector
    allocate (choleskyVector(dim * (dim + 1) / 2))
    k = 0
    do j = 1, dim
      do i = j, dim
        k = k + 1
        choleskyVector(k) = mMat(i, j)
      end do
    end do

    deallocate (mMat)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmCholeskyExtract_FillRdm3aabFromCholeskyVector(rdm3aab, choleskyVector)
    !> Reconstructs a 3-RDM with (aab,aab) spin signature from a
    !> linear vector of its Cholesky decomposition factors.
    implicit none

    complex(R64), intent(out), contiguous :: rdm3aab(:, :, :, :, :, :)
    complex(R64), intent(in), contiguous :: choleskyVector(:)

    integer(I32) :: dim, i1, i2, i3, j1, j2, j3, i, j, k
    complex(R64), allocatable :: mMat(:, :), lMat(:, :)
    integer(I32) :: nOa, nOb
    integer(I32) :: choleskyVectorSize

    nOa = size(rdm3aab, 1)
    nOb = size(rdm3aab, 3)

    dim = (nOa * (nOa - 1) / 2) * nOb
    choleskyVectorSize = dim * (dim + 1) / 2
    if (size(choleskyVector) .ne. choleskyVectorSize) then
      ! Handle error: input vector size is incorrect
      return
    end if

    allocate (lMat(dim, dim), mMat(dim, dim))

    ! Zero diagonal alpha-pair blocks:
    do k = 1, nOa
      rdm3aab(k, k, :, :, :, :) = (0.0_R64, 0.0_R64)
      rdm3aab(:, :, :, k, k, :) = (0.0_R64, 0.0_R64)
    end do

    ! Reconstruct lower triangular matrix L from vector
    lMat = (0.0_R64, 0.0_R64)
    k = 0
    do j = 1, dim
      do i = j, dim
        k = k + 1
        lMat(i, j) = choleskyVector(k)
      end do
    end do

    ! Reconstruct the matrix mMat = L * L^H
    mMat = matmul(lMat, conjg(transpose(lMat)))

    ! Reshape mMat back to 3-RDM
    i = 0
    do i3 = 1, nOb
      do i2 = 2, nOa
        do i1 = 1, i2 - 1
          i = i + 1
          j = 0
          do j3 = 1, nOb
            do j2 = 2, nOa
              do j1 = 1, j2 - 1
                j = j + 1
                rdm3aab(i1, i2, i3, j1, j2, j3) = mMat(i, j)
                ! Antisymmetry in first two alpha indices (bra)
                rdm3aab(i2, i1, i3, j1, j2, j3) = -mMat(i, j)
                ! Antisymmetry in first two alpha indices (ket)
                rdm3aab(i1, i2, i3, j2, j1, j3) = -mMat(i, j)
                ! Both swapped
                rdm3aab(i2, i1, i3, j2, j1, j3) = mMat(i, j)
              end do
            end do
          end do
        end do
      end do
    end do

    deallocate (mMat, lMat)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmCholeskyExtract_FillCholeskyVectorFromRdm2ab(choleskyVector, rdm2ab)
    !> Converts a 2-RDM with (ab,ab) spin signature to its Cholesky
    !> decomposition factors, returned as a linear vector of the lower
    !> triangular matrix.
    use M_Utils_LapackLib
    implicit none
    complex(R64), intent(out), allocatable :: choleskyVector(:)
    complex(R64), intent(in), contiguous :: rdm2ab(:, :, :, :)
    integer(I32) :: dim, i1, i2, j1, j2, i, j, k
    integer(I32) :: nOa, nOb
    complex(R64), allocatable :: mMat(:, :)

    nOa = size(rdm2ab, 1)
    nOb = size(rdm2ab, 2)
    dim = nOa * nOb
    allocate (mMat(dim, dim))

    ! Reshape 2-RDM into a 2D matrix mMat
    i = 0
    do i2 = 1, nOb
      do i1 = 1, nOa
        i = i + 1
        j = 0
        do j2 = 1, nOb
          do j1 = 1, nOa
            j = j + 1
            mMat(i, j) = rdm2ab(i1, i2, j1, j2)
          end do
        end do
      end do
    end do

    ! Cholesky decomposition mMat = L * L^H
    call LapackLib_FactorizeCholesky(mMat, 'L')

    ! Extract lower triangular elements into a vector
    allocate (choleskyVector(dim * (dim + 1) / 2))
    k = 0
    do j = 1, dim
      do i = j, dim
        k = k + 1
        choleskyVector(k) = mMat(i, j)
      end do
    end do

    deallocate (mMat)
  end subroutine

  subroutine RdmCholeskyExtract_FillRdm2abFromCholeskyVector(rdm2ab, choleskyVector)
    !> Reconstructs a 2-RDM with (ab,ab) spin signature from a
    !> linear vector of its Cholesky decomposition factors.
    implicit none
    complex(R64), intent(out), contiguous :: rdm2ab(:, :, :, :)
    complex(R64), intent(in), contiguous :: choleskyVector(:)
    integer(I32) :: dim, i1, i2, j1, j2, i, j, k
    complex(R64), allocatable :: mMat(:, :), lMat(:, :)
    integer(I32) :: nOa, nOb
    integer(I32) :: choleskyVectorSize

    nOa = size(rdm2ab, 1)
    nOb = size(rdm2ab, 2)
    dim = nOa * nOb
    choleskyVectorSize = dim * (dim + 1) / 2
    if (size(choleskyVector) .ne. choleskyVectorSize) then
      return
    end if

    allocate (lMat(dim, dim), mMat(dim, dim))

    ! Reconstruct lower triangular matrix L from vector
    lMat = (0.0_R64, 0.0_R64)
    k = 0
    do j = 1, dim
      do i = j, dim
        k = k + 1
        lMat(i, j) = choleskyVector(k)
      end do
    end do

    ! Reconstruct the matrix mMat = L * L^H
    mMat = matmul(lMat, conjg(transpose(lMat)))

    ! Reshape mMat back to 2-RDM
    i = 0
    do i2 = 1, nOb
      do i1 = 1, nOa
        i = i + 1
        j = 0
        do j2 = 1, nOb
          do j1 = 1, nOa
            j = j + 1
            rdm2ab(i1, i2, j1, j2) = mMat(i, j)
          end do
        end do
      end do
    end do

    deallocate (mMat, lMat)
  end subroutine

end module
