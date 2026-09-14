! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Utilities for packing/unpacking RDMs with symmetries.
module M_Utils_RdmPacking
  use M_Utils_Types
  implicit none

contains

  !> Calculate the number of floats in the packed representation
  !> Takes into account that diagonal elements store only 1 float.
  pure function RdmPacking_GetFloatCount(nOa) result(nFloats)
    integer(I32), intent(in) :: nOa
    integer(I32) :: nFloats
    integer(I32) :: i1, i2, j1, j2

    nFloats = 0
    do j2 = 1, nOa
      do j1 = 1, nOa
        do i2 = 1, nOa
          do i1 = 1, nOa
            if (isCanonical(i1, i2, j1, j2)) then
              if (isDiagonal(i1, i2, j1, j2)) then
                nFloats = nFloats + 1
              else
                nFloats = nFloats + 2
              end if
            end if
          end do
        end do
      end do
    end do
  end function RdmPacking_GetFloatCount

  !> Pack a Hermitian tensor, omitting imaginary part of diagonals
  subroutine RdmPacking_PackHermitian(packed, m2ab)
    real(R32), intent(out) :: packed(:)
    complex(R64), intent(in) :: m2ab(:, :, :, :)
    integer(I32) :: nOa

    integer(I32) :: i1, i2, j1, j2, idx

    nOa = size(m2ab, 1)
    idx = 1
    do j2 = 1, nOa
      do j1 = 1, nOa
        do i2 = 1, nOa
          do i1 = 1, nOa
            if (isCanonical(i1, i2, j1, j2)) then
              packed(idx) = real(m2ab(i1, i2, j1, j2), R32)
              idx = idx + 1
              if (.not. isDiagonal(i1, i2, j1, j2)) then
                packed(idx) = real(aimag(m2ab(i1, i2, j1, j2)), R32)
                idx = idx + 1
              end if
            end if
          end do
        end do
      end do
    end do
  end subroutine RdmPacking_PackHermitian

  !> Pack an Anti-Hermitian tensor, omitting real part of diagonals
  subroutine RdmPacking_PackAntiHermitian(packed, m2ab)
    real(R32), intent(out) :: packed(:)
    complex(R64), intent(in) :: m2ab(:, :, :, :)
    integer(I32) :: nOa

    integer(I32) :: i1, i2, j1, j2, idx

    nOa = size(m2ab, 1)
    idx = 1
    do j2 = 1, nOa
      do j1 = 1, nOa
        do i2 = 1, nOa
          do i1 = 1, nOa
            if (isCanonical(i1, i2, j1, j2)) then
              if (isDiagonal(i1, i2, j1, j2)) then
                packed(idx) = real(aimag(m2ab(i1, i2, j1, j2)), R32)
                idx = idx + 1
              else
                packed(idx) = real(m2ab(i1, i2, j1, j2), R32)
                packed(idx + 1) = real(aimag(m2ab(i1, i2, j1, j2)), R32)
                idx = idx + 2
              end if
            end if
          end do
        end do
      end do
    end do
  end subroutine RdmPacking_PackAntiHermitian

  !> Unpack a flat array of reals into a 4-index tensor, enforcing Hermitian symmetries
  !> Assumes diagonals only store real part
  subroutine RdmPacking_UnpackHermitian(m2ab, packed)
    complex(R64), intent(out) :: m2ab(:, :, :, :)
    real(R32), intent(in) :: packed(:)
    integer(I32) :: nOa

    integer(I32) :: i1, i2, j1, j2, idx
    complex(R64) :: val

    m2ab = 0.0_R64
    nOa = size(m2ab, 1)
    idx = 1
    do j2 = 1, nOa
      do j1 = 1, nOa
        do i2 = 1, nOa
          do i1 = 1, nOa
            if (isCanonical(i1, i2, j1, j2)) then
              if (isDiagonal(i1, i2, j1, j2)) then
                val = cmplx(packed(idx), 0.0_R64, R64)
                idx = idx + 1
              else
                val = cmplx(packed(idx), packed(idx + 1), R64)
                idx = idx + 2
              end if

              ! Assign to all symmetric equivalents
              m2ab(i1, i2, j1, j2) = val
              m2ab(i2, i1, j2, j1) = val

              m2ab(j1, j2, i1, i2) = conjg(val)
              m2ab(j2, j1, i2, i1) = conjg(val)
            end if
          end do
        end do
      end do
    end do
  end subroutine RdmPacking_UnpackHermitian

  !> Unpack a flat array of reals into a 4-index tensor, enforcing anti-Hermitian symmetries
  subroutine RdmPacking_UnpackAntiHermitian(m2ab, packed)
    complex(R64), intent(out) :: m2ab(:, :, :, :)
    real(R32), intent(in) :: packed(:)
    integer(I32) :: nOa

    integer(I32) :: i1, i2, j1, j2, idx
    complex(R64) :: val

    m2ab = 0.0_R64
    nOa = size(m2ab, 1)
    idx = 1
    do j2 = 1, nOa
      do j1 = 1, nOa
        do i2 = 1, nOa
          do i1 = 1, nOa
            if (isCanonical(i1, i2, j1, j2)) then
              if (isDiagonal(i1, i2, j1, j2)) then
                val = cmplx(0.0_R64, packed(idx), R64)
                idx = idx + 1
              else
                val = cmplx(packed(idx), packed(idx + 1), R64)
                idx = idx + 2
              end if

              ! Assign to all symmetric equivalents
              m2ab(i1, i2, j1, j2) = val
              m2ab(i2, i1, j2, j1) = val

              m2ab(j1, j2, i1, i2) = -conjg(val)
              m2ab(j2, j1, i2, i1) = -conjg(val)
            end if
          end do
        end do
      end do
    end do
  end subroutine RdmPacking_UnpackAntiHermitian

  !> Check if indices are the canonical representative of their orbit
  pure function isCanonical(i1, i2, j1, j2) result(isCan)
    integer(I32), intent(in) :: i1, i2, j1, j2
    logical :: isCan

    ! Symmetry 1: (i1, i2, j1, j2) ~ (i2, i1, j2, j1)
    ! We choose the one that is lexicographically smaller
    if (compareIndices(i1, i2, j1, j2, i2, i1, j2, j1) > 0) then
      isCan = .false.
      return
    end if

    ! Symmetry 2: (i1, i2, j1, j2) ~ (j1, j2, i1, i2)*
    if (compareIndices(i1, i2, j1, j2, j1, j2, i1, i2) > 0) then
      isCan = .false.
      return
    end if

    ! Combined: (i1, i2, j1, j2) ~ (j2, j1, i2, i1)*
    if (compareIndices(i1, i2, j1, j2, j2, j1, i2, i1) > 0) then
      isCan = .false.
      return
    end if

    isCan = .true.
  end function isCanonical

  pure function compareIndices(i1, i2, j1, j2, k1, k2, l1, l2) result(cmp)
    integer(I32), intent(in) :: i1, i2, j1, j2, k1, k2, l1, l2
    integer(I32) :: cmp

    if (j2 < l2) then
      cmp = -1
      return
    else if (j2 > l2) then
      cmp = 1
      return
    end if

    if (j1 < l1) then
      cmp = -1
      return
    else if (j1 > l1) then
      cmp = 1
      return
    end if

    if (i2 < k2) then
      cmp = -1
      return
    else if (i2 > k2) then
      cmp = 1
      return
    end if

    if (i1 < k1) then
      cmp = -1
      return
    else if (i1 > k1) then
      cmp = 1
      return
    end if

    cmp = 0
  end function compareIndices

  pure function isDiagonal(i1, i2, j1, j2) result(isDiag)
    integer(I32), intent(in) :: i1, i2, j1, j2
    logical :: isDiag
    isDiag = (i1 == j1 .and. i2 == j2) .or. (i1 == j2 .and. i2 == j1)
  end function isDiagonal

end module M_Utils_RdmPacking
