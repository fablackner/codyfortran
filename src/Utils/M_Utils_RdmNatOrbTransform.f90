! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
module M_Utils_RdmNatOrbTransform
  use M_Utils_Types

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmNatOrbTransform_CalculateNatOrbs(natorbs, rdm1)
    use M_Utils_RdmDiagonalize

    complex(R64), intent(out), allocatable :: natorbs(:, :)
    complex(R64), intent(in), contiguous :: rdm1(:, :)
    integer(I32) :: nO
    real(R64), allocatable :: natocc(:)

    nO = size(rdm1, 1)
    allocate (natocc(nO))
    allocate (natorbs(nO, nO))

    call RdmDiagonalize_Rdm1(natocc, natorbs, rdm1)
  end subroutine

  subroutine RdmNatOrbTransform_TransformRdm1(rdm1, natorbs)
    complex(R64), intent(inout), contiguous :: rdm1(:, :)
    complex(R64), intent(in), contiguous :: natorbs(:, :)
    integer(I32) :: nO
    complex(R64), allocatable :: temp1(:, :), natorbs_H(:, :)

    nO = size(rdm1, 1)
    allocate (natorbs_H(nO, nO))
    natorbs_H = conjg(transpose(natorbs))

    allocate (temp1(nO, nO))
    temp1 = matmul(rdm1, natorbs)
    rdm1 = matmul(natorbs_H, temp1)
    deallocate (temp1, natorbs_H)
  end subroutine

  subroutine RdmNatOrbTransform_TransformRdm2(rdm2, natorbs)
    complex(R64), intent(inout), contiguous :: rdm2(:, :, :, :)
    complex(R64), intent(in), contiguous :: natorbs(:, :)
    integer(I32) :: nO, i1, i2, j1, j2, k
    complex(R64), allocatable :: temp2(:, :, :, :), natorbs_H(:, :)

    nO = size(rdm2, 1)
    allocate (natorbs_H(nO, nO))
    natorbs_H = conjg(transpose(natorbs))

    allocate (temp2(nO, nO, nO, nO))
    temp2 = 0.0_R64
    do j1 = 1, nO; do i2 = 1, nO; do i1 = 1, nO; do k = 1, nO
        temp2(i1, i2, j1, :) = temp2(i1, i2, j1, :) + rdm2(i1, i2, j1, k) * natorbs(k, :)
      end do; end do; end do; end do
    rdm2 = 0.0_R64
    do i2 = 1, nO; do i1 = 1, nO; do j2 = 1, nO; do k = 1, nO
        rdm2(i1, i2, :, j2) = rdm2(i1, i2, :, j2) + temp2(i1, i2, k, j2) * natorbs(k, :)
      end do; end do; end do; end do
    temp2 = 0.0_R64
    do i1 = 1, nO; do j1 = 1, nO; do j2 = 1, nO; do k = 1, nO
        temp2(i1, :, j1, j2) = temp2(i1, :, j1, j2) + natorbs_H(:, k) * rdm2(i1, k, j1, j2)
      end do; end do; end do; end do
    rdm2 = 0.0_R64
    do i2 = 1, nO; do j1 = 1, nO; do j2 = 1, nO; do k = 1, nO
        rdm2(:, i2, j1, j2) = rdm2(:, i2, j1, j2) + natorbs_H(:, k) * temp2(k, i2, j1, j2)
      end do; end do; end do; end do
    deallocate (temp2, natorbs_H)
  end subroutine

  subroutine RdmNatOrbTransform_InvTransformRdm2(rdm2, natorbs)
    complex(R64), intent(inout), contiguous :: rdm2(:, :, :, :)
    complex(R64), intent(in), contiguous :: natorbs(:, :)
    integer(I32) :: nO, i1, i2, j1, j2, k
    complex(R64), allocatable :: temp2(:, :, :, :), natorbs_H(:, :)

    nO = size(rdm2, 1)
    allocate (natorbs_H(nO, nO))
    natorbs_H = conjg(transpose(natorbs))

    allocate (temp2(nO, nO, nO, nO))
    temp2 = 0.0_R64
    do j1 = 1, nO; do i2 = 1, nO; do i1 = 1, nO; do k = 1, nO
        temp2(i1, i2, j1, :) = temp2(i1, i2, j1, :) + rdm2(i1, i2, j1, k) * natorbs_H(k, :)
      end do; end do; end do; end do
    rdm2 = 0.0_R64
    do i2 = 1, nO; do i1 = 1, nO; do j2 = 1, nO; do k = 1, nO
        rdm2(i1, i2, :, j2) = rdm2(i1, i2, :, j2) + temp2(i1, i2, k, j2) * natorbs_H(k, :)
      end do; end do; end do; end do
    temp2 = 0.0_R64
    do i1 = 1, nO; do j1 = 1, nO; do j2 = 1, nO; do k = 1, nO
        temp2(i1, :, j1, j2) = temp2(i1, :, j1, j2) + natorbs(:, k) * rdm2(i1, k, j1, j2)
      end do; end do; end do; end do
    rdm2 = 0.0_R64
    do i2 = 1, nO; do j1 = 1, nO; do j2 = 1, nO; do k = 1, nO
        rdm2(:, i2, j1, j2) = rdm2(:, i2, j1, j2) + natorbs(:, k) * temp2(k, i2, j1, j2)
      end do; end do; end do; end do
    deallocate (temp2, natorbs_H)
  end subroutine

  subroutine RdmNatOrbTransform_TransformRdm3(rdm3, natorbs)
    complex(R64), intent(inout), contiguous :: rdm3(:, :, :, :, :, :)
    complex(R64), intent(in), contiguous :: natorbs(:, :)
    integer(I32) :: nO, i1, i2, i3, j1, j2, j3, k
    complex(R64), allocatable :: temp3(:, :, :, :, :, :), natorbs_H(:, :)

    nO = size(rdm3, 1)
    allocate (natorbs_H(nO, nO))
    natorbs_H = conjg(transpose(natorbs))

    allocate (temp3(nO, nO, nO, nO, nO, nO))
    temp3 = 0.0_R64
    do j2 = 1, nO; do j1 = 1, nO; do i3 = 1, nO; do i2 = 1, nO; do i1 = 1, nO; do k = 1, nO
        temp3(i1, i2, i3, j1, j2, :) = temp3(i1, i2, i3, j1, j2, :) + rdm3(i1, i2, i3, j1, j2, k) * natorbs(k, :)
      end do; end do; end do; end do; end do; end do
    rdm3 = 0.0_R64
    do j1 = 1, nO; do i3 = 1, nO; do i2 = 1, nO; do i1 = 1, nO; do j3 = 1, nO; do k = 1, nO
        rdm3(i1, i2, i3, j1, :, j3) = rdm3(i1, i2, i3, j1, :, j3) + temp3(i1, i2, i3, j1, k, j3) * natorbs(k, :)
      end do; end do; end do; end do; end do; end do
    temp3 = 0.0_R64
    do i3 = 1, nO; do i2 = 1, nO; do i1 = 1, nO; do j2 = 1, nO; do j3 = 1, nO; do k = 1, nO
        temp3(i1, i2, i3, :, j2, j3) = temp3(i1, i2, i3, :, j2, j3) + rdm3(i1, i2, i3, k, j2, j3) * natorbs(k, :)
      end do; end do; end do; end do; end do; end do
    rdm3 = 0.0_R64
    do i2 = 1, nO; do i1 = 1, nO; do j1 = 1, nO; do j2 = 1, nO; do j3 = 1, nO; do k = 1, nO
        rdm3(i1, i2, :, j1, j2, j3) = rdm3(i1, i2, :, j1, j2, j3) + natorbs_H(:, k) * temp3(i1, i2, k, j1, j2, j3)
      end do; end do; end do; end do; end do; end do
    temp3 = 0.0_R64
    do i1 = 1, nO; do j1 = 1, nO; do j2 = 1, nO; do j3 = 1, nO; do i3 = 1, nO; do k = 1, nO
        temp3(i1, :, i3, j1, j2, j3) = temp3(i1, :, i3, j1, j2, j3) + natorbs_H(:, k) * rdm3(i1, k, i3, j1, j2, j3)
      end do; end do; end do; end do; end do; end do
    rdm3 = 0.0_R64
    do i2 = 1, nO; do j1 = 1, nO; do j2 = 1, nO; do j3 = 1, nO; do i3 = 1, nO; do k = 1, nO
        rdm3(:, i2, i3, j1, j2, j3) = rdm3(:, i2, i3, j1, j2, j3) + natorbs_H(:, k) * temp3(k, i2, i3, j1, j2, j3)
      end do; end do; end do; end do; end do; end do
    deallocate (temp3, natorbs_H)
  end subroutine

  subroutine RdmNatOrbTransform_InvTransformRdm3(rdm3, natorbs)
    complex(R64), intent(inout), contiguous :: rdm3(:, :, :, :, :, :)
    complex(R64), intent(in), contiguous :: natorbs(:, :)
    integer(I32) :: nO, i1, i2, i3, j1, j2, j3, k
    complex(R64), allocatable :: temp3(:, :, :, :, :, :), natorbs_H(:, :)

    nO = size(rdm3, 1)
    allocate (natorbs_H(nO, nO))
    natorbs_H = conjg(transpose(natorbs))

    allocate (temp3(nO, nO, nO, nO, nO, nO))
    ! Inverse transform: swap natorbs and natorbs_H
    ! Ket side (indices 4,5,6) multiplied by natorbs_H
    temp3 = 0.0_R64
    do j2 = 1, nO; do j1 = 1, nO; do i3 = 1, nO; do i2 = 1, nO; do i1 = 1, nO; do k = 1, nO
        temp3(i1, i2, i3, j1, j2, :) = temp3(i1, i2, i3, j1, j2, :) + rdm3(i1, i2, i3, j1, j2, k) * natorbs_H(k, :)
      end do; end do; end do; end do; end do; end do
    rdm3 = 0.0_R64
    do j1 = 1, nO; do i3 = 1, nO; do i2 = 1, nO; do i1 = 1, nO; do j3 = 1, nO; do k = 1, nO
        rdm3(i1, i2, i3, j1, :, j3) = rdm3(i1, i2, i3, j1, :, j3) + temp3(i1, i2, i3, j1, k, j3) * natorbs_H(k, :)
      end do; end do; end do; end do; end do; end do
    temp3 = 0.0_R64
    do i3 = 1, nO; do i2 = 1, nO; do i1 = 1, nO; do j2 = 1, nO; do j3 = 1, nO; do k = 1, nO
        temp3(i1, i2, i3, :, j2, j3) = temp3(i1, i2, i3, :, j2, j3) + rdm3(i1, i2, i3, k, j2, j3) * natorbs_H(k, :)
      end do; end do; end do; end do; end do; end do

    ! Bra side (indices 1,2,3) multiplied by natorbs
    rdm3 = 0.0_R64
    do i2 = 1, nO; do i1 = 1, nO; do j1 = 1, nO; do j2 = 1, nO; do j3 = 1, nO; do k = 1, nO
        rdm3(i1, i2, :, j1, j2, j3) = rdm3(i1, i2, :, j1, j2, j3) + natorbs(:, k) * temp3(i1, i2, k, j1, j2, j3)
      end do; end do; end do; end do; end do; end do
    temp3 = 0.0_R64
    do i1 = 1, nO; do j1 = 1, nO; do j2 = 1, nO; do j3 = 1, nO; do i3 = 1, nO; do k = 1, nO
        temp3(i1, :, i3, j1, j2, j3) = temp3(i1, :, i3, j1, j2, j3) + natorbs(:, k) * rdm3(i1, k, i3, j1, j2, j3)
      end do; end do; end do; end do; end do; end do
    rdm3 = 0.0_R64
    do i2 = 1, nO; do j1 = 1, nO; do j2 = 1, nO; do j3 = 1, nO; do i3 = 1, nO; do k = 1, nO
        rdm3(:, i2, i3, j1, j2, j3) = rdm3(:, i2, i3, j1, j2, j3) + natorbs(:, k) * temp3(k, i2, i3, j1, j2, j3)
      end do; end do; end do; end do; end do; end do
    deallocate (temp3, natorbs_H)
  end subroutine

end module
