! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Trace and contraction utilities for reduced density matrices.
!>
!> Implements trace identities, consistency checks, and partial-trace
!> operations used to reduce higher-order RDMs and verify constraints.
module M_Utils_RdmHole
  use M_Utils_Types
  use M_Utils_RdmTrace

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmHole_FillQrdm2ab(qrdm2ab, rdm2ab)

    complex(R64), intent(out), allocatable :: qrdm2ab(:, :, :, :)
    complex(R64), intent(in), contiguous   :: rdm2ab(:, :, :, :)

    integer, allocatable :: delta_a(:, :), delta_b(:, :)
    complex(R64), allocatable :: rdm1a(:, :), rdm1b(:, :)
    integer(I32) :: i1, i2, j1, j2, nOa, nOb

    nOa = size(rdm2ab, 1)
    nOb = size(rdm2ab, 3)

    if (.not. allocated(qrdm2ab)) allocate (qrdm2ab(nOa, nOb, nOa, nOb))
    allocate (delta_a(nOa, nOa))
    allocate (delta_b(nOb, nOb))

    call RdmTrace_FillRdm1aFromRdm2ab(rdm1a, rdm2ab)
    call RdmTrace_FillRdm1bFromRdm2ab(rdm1b, rdm2ab)

    delta_a(:, :) = 0
    Do i1 = 1, nOa
      delta_a(i1, i1) = 1
    end do

    delta_b(:, :) = 0
    Do i2 = 1, nOb
      delta_b(i2, i2) = 1
    end do

    do j2 = 1, nOb
      do j1 = 1, nOa
        do i2 = 1, nOb
          do i1 = 1, nOa

            qrdm2ab(i1, i2, j1, j2) = rdm2ab(i1, i2, j1, j2) - delta_a(i1, j1) * rdm1b(i2, j2) - &
                                      rdm1a(i1, j1) * delta_b(i2, j2) + delta_a(i1, j1) * delta_b(i2, j2)

          end do
        end do
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmHole_FillGrdm2ab(grdm2ab, rdm2ab)

    complex(R64), intent(out), allocatable :: grdm2ab(:, :, :, :)
    complex(R64), intent(in), contiguous   :: rdm2ab(:, :, :, :)

    integer, allocatable :: delta_a(:, :)
    complex(R64), allocatable :: rdm1b(:, :)
    integer(I32) :: i1, i2, j1, j2, nOa, nOb

    nOa = size(rdm2ab, 1)
    nOb = size(rdm2ab, 3)

    if (.not. allocated(grdm2ab)) allocate (grdm2ab(nOa, nOb, nOa, nOb))
    allocate (delta_a(nOa, nOa))

    call RdmTrace_FillRdm1bFromRdm2ab(rdm1b, rdm2ab)

    delta_a(:, :) = 0
    Do i1 = 1, nOa
      delta_a(i1, i1) = 1
    end do

    do j2 = 1, nOb
      do j1 = 1, nOa
        do i2 = 1, nOb
          do i1 = 1, nOa
            ! G(i1, i2, j1, j2) = delta(i1, j1) * rdm1b(i2, j2) + rdm2ab(j1, i2, i1, j2)
            grdm2ab(i1, i2, j1, j2) = delta_a(i1, j1) * rdm1b(i2, j2) + rdm2ab(j1, i2, i1, j2)
          end do
        end do
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmHole_FillQrdm2(qrdm2, rdm2)

    complex(R64), intent(out), allocatable :: qrdm2(:, :, :, :)
    complex(R64), intent(in), contiguous   :: rdm2(:, :, :, :)

    integer, allocatable :: delta(:, :)
    complex(R64), allocatable :: rdm1(:, :)
    integer(I32) :: i1, i2, j1, j2, nO

    nO = size(rdm2, 1)

    if (.not. allocated(qrdm2)) allocate (qrdm2(nO, nO, nO, nO))
    allocate (delta(nO, nO))

    call RdmTrace_FillRdm1FromRdm2(rdm1, rdm2)

    delta(:, :) = 0
    Do i1 = 1, nO
      delta(i1, i1) = 1
    end do

    do j2 = 1, nO
      do j1 = 1, nO
        do i2 = 1, nO
          do i1 = 1, nO
            ! Q^{i1 i2}_{j1 j2} = (delta_{i1 j1}delta_{i2 j2} - delta_{i1 j2}delta_{i2 j1})
            !             - delta_{i2 j2}gamma^{j1}_{i1} - delta_{i1 j1}gamma^{j2}_{i2}
            !             + delta_{i2 j1}gamma^{j2}_{i1} + delta_{i1 j2}gamma^{j1}_{i2}
            !             + D^{j1 j2}_{i1 i2}

            qrdm2(i1, i2, j1, j2) = (delta(i1, j1) * delta(i2, j2) - delta(i1, j2) * delta(i2, j1)) &
                                    - delta(i2, j2) * rdm1(j1, i1) - delta(i1, j1) * rdm1(j2, i2) &
                                    + delta(i2, j1) * rdm1(j2, i1) + delta(i1, j2) * rdm1(j1, i2) &
                                    + rdm2(j1, j2, i1, i2)
          end do
        end do
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmHole_FillGrdm2(grdm2, rdm2)

    complex(R64), intent(out), allocatable :: grdm2(:, :, :, :)
    complex(R64), intent(in), contiguous   :: rdm2(:, :, :, :)

    integer, allocatable :: delta(:, :)
    complex(R64), allocatable :: rdm1(:, :)
    integer(I32) :: i1, i2, j1, j2, nO

    nO = size(rdm2, 1)

    if (.not. allocated(grdm2)) allocate (grdm2(nO, nO, nO, nO))
    allocate (delta(nO, nO))

    call RdmTrace_FillRdm1FromRdm2(rdm1, rdm2)

    delta(:, :) = 0
    Do i1 = 1, nO
      delta(i1, i1) = 1
    end do

    do j2 = 1, nO
      do j1 = 1, nO
        do i2 = 1, nO
          do i1 = 1, nO
            ! G^{i1 i2}_{j1 j2} = delta_{i2 j2} gamma^{i1}_{j1} + D^{i1 j2}_{j1 i2}
            grdm2(i1, i2, j1, j2) = delta(i2, j2) * rdm1(i1, j1) + rdm2(i1, j2, j1, i2)
          end do
        end do
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmHole_FillQrdm1a(qrdm1a, rdm1a)

    complex(R64), intent(out), allocatable :: qrdm1a(:, :)
    complex(R64), intent(in), contiguous   :: rdm1a(:, :)

    integer(I32) :: i, nOa

    nOa = size(rdm1a, 1)

    if (.not. allocated(qrdm1a)) allocate (qrdm1a(nOa, nOa))

    qrdm1a(:, :) = -rdm1a(:, :)
    do i = 1, nOa
      qrdm1a(i, i) = qrdm1a(i, i) + 1.0_R64
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmHole_FillQrdm1b(qrdm1b, rdm1b)

    complex(R64), intent(out), allocatable :: qrdm1b(:, :)
    complex(R64), intent(in), contiguous   :: rdm1b(:, :)

    integer(I32) :: i, nOb

    nOb = size(rdm1b, 1)

    if (.not. allocated(qrdm1b)) allocate (qrdm1b(nOb, nOb))

    qrdm1b(:, :) = -rdm1b(:, :)
    do i = 1, nOb
      qrdm1b(i, i) = qrdm1b(i, i) + 1.0_R64
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmHole_FillQrdm1(qrdm1, rdm1)

    complex(R64), intent(out), allocatable :: qrdm1(:, :)
    complex(R64), intent(in), contiguous   :: rdm1(:, :)

    integer(I32) :: i, nO

    nO = size(rdm1, 1)

    if (.not. allocated(qrdm1)) allocate (qrdm1(nO, nO))

    qrdm1(:, :) = -rdm1(:, :)
    do i = 1, nO
      qrdm1(i, i) = qrdm1(i, i) + 1.0_R64
    end do

  end subroutine

end module
