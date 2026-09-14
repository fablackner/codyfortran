! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Trace and contraction utilities for reduced density matrices.
!>
!> Implements trace identities, consistency checks, and partial-trace
!> operations used to reduce higher-order RDMs and verify constraints.
module M_Utils_RdmTrace
  use M_Utils_Types

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmTrace_FillRdm1aFromRdm2(rdm1a, rdm2, nOa, nOb)

    complex(R64), intent(out), allocatable :: rdm1a(:, :)
    complex(R64), intent(in), contiguous                :: rdm2(:, :, :, :)
    integer(I32), intent(in) :: nOa, nOb

    complex(R64), allocatable :: rdm2ab(:, :, :, :)

    call RdmTrace_FillRdm2abFromRdm2(rdm2ab, rdm2, nOa, nOb)
    call RdmTrace_FillRdm1aFromRdm2ab(rdm1a, rdm2ab)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmTrace_FillRdm1aFromRdm2ab(rdm1a, rdm2ab)
    use M_Method_Mb

    complex(R64), intent(out), allocatable :: rdm1a(:, :)
    complex(R64), intent(in), contiguous                :: rdm2ab(:, :, :, :)

    integer(I32) :: nOa, nOb
    integer(I32) :: i, j, k

    nOa = size(rdm2ab, 1)
    nOb = size(rdm2ab, 3)

    if (.not. allocated(rdm1a)) allocate (rdm1a(nOa, nOa))

    rdm1a(:, :) = 0.0_R64
    do j = 1, nOa
      do i = 1, nOa

        do k = 1, nOb
          rdm1a(i, j) = rdm1a(i, j) + rdm2ab(i, k, j, k)
        end do

      end do
    end do

    rdm1a(:, :) = rdm1a(:, :) / Method_Mb_nBodies(2)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmTrace_FillRdm1bFromRdm2(rdm1b, rdm2, nOa, nOb)

    complex(R64), intent(out), allocatable :: rdm1b(:, :)
    complex(R64), intent(in), contiguous                :: rdm2(:, :, :, :)
    integer(I32), intent(in) :: nOa, nOb

    complex(R64), allocatable :: rdm2ab(:, :, :, :)

    call RdmTrace_FillRdm2abFromRdm2(rdm2ab, rdm2, nOa, nOb)
    call RdmTrace_FillRdm1bFromRdm2ab(rdm1b, rdm2ab)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmTrace_FillRdm1bFromRdm2ab(rdm1b, rdm2ab)
    use M_Method_Mb

    complex(R64), intent(out), allocatable :: rdm1b(:, :)
    complex(R64), intent(in), contiguous                :: rdm2ab(:, :, :, :)

    integer(I32) :: nOa, nOb
    integer(I32) :: i, j, k

    nOa = size(rdm2ab, 1)
    nOb = size(rdm2ab, 3)

    if (.not. allocated(rdm1b)) allocate (rdm1b(nOb, nOb))

    rdm1b(:, :) = 0.0_R64
    do j = 1, nOb
      do i = 1, nOb

        do k = 1, nOa
          rdm1b(i, j) = rdm1b(i, j) + rdm2ab(k, i, k, j)
        end do

      end do
    end do

    rdm1b(:, :) = rdm1b(:, :) / Method_Mb_nBodies(1)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmTrace_FillRdm1FromRdm2(rdm1, rdm2)
    use M_Method_Mb

    complex(R64), intent(out), allocatable :: rdm1(:, :)
    complex(R64), intent(in), contiguous                :: rdm2(:, :, :, :)

    integer(I32) :: nO
    integer(I32) :: i, j, k

    nO = size(rdm2, 1)

    if (.not. allocated(rdm1)) allocate (rdm1(nO, nO))

    rdm1(:, :) = 0.0_R64
    do j = 1, nO
      do i = 1, nO

        do k = 1, nO
          rdm1(i, j) = rdm1(i, j) + rdm2(i, k, j, k)
        end do

      end do
    end do

    rdm1(:, :) = rdm1(:, :) / (Method_Mb_nBodiesSum - 1)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmTrace_FillRdm1FromRdm2ab(rdm1, rdm2ab)

    complex(R64), intent(out), allocatable :: rdm1(:, :)
    complex(R64), intent(in), contiguous                :: rdm2ab(:, :, :, :)

    complex(R64), allocatable :: rdm1a(:, :)
    complex(R64), allocatable :: rdm1b(:, :)

    integer(I32) :: nOa, nOb, nO
    integer(I32) :: i, j

    nOa = size(rdm2ab, 1)
    nOb = size(rdm2ab, 3)
    nO = nOa + nOb

    call RdmTrace_FillRdm1aFromRdm2ab(rdm1a, rdm2ab)
    call RdmTrace_FillRdm1bFromRdm2ab(rdm1b, rdm2ab)

    if (.not. allocated(rdm1)) allocate (rdm1(nO, nO))

    rdm1(:, :) = 0.0_R64

    do j = 1, nOa
      do i = 1, nOa
        rdm1(i, j) = rdm1a(i, j)
      end do
    end do

    do j = 1, nOb
      do i = 1, nOb
        rdm1(i + nOa, j + nOa) = rdm1b(i, j)
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmTrace_FillRdm2abFromRdm2(rdm2ab, rdm2, nOa, nOb)

    complex(R64), intent(out), allocatable :: rdm2ab(:, :, :, :)
    complex(R64), intent(in), contiguous                :: rdm2(:, :, :, :)
    integer(I32), intent(in) :: nOa, nOb

    integer(I32) :: i1, i2, j1, j2

    if (.not. allocated(rdm2ab)) allocate (rdm2ab(nOa, nOb, nOa, nOb))

    do j2 = 1, nOb
      do j1 = 1, nOa
        do i2 = 1, nOb
          do i1 = 1, nOa

            rdm2ab(i1, i2, j1, j2) = rdm2(i1, i2 + nOa, j1, j2 + nOa)

          end do
        end do
      end do
    end do

  end subroutine

end module
