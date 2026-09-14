! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Relations linking spin-singlet components of 2-body RDMs.
!>
!> Implements helpers to reconstruct same-spin blocks from opposite-spin
!> ones under singlet symmetry assumptions.
module M_Utils_RdmSingletRelations
  use M_Utils_Types

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmSingletRelations_FillRdm2aaFromRdm2ab(rdm2aa, rdm2ab)

    complex(R64), intent(out), allocatable :: rdm2aa(:, :, :, :)
    complex(R64), intent(in), contiguous                :: rdm2ab(:, :, :, :)

    integer(I32) :: nOa
    integer(I32) :: i1, i2, j1, j2

    nOa = size(rdm2ab, 1)

    if (.not. allocated(rdm2aa)) allocate (rdm2aa(nOa, nOa, nOa, nOa))

    do j2 = 1, nOa
      do j1 = 1, nOa
        do i2 = 1, nOa
          do i1 = 1, nOa

            rdm2aa(i1, i2, j1, j2) = rdm2ab(i1, i2, j1, j2) - rdm2ab(i2, i1, j1, j2)

          end do
        end do
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmSingletRelations_FillRdm2FromRdm2ab(rdm2, rdm2ab)

    complex(R64), intent(out), allocatable :: rdm2(:, :, :, :)
    complex(R64), intent(in), contiguous                :: rdm2ab(:, :, :, :)

    complex(R64) :: cTmp
    integer(I32) :: nOa, nOb, nO
    integer(I32) :: i1, i2, j1, j2

    nOa = size(rdm2ab, 1)
    nOb = size(rdm2ab, 2)
    nO = nOa + nOb

    if (.not. allocated(rdm2)) allocate (rdm2(nO, nO, nO, nO))

    rdm2(:, :, :, :) = 0.0_R64

    do j2 = 1, nOb
      do j1 = 1, nOa
        do i2 = 1, nOb
          do i1 = 1, nOa

            cTmp = rdm2ab(i1, i2, j1, j2)

            rdm2(i1, i2 + nOa, j1, j2 + nOa) = cTmp
            rdm2(i2 + nOa, i1, j1, j2 + nOa) = -cTmp
            rdm2(i1, i2 + nOa, j2 + nOa, j1) = -cTmp
            rdm2(i2 + nOa, i1, j2 + nOa, j1) = cTmp

          end do
        end do
      end do
    end do

    do j2 = 1, nOa
      do j1 = 1, nOa
        do i2 = 1, nOa
          do i1 = 1, nOa

            cTmp = rdm2ab(i1, i2, j1, j2) - rdm2ab(i2, i1, j1, j2)

            rdm2(i1, i2, j1, j2) = cTmp
            rdm2(i1 + nOa, i2 + nOa, j1 + nOa, j2 + nOa) = cTmp

          end do
        end do
      end do
    end do

  end subroutine

end module
