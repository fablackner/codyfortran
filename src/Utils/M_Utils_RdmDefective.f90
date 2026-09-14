! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Core routines for assembling and manipulating reduced density matrices.
!>
!> Provides 1/2/3-body RDM construction, projections, contractions, and
!> normalization/trace utilities used throughout the analysis pipeline.
module M_Utils_RdmDefective
  use M_Utils_Types

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmDefective_FillDefectiveRdm2sym(rdm2Defective, minEval, nFound, rdm2, threshhold)
    use M_Utils_RdmDiagonalize

    complex(R64), intent(out), allocatable :: rdm2Defective(:, :, :, :)
    real(R64), intent(out)                 :: minEval
    integer(I32), intent(out)              :: nFound
    complex(R64), intent(in), contiguous   :: rdm2(:, :, :, :)
    real(R64), intent(in)                  :: threshhold

    integer(I32) :: nO
    integer(I32) :: i1, i2, j1, j2, i

    real(R64), allocatable :: gemocc(:)
    complex(R64), allocatable :: gems(:, :, :)

    nO = size(rdm2, 1)

    if (.not. allocated(rdm2Defective)) allocate (rdm2Defective, mold=rdm2)
    rdm2Defective(:, :, :, :) = 0.0_R64

    call RdmDiagonalize_Rdm2sym(gemocc, gems, rdm2)

    nFound = 0

    do i = 1, size(gemocc)
      if (gemocc(i) > threshhold) cycle

      nFound = nFound + 1

      do j2 = 1, nO
        do j1 = 1, nO
          do i2 = 1, nO
            do i1 = 1, nO

              rdm2Defective(i1, i2, j1, j2) = rdm2Defective(i1, i2, j1, j2) + &
                                              gemocc(i) * gems(i1, i2, i) * conjg(gems(j1, j2, i))

            end do
          end do
        end do
      end do

    end do

    minEval = gemocc(size(gemocc))

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmDefective_FillDefectiveRdm2antisym(rdm2Defective, minEval, nFound, rdm2, threshhold)
    use M_Utils_RdmDiagonalize

    complex(R64), intent(out), allocatable :: rdm2Defective(:, :, :, :)
    real(R64), intent(out)                 :: minEval
    integer(I32), intent(out)              :: nFound
    complex(R64), intent(in), contiguous                :: rdm2(:, :, :, :)
    real(R64), intent(in)                  :: threshhold

    integer(I32) :: nO
    integer(I32) :: i1, i2, j1, j2, i

    real(R64), allocatable :: gemocc(:)
    complex(R64), allocatable :: gems(:, :, :)

    nO = size(rdm2, 1)

    if (.not. allocated(rdm2Defective)) allocate (rdm2Defective, mold=rdm2)
    rdm2Defective(:, :, :, :) = 0.0_R64

    call RdmDiagonalize_Rdm2antisym(gemocc, gems, rdm2)

    nFound = 0

    do i = 1, size(gemocc)
      if (gemocc(i) > threshhold) cycle

      nFound = nFound + 1

      do j2 = 1, nO
        do j1 = 1, nO
          do i2 = 1, nO
            do i1 = 1, nO

              rdm2Defective(i1, i2, j1, j2) = rdm2Defective(i1, i2, j1, j2) + &
                                              gemocc(i) * gems(i1, i2, i) * conjg(gems(j1, j2, i))

            end do
          end do
        end do
      end do

    end do

    minEval = gemocc(size(gemocc))

  end subroutine

end module
