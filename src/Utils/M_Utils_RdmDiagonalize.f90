! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Core routines for assembling and manipulating reduced density matrices.
!>
!> Provides 1/2/3-body RDM construction, projections, contractions, and
!> normalization/trace utilities used throughout the analysis pipeline.
module M_Utils_RdmDiagonalize
  use M_Utils_Types

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmDiagonalize_Rdm1(natocc, natorbs, rdm1)
    use M_Utils_LapackLib

    real(R64), intent(out), allocatable    :: natocc(:)
    complex(R64), intent(out), allocatable :: natorbs(:, :)
    complex(R64), intent(in), contiguous                :: rdm1(:, :)

    integer(I32) :: nO

    nO = size(rdm1, 1)

    if (.not. allocated(natocc)) allocate (natocc(nO))
    if (.not. allocated(natorbs)) allocate (natorbs(1:nO, 1:nO))

    call LapackLib_DiagonalizeGeneric(natocc, natorbs, rdm1, .true.)

    natocc(:) = natocc(nO:1:-1) ! invert the vals into descending order
    natorbs(:, :) = natorbs(:, nO:1:-1) ! invert the vecs into descending order

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmDiagonalize_Rdm2(gemocc, gems, rdm2)
    use M_Utils_LapackLib

    real(R64), intent(out), allocatable          :: gemocc(:)
    complex(R64), intent(out), allocatable       :: gems(:, :, :)
    complex(R64), intent(in), contiguous, target :: rdm2(:, :, :, :)

    complex(R64), contiguous, pointer :: rdm2Mat(:, :)

    complex(R64), allocatable, target :: gemsMat(:, :)

    integer(I32) :: iGem1, nO, nDim
    integer(I32) :: iOrb1, iOrb2
    integer(I32) :: count1

    nO = size(rdm2, 1)
    nDim = nO**2

    rdm2Mat(1:nDim, 1:nDim) => rdm2(:, :, :, :)

    allocate (gemsMat(1:nDim, 1:nDim))

    if (.not. allocated(gemocc)) allocate (gemocc(nDim))
    if (.not. allocated(gems)) allocate (gems(nO, nO, nDim))

    call LapackLib_DiagonalizeGeneric(gemocc, gemsMat, rdm2Mat, .true.)

    gemocc(:) = gemocc(nDim:1:-1) ! invert the array into descending order

    do iGem1 = 1, nDim

      count1 = 1
      do iOrb2 = 1, nO
        do iOrb1 = 1, nO

          gems(iOrb1, iOrb2, iGem1) = gemsMat(count1, nDim - iGem1 + 1)
          count1 = count1 + 1
        end do
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmDiagonalize_Rdm2sym(gemocc, gems, rdm2)
    use M_Utils_LapackLib

    real(R64), intent(out), allocatable    :: gemocc(:)
    complex(R64), intent(out), allocatable :: gems(:, :, :)
    complex(R64), intent(in), contiguous                :: rdm2(:, :, :, :)

    complex(R64), allocatable :: rdm2Mat(:, :)

    complex(R64), allocatable :: gemsMat(:, :)

    integer(I32) :: nO, nDim
    integer(I32) :: i1, i2, j1, j2, iGem1
    integer(I32) :: count1, count2

    nO = size(rdm2, 1)
    nDim = nO * (nO + 1) / 2

    allocate (rdm2Mat(nDim, nDim))
    rdm2Mat(:, :) = 0.0_R64

    count2 = 0
    do j2 = 1, nO
      do j1 = 1, nO
        if (j1 <= j2) then
          count2 = count2 + 1
          count1 = 0
          do i2 = 1, nO
            do i1 = 1, nO
              if (i1 <= i2) then
                count1 = count1 + 1

                rdm2Mat(count1, count2) = 2.0_R64 * rdm2(i1, i2, j1, j2)

                if (i1 .eq. i2) rdm2Mat(count1, count2) = sqrt(0.5_R64) * rdm2Mat(count1, count2)
                if (j1 .eq. j2) rdm2Mat(count1, count2) = sqrt(0.5_R64) * rdm2Mat(count1, count2)

              end if
            end do
          end do
        end if
      end do
    end do

    allocate (gemsMat(1:nDim, 1:nDim))

    if (.not. allocated(gemocc)) allocate (gemocc(nDim))
    if (.not. allocated(gems)) allocate (gems(nO, nO, nDim))

    call LapackLib_DiagonalizeGeneric(gemocc, gemsMat, rdm2Mat, .true.)

    gemocc(:) = gemocc(nDim:1:-1) ! invert the array into descending order

    gems(:, :, :) = 0.0_R64

    do iGem1 = 1, nDim

      count1 = 0
      do i2 = 1, nO
        do i1 = 1, nO
          if (i1 < i2) then
            count1 = count1 + 1

            gems(i1, i2, iGem1) = sqrt(0.5_R64) * gemsMat(count1, nDim - iGem1 + 1)
            gems(i2, i1, iGem1) = sqrt(0.5_R64) * gemsMat(count1, nDim - iGem1 + 1)

          else if (i1 .eq. i2) then
            count1 = count1 + 1

            gems(i1, i1, iGem1) = gemsMat(count1, nDim - iGem1 + 1)

          end if
        end do
      end do

    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmDiagonalize_Rdm2antisym(gemocc, gems, rdm2)
    use M_Utils_LapackLib

    real(R64), intent(out), allocatable    :: gemocc(:)
    complex(R64), intent(out), allocatable :: gems(:, :, :)
    complex(R64), intent(in), contiguous                :: rdm2(:, :, :, :)

    complex(R64), allocatable :: rdm2Mat(:, :)

    complex(R64), allocatable, target :: gemsMat(:, :)

    integer(I32) :: nO, nDim
    integer(I32) :: i1, i2, j1, j2, iGem1
    integer(I32) :: count1, count2

    nO = size(rdm2, 1)
    nDim = nO * (nO - 1) / 2

    allocate (rdm2Mat(nDim, nDim))
    rdm2Mat(:, :) = 0.0_R64

    count2 = 0
    do j2 = 1, nO
      do j1 = 1, nO
        if (j1 < j2) then
          count2 = count2 + 1
          count1 = 0
          do i2 = 1, nO
            do i1 = 1, nO
              if (i1 < i2) then
                count1 = count1 + 1

                rdm2Mat(count1, count2) = 2.0_R64 * rdm2(i1, i2, j1, j2)

              end if
            end do
          end do
        end if
      end do
    end do

    allocate (gemsMat(1:nDim, 1:nDim))

    if (.not. allocated(gemocc)) allocate (gemocc(nDim))
    if (.not. allocated(gems)) allocate (gems(nO, nO, nDim))

    call LapackLib_DiagonalizeGeneric(gemocc, gemsMat, rdm2Mat, .true.)

    gemocc(:) = gemocc(nDim:1:-1) ! invert the array into descending order

    gems(:, :, :) = 0.0_R64

    do iGem1 = 1, nDim

      count1 = 0
      do i2 = 1, nO
        do i1 = 1, nO
          if (i1 < i2) then
            count1 = count1 + 1

            gems(i1, i2, iGem1) = sqrt(0.5_R64) * gemsMat(count1, nDim - iGem1 + 1)
            gems(i2, i1, iGem1) = -sqrt(0.5_R64) * gemsMat(count1, nDim - iGem1 + 1)

          end if
        end do
      end do
    end do

  end subroutine

end module
