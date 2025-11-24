! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Core routines for assembling and manipulating reduced density matrices.
!>
!> Provides 1/2/3-body RDM construction, projections, contractions, and
!> normalization/trace utilities used throughout the analysis pipeline.
module M_Utils_RdmBasics
  use M_Utils_Types

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function RdmBasics_CalculateEnergy(rdm1, h1, rdm2, h2) result(res)

    real(R64)                :: res
    complex(R64), intent(in), contiguous  :: rdm1(:, :)
    complex(R64), intent(in), contiguous  :: h1(:, :)
    complex(R64), intent(in), contiguous  :: rdm2(:, :, :, :)
    complex(R64), intent(in), contiguous  :: h2(:, :, :, :)

    complex(R64) :: cTmp

    integer(I32) :: k1, l1, k2, l2, nO

    nO = size(h1, 1)

    cTmp = 0.0_R64
    do k1 = 1, nO
      do l1 = 1, nO
        cTmp = cTmp + h1(k1, l1) * rdm1(l1, k1)
      end do
    end do

    do k1 = 1, nO
      do k2 = 1, nO
        do l1 = 1, nO
          do l2 = 1, nO

            cTmp = cTmp + h2(k1, k2, l1, l2) * rdm2(l1, l2, k1, k2)

          end do
        end do
      end do
    end do

    res = real(cTmp, kind=R64)

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmBasics_Antisymmetrize(A2)

    complex(R64), intent(inout), contiguous :: A2(:, :, :, :)

    integer(I32) :: nO
    integer(I32) :: i1, i2, j1, j2
    complex(R64), allocatable :: A2Tmp(:, :, :, :)

    nO = size(A2, 1)
    allocate (A2Tmp, source=A2)

    do j2 = 1, nO
      do j1 = 1, nO

        do i2 = 1, nO
          do i1 = 1, nO

            A2(i1, i2, j1, j2) = 0.25_R64 * (A2Tmp(i1, i2, j1, j2) - &
                                             A2Tmp(i2, i1, j1, j2) - &
                                             A2Tmp(i1, i2, j2, j1) + &
                                             A2Tmp(i2, i1, j2, j1))

          end do
        end do
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmBasics_Symmetrize(A2)

    complex(R64), intent(inout), contiguous :: A2(:, :, :, :)

    integer(I32) :: nO
    integer(I32) :: i1, i2, j1, j2
    complex(R64), allocatable :: A2Tmp(:, :, :, :)

    nO = size(A2, 1)
    allocate (A2Tmp, source=A2)

    do j2 = 1, nO
      do j1 = 1, nO

        do i2 = 1, nO
          do i1 = 1, nO

            A2(i1, i2, j1, j2) = 0.25_R64 * (A2Tmp(i1, i2, j1, j2) + &
                                             A2Tmp(i2, i1, j1, j2) + &
                                             A2Tmp(i1, i2, j2, j1) + &
                                             A2Tmp(i2, i1, j2, j1))

          end do
        end do
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmBasics_DiagonalizeRdm1(natocc, natorbs, rdm1)
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
  subroutine RdmBasics_DiagonalizeRdm2(gemocc, gems, rdm2)
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
  subroutine RdmBasics_DiagonalizeRdm2sym(gemocc, gems, rdm2)
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
  subroutine RdmBasics_DiagonalizeRdm2antisym(gemocc, gems, rdm2)
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

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmBasics_FillDefectiveRdm2sym(rdm2Defective, minEval, nFound, rdm2, threshhold)

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

    call RdmBasics_DiagonalizeRdm2sym(gemocc, gems, rdm2)

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
  subroutine RdmBasics_FillDefectiveRdm2antisym(rdm2Defective, minEval, nFound, rdm2, threshhold)

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

    call RdmBasics_DiagonalizeRdm2antisym(gemocc, gems, rdm2)

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
