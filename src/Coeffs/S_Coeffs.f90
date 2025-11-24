! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Coeffs) S_Coeffs

  implicit none

!=============================================================================
! local procedures kinds
!=============================================================================

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Coeffs_Fabricate
    use M_Utils_Say
    use M_Utils_Json
    use M_Coeffs_Generic
    use M_Coeffs_Hubbard

    call Say_Fabricate("coeffs")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    Coeffs_Normalize => Normalize
    Coeffs_ProjectOnSubspace => ProjectOnSubspace
    Coeffs_SaveCoeffs => SaveCoeffs
    Coeffs_SaveTwoRdm => SaveTwoRdm
    Coeffs_FillRdm3Bt => FillRdm3Bt
    Coeffs_FillRdm2Bt => FillRdm2Bt
    Coeffs_FillRdm1Bt => FillRdm1Bt

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("coeffs.generic")) then
      call Coeffs_Generic_Fabricate

    else if (Json_GetExistence("coeffs.hubbard")) then
      call Coeffs_Hubbard_Fabricate

    else
      error stop "coeffs is missing one of: generic, hubbard"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Normalize(coeffs)
    use M_Utils_BlasLib

    complex(R64), intent(inout), contiguous  :: coeffs(:)

    real(R64) :: norm

    norm = BlasLib_CalcNorm(coeffs)
    coeffs(:) = coeffs(:) / norm

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ProjectOnSubspace(dCoeffs, coeffs)
    complex(R64), intent(inout), contiguous  :: dCoeffs(:)
    complex(R64), intent(in), contiguous     :: coeffs(:)

    dCoeffs(:) = dCoeffs(:) - dot_product(coeffs, dCoeffs) * coeffs(:)

  end subroutine

!------------------------------------------------------------------------------
! Helper: group excitations by body type and apply per body type in one call
!------------------------------------------------------------------------------
  subroutine ApplyGroupedExcitations(coeffsTmp, iExc, jExc, btExc, nExc)
    complex(R64), intent(inout), contiguous :: coeffsTmp(:)
    integer(I32), intent(in) :: iExc(:)
    integer(I32), intent(in) :: jExc(:)
    integer(I32), intent(in) :: btExc(:)
    integer, intent(in) :: nExc

    integer(I32) :: uniqBts(nExc)
    integer(I32) :: iBuf(nExc), jBuf(nExc)
    integer      :: nGroups, g, k, cnt

    nGroups = 0
    do k = 1, nExc
      if (.not. any(uniqBts(1:nGroups) .eq. btExc(k))) then
        nGroups = nGroups + 1
        uniqBts(nGroups) = btExc(k)
      end if
    end do

    do g = 1, nGroups
      cnt = 0
      do k = 1, nExc
        if (btExc(k) .eq. uniqBts(g)) then
          cnt = cnt + 1
          iBuf(cnt) = iExc(k)
          jBuf(cnt) = jExc(k)
        end if
      end do
      call Coeffs_ApplyExcitation(coeffsTmp, iBuf(1:cnt), jBuf(1:cnt), uniqBts(g))
    end do
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine SaveCoeffs(coeffs)
    use M_Utils_DataStorage

    complex(R64), intent(in), contiguous  :: coeffs(:)

    call SaveData(coeffs, 'coeffs.dat', storage_size(coeffs))

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine SaveTwoRdm(coeffs)
    use M_Utils_DataStorage
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    complex(R64), intent(in), contiguous  :: coeffs(:)

    complex(R64), allocatable :: rdm2(:, :, :, :)
    complex(R64), allocatable :: rdm2Bt(:, :, :, :)
    character(len=256) :: filename
    integer(I32)      :: ibt1, ibt2

    call Coeffs_ApplyH2FillRdm2(coeffs, rdm2_=rdm2)

    do ibt1 = 1, Method_Mb_nBodyTypes
      do ibt2 = 1, Method_Mb_nBodyTypes

        write (filename, '("rdm2Bt",I2.2,"_",I2.2,".dat")') ibt1, ibt2

        rdm2Bt = rdm2(Method_Mb_OrbBased_nOrbsStart(ibt1):Method_Mb_OrbBased_nOrbsEnd(ibt1), &
                      Method_Mb_OrbBased_nOrbsStart(ibt2):Method_Mb_OrbBased_nOrbsEnd(ibt2), &
                      Method_Mb_OrbBased_nOrbsStart(ibt1):Method_Mb_OrbBased_nOrbsEnd(ibt1), &
                      Method_Mb_OrbBased_nOrbsStart(ibt2):Method_Mb_OrbBased_nOrbsEnd(ibt2))

        call SaveData(rdm2Bt, trim(filename), storage_size(rdm2Bt))

      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillRdm3Bt(rdm3Bt, coeffs, bt1, bt2, bt3)
    use M_Method_Mb_OrbBased

    complex(R64), intent(out), allocatable :: rdm3Bt(:, :, :, :, :, :)
    complex(R64), intent(in), contiguous, target     :: coeffs(:)
    integer(I32), intent(in) :: bt1
    integer(I32), intent(in) :: bt2
    integer(I32), intent(in) :: bt3

    complex(R64), allocatable  :: coeffsTmp(:)

    integer(I32) :: i1, i2, i3, j1, j2, j3

    if (.not. allocated(rdm3Bt)) allocate (rdm3Bt(Method_Mb_OrbBased_nOrbs(bt1), &
                                                  Method_Mb_OrbBased_nOrbs(bt2), &
                                                  Method_Mb_OrbBased_nOrbs(bt3), &
                                                  Method_Mb_OrbBased_nOrbs(bt1), &
                                                  Method_Mb_OrbBased_nOrbs(bt2), &
                                                  Method_Mb_OrbBased_nOrbs(bt3)))

    allocate (coeffsTmp, source=coeffs)

    do j3 = 1, size(rdm3Bt, 6)
      do j2 = 1, size(rdm3Bt, 5)
        do j1 = 1, size(rdm3Bt, 4)
          do i3 = 1, size(rdm3Bt, 3)
            do i2 = 1, size(rdm3Bt, 2)
              do i1 = 1, size(rdm3Bt, 1)

                coeffsTmp = coeffs

                call ApplyGroupedExcitations(coeffsTmp, [i1, i2, i3], [j1, j2, j3], [bt1, bt2, bt3], 3)

                rdm3Bt(i1, i2, i3, j1, j2, j3) = dot_product(coeffs, coeffsTmp)
              end do
            end do
          end do
        end do
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillRdm2Bt(rdm2Bt, coeffs, bt1, bt2)
    use M_Method_Mb_OrbBased

    complex(R64), intent(out), allocatable :: rdm2Bt(:, :, :, :)
    complex(R64), intent(in), contiguous, target     :: coeffs(:)
    integer(I32), intent(in) :: bt1
    integer(I32), intent(in) :: bt2

    complex(R64), allocatable  :: coeffsTmp(:)

    integer(I32) :: i1, i2, j1, j2

    if (.not. allocated(rdm2Bt)) allocate (rdm2Bt(Method_Mb_OrbBased_nOrbs(bt1), &
                                                  Method_Mb_OrbBased_nOrbs(bt2), &
                                                  Method_Mb_OrbBased_nOrbs(bt1), &
                                                  Method_Mb_OrbBased_nOrbs(bt2)))

    allocate (coeffsTmp, source=coeffs)

    do j2 = 1, size(rdm2Bt, 4)
      do j1 = 1, size(rdm2Bt, 3)
        do i2 = 1, size(rdm2Bt, 2)
          do i1 = 1, size(rdm2Bt, 1)

            coeffsTmp = coeffs

            call ApplyGroupedExcitations(coeffsTmp, [i1, i2], [j1, j2], [bt1, bt2], 2)

            rdm2Bt(i1, i2, j1, j2) = dot_product(coeffs, coeffsTmp)
          end do
        end do
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillRdm1Bt(rdm1Bt, coeffs, bt1)
    use M_Method_Mb_OrbBased

    complex(R64), intent(out), allocatable :: rdm1Bt(:, :)
    complex(R64), intent(in), contiguous, target     :: coeffs(:)
    integer(I32), intent(in) :: bt1

    complex(R64), allocatable  :: coeffsTmp(:)

    integer(I32) :: i1, j1

    if (.not. allocated(rdm1Bt)) allocate (rdm1Bt(Method_Mb_OrbBased_nOrbs(bt1), &
                                                  Method_Mb_OrbBased_nOrbs(bt1)))

    allocate (coeffsTmp, source=coeffs)

    do j1 = 1, size(rdm1Bt, 2)
      do i1 = 1, size(rdm1Bt, 1)

        coeffsTmp = coeffs

        call ApplyGroupedExcitations(coeffsTmp, [i1], [j1], [bt1], 1)

        rdm1Bt(i1, j1) = dot_product(coeffs, coeffsTmp)

      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end submodule
