! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Cumulant-based decompositions and reconstructions of RDMs.
!>
!> Provides routines to separate connected and disconnected contributions of
!> higher-order RDMs and to build approximations via cumulant expansions.
module M_Utils_RdmCumulants
  use M_Utils_Types

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmCumulants_FillCumul2(cumul2, rdm1, rdm2)

    complex(R64), intent(out), allocatable :: cumul2(:, :, :, :)
    complex(R64), intent(in), contiguous                :: rdm1(:, :)
    complex(R64), intent(in), contiguous                :: rdm2(:, :, :, :)

    call RdmCumulants_FillCumulFreePartM2(cumul2, rdm1)

    cumul2(:, :, :, :) = rdm2(:, :, :, :) - cumul2(:, :, :, :)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmCumulants_FillCumul2ab(cumul2ab, rdm1a, rdm1b, rdm2ab)

    complex(R64), intent(out), allocatable :: cumul2ab(:, :, :, :)
    complex(R64), intent(in), contiguous                :: rdm1a(:, :)
    complex(R64), intent(in), contiguous                :: rdm1b(:, :)
    complex(R64), intent(in), contiguous                :: rdm2ab(:, :, :, :)

    integer(I32) :: i1, i2, j1, j2, nOa, nOb

    nOa = size(rdm1a, 1)
    nOb = size(rdm1b, 1)

    if (.not. allocated(cumul2ab)) allocate (cumul2ab(nOa, nOb, nOa, nOb))

    do j2 = 1, nOb
      do j1 = 1, nOa
        do i2 = 1, nOb
          do i1 = 1, nOa
            cumul2ab(i1, i2, j1, j2) = rdm2ab(i1, i2, j1, j2) - rdm1a(i1, j1) * rdm1b(i2, j2)
          end do
        end do
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmCumulants_FillCumul3(cumul3, rdm1, rdm2, rdm3)

    complex(R64), intent(out), allocatable :: cumul3(:, :, :, :, :, :)
    complex(R64), intent(in), contiguous                :: rdm1(:, :)
    complex(R64), intent(in), contiguous                :: rdm2(:, :, :, :)
    complex(R64), intent(in), contiguous                :: rdm3(:, :, :, :, :, :)

    integer(I32) :: nO

    nO = size(rdm2, 1)

    call RdmCumulants_FillCumulFreePartM3(cumul3, rdm1, rdm2)

    cumul3(:, :, :, :, :, :) = rdm3(:, :, :, :, :, :) - cumul3(:, :, :, :, :, :)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmCumulants_FillCumulFreePartM2(rdm2CumulFree, rdm1)

    complex(R64), intent(out), allocatable :: rdm2CumulFree(:, :, :, :)
    complex(R64), intent(in), contiguous                :: rdm1(:, :)

    integer(I32) :: i1, i2, j1, j2, nO

    nO = size(rdm1, 1)

    if (.not. allocated(rdm2CumulFree)) allocate (rdm2CumulFree(nO, nO, nO, nO))
    rdm2CumulFree(:, :, :, :) = 0.0_R64

    do j2 = 1, nO
      do j1 = 1, nO

        do i2 = 1, nO
          do i1 = 1, nO

            rdm2CumulFree(i1, i2, j1, j2) = rdm1(i1, j1) * rdm1(i2, j2) - rdm1(i1, j2) * rdm1(i2, j1)

          end do
        end do
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmCumulants_FillCumulFreePartM3(rdm3CumulFree, rdm1, rdm2)

    complex(R64), intent(out), allocatable :: rdm3CumulFree(:, :, :, :, :, :)
    complex(R64), intent(in), contiguous                :: rdm1(:, :)
    complex(R64), intent(in), contiguous                :: rdm2(:, :, :, :)

    complex(R64), allocatable :: cumul2(:, :, :, :)
    integer(I32) :: nO
    integer(I32) :: i1, i2, i3, j1, j2, j3

    nO = size(rdm2, 1)

    if (.not. allocated(rdm3CumulFree)) allocate (rdm3CumulFree(nO, nO, nO, nO, nO, nO))

    rdm3CumulFree(:, :, :, :, :, :) = 0.0_R64

    call RdmCumulants_FillCumul2(cumul2, rdm1, rdm2)

    ! 0th order

    do j3 = 1, nO
      do j2 = 1, nO
        do j1 = 1, nO

          do i3 = 1, nO
            do i2 = 1, nO
              do i1 = 1, nO

                rdm3CumulFree(i1, i2, i3, j1, j2, j3) = rdm3CumulFree(i1, i2, i3, j1, j2, j3) + &
                                                        rdm1(i1, j1) * rdm1(i2, j2) * rdm1(i3, j3) - &
                                                        rdm1(i1, j1) * rdm1(i2, j3) * rdm1(i3, j2) + &
                                                        rdm1(i1, j2) * rdm1(i2, j3) * rdm1(i3, j1) - &
                                                        rdm1(i1, j2) * rdm1(i2, j1) * rdm1(i3, j3) + &
                                                        rdm1(i1, j3) * rdm1(i2, j1) * rdm1(i3, j2) - &
                                                        rdm1(i1, j3) * rdm1(i2, j2) * rdm1(i3, j1)

              end do
            end do
          end do
        end do
      end do
    end do

    ! 1th order

    do j3 = 1, nO
      do j2 = 1, nO
        do j1 = 1, nO

          do i3 = 1, nO
            do i2 = 1, nO
              do i1 = 1, nO

                rdm3CumulFree(i1, i2, i3, j1, j2, j3) = rdm3CumulFree(i1, i2, i3, j1, j2, j3) + &
                                                        cumul2(i1, i2, j1, j2) * rdm1(i3, j3) - &
                                                        cumul2(i1, i3, j1, j2) * rdm1(i2, j3) - &
                                                        cumul2(i3, i2, j1, j2) * rdm1(i1, j3) - &
                                                        cumul2(i1, i2, j1, j3) * rdm1(i3, j2) + &
                                                        cumul2(i1, i3, j1, j3) * rdm1(i2, j2) - &
                                                        cumul2(i3, i2, j1, j3) * rdm1(i1, j2) - &
                                                        cumul2(i1, i2, j3, j2) * rdm1(i3, j1) + &
                                                        cumul2(i1, i3, j3, j2) * rdm1(i2, j1) + &
                                                        cumul2(i3, i2, j3, j2) * rdm1(i1, j1)

              end do
            end do
          end do
        end do
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmCumulants_FillCumulFreePartM3aab(rdm3aabCumulFree, rdm1a, rdm2ab)

    complex(R64), intent(out), allocatable :: rdm3aabCumulFree(:, :, :, :, :, :)
    complex(R64), intent(in), contiguous                :: rdm1a(:, :)
    complex(R64), intent(in), contiguous                :: rdm2ab(:, :, :, :)

    complex(R64), allocatable :: rdm2aa(:, :, :, :)
    integer(I32) :: nOa
    integer(I32) :: i1, i2, i3, j1, j2, j3

    nOa = size(rdm2ab, 1)

    if (.not. allocated(rdm3aabCumulFree)) allocate (rdm3aabCumulFree(nOa, nOa, nOa, nOa, nOa, nOa))

    rdm3aabCumulFree(:, :, :, :, :, :) = 0.0_R64

    allocate (rdm2aa(nOa, nOa, nOa, nOa))

    do j2 = 1, nOa
      do j1 = 1, nOa
        do i2 = 1, nOa
          do i1 = 1, nOa

            rdm2aa(i1, i2, j1, j2) = rdm2ab(i1, i2, j1, j2) - rdm2ab(i2, i1, j1, j2)

          end do
        end do
      end do
    end do

    ! the 3-RDM has to be initialized before calling the subroutine!
    ! this subroutine adds the valdemoro Reconstruction to the input 3-RDM

    do j2 = 1, nOa
      do j1 = 1, nOa
        if (j1 < j2) then

          do i2 = 1, nOa
            do i1 = 1, nOa
              if (i1 < i2) then

                do i3 = 1, nOa
                  do j3 = 1, nOa

                    rdm3aabCumulFree(i1, i2, i3, j1, j2, j3) = rdm3aabCumulFree(i1, i2, i3, j1, j2, j3) + &
                                                               rdm2aa(i1, i2, j1, j2) * rdm1a(i3, j3) + &
                                                               rdm2ab(i2, i3, j2, j3) * rdm1a(i1, j1) + &
                                                               rdm2ab(i1, i3, j1, j3) * rdm1a(i2, j2) - &
                                                               rdm2ab(i1, i3, j2, j3) * rdm1a(i2, j1) - &
                                                               rdm2ab(i2, i3, j1, j3) * rdm1a(i1, j2) - &
                                                               2.0_R64 * ( &
                                                               rdm1a(i1, j1) * rdm1a(i2, j2) - &
                                                               rdm1a(i2, j1) * rdm1a(i1, j2) &
                                                               ) * &
                                                               rdm1a(i3, j3)

                    rdm3aabCumulFree(i2, i1, i3, j1, j2, j3) = -rdm3aabCumulFree(i1, i2, i3, j1, j2, j3)
                    rdm3aabCumulFree(i1, i2, i3, j2, j1, j3) = -rdm3aabCumulFree(i1, i2, i3, j1, j2, j3)
                    rdm3aabCumulFree(i2, i1, i3, j2, j1, j3) = rdm3aabCumulFree(i1, i2, i3, j1, j2, j3)

                  end do
                end do
              end if
            end do
          end do
        end if
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmCumulants_AddNyCumulantM3aab(rdm3aab, nBodiesUP, rdm1a, rdm2ab)
    use M_Utils_RdmDiagonalize

    complex(R64), intent(inout), contiguous  :: rdm3aab(:, :, :, :, :, :)
    integer(I32), intent(in)    :: nBodiesUP
    complex(R64), intent(in), contiguous     :: rdm1a(:, :)
    complex(R64), intent(in), contiguous     :: rdm2ab(:, :, :, :)

    complex(R64), allocatable :: cumul2ab(:, :, :, :)
    complex(R64), allocatable :: cumul2aa(:, :, :, :)
    complex(R64), allocatable :: Rxcumul2ab(:, :, :, :)
    complex(R64), allocatable :: Rxcumul2aa(:, :, :, :)
    complex(R64), allocatable :: R(:, :)
    real(R64), allocatable    :: natocc(:)
    complex(R64), allocatable :: natorbs(:, :)

    integer(I32) :: nOa
    integer(I32) :: i1, i2, i3, j1, j2, j3, i, k

    nOa = size(rdm1a, 1)

    allocate (cumul2ab(nOa, nOa, nOa, nOa))
    allocate (cumul2aa(nOa, nOa, nOa, nOa))
    allocate (Rxcumul2ab(nOa, nOa, nOa, nOa))
    allocate (Rxcumul2aa(nOa, nOa, nOa, nOa))
    allocate (R(nOa, nOa))

    do j2 = 1, nOa
      do j1 = 1, nOa
        do i2 = 1, nOa
          do i1 = 1, nOa

            cumul2ab(i1, i2, j1, j2) = rdm2ab(i1, i2, j1, j2) - rdm1a(i1, j1) * rdm1a(i2, j2)

          end do
        end do
      end do
    end do

    do j2 = 1, nOa
      do j1 = 1, nOa
        do i2 = 1, nOa
          do i1 = 1, nOa

            cumul2aa(i1, i2, j1, j2) = cumul2ab(i1, i2, j1, j2) - cumul2ab(i2, i1, j1, j2)

          end do
        end do
      end do
    end do

    call RdmDiagonalize_Rdm1(natocc, natorbs, rdm1a)

    do i = 1, nBodiesUP
      natocc(i) = 1.0_R64
    end do

    do i = nBodiesUP + 1, nOa
      natocc(i) = -1.0_R64
    end do

    R = 0.0_R64
    do j1 = 1, nOa
      do i1 = 1, nOa
        do j2 = 1, nOa

          R(i1, j1) = R(i1, j1) + natorbs(i1, j2) * natocc(j2) * conjg(natorbs(j1, j2))

        end do
      end do
    end do

    Rxcumul2ab(:, :, :, :) = 0.0_R64
    Rxcumul2aa(:, :, :, :) = 0.0_R64

    do j2 = 1, nOa
      do j1 = 1, nOa
        do i2 = 1, nOa
          do i1 = 1, nOa

            do k = 1, nOa

              Rxcumul2ab(i1, i2, j1, j2) = Rxcumul2ab(i1, i2, j1, j2) + R(i1, k) * cumul2ab(k, i2, j1, j2)
              Rxcumul2aa(i1, i2, j1, j2) = Rxcumul2aa(i1, i2, j1, j2) + R(i1, k) * cumul2aa(k, i2, j1, j2)

            end do
          end do
        end do
      end do
    end do

    ! the 3-RDM has to be initialized before calling the subroutine!
    !this subroutine adds the valdemoro Reconstruction to the input 3-RDM

    do i2 = 1, nOa
      do i1 = 1, nOa
        if (i1 < i2) then

          do j2 = 1, nOa
            do j1 = 1, nOa
              if (j1 < j2) then

                do j3 = 1, nOa
                  do i3 = 1, nOa

                    do k = 1, nOa

                      rdm3aab(i1, i2, i3, j1, j2, j3) = rdm3aab(i1, i2, i3, j1, j2, j3) + &
                                                        cumul2aa(i1, i2, j1, k) * Rxcumul2ab(k, i3, j2, j3) - &
                                                        cumul2aa(i1, i2, j2, k) * Rxcumul2ab(k, i3, j1, j3) + &
                                                        cumul2ab(i1, i3, j1, k) * Rxcumul2ab(k, i2, j3, j2) - &
                                                        cumul2ab(i1, i3, j2, k) * Rxcumul2ab(k, i2, j3, j1) + &
                                                        cumul2ab(i1, i3, k, j3) * Rxcumul2aa(k, i2, j1, j2) - &
                                                        cumul2ab(i2, i3, j1, k) * Rxcumul2ab(k, i1, j3, j2) + &
                                                        cumul2ab(i2, i3, j2, k) * Rxcumul2ab(k, i1, j3, j1) - &
                                                        cumul2ab(i2, i3, k, j3) * Rxcumul2aa(k, i1, j1, j2)

                    end do

                    rdm3aab(i2, i1, i3, j1, j2, j3) = -rdm3aab(i1, i2, i3, j1, j2, j3)
                    rdm3aab(i1, i2, i3, j2, j1, j3) = -rdm3aab(i1, i2, i3, j1, j2, j3)
                    rdm3aab(i2, i1, i3, j2, j1, j3) = rdm3aab(i1, i2, i3, j1, j2, j3)

                  end do
                end do
              end if
            end do
          end do
        end if
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end module
