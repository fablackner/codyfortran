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
  !> Free-cumulant (Valdemoro) part of the up-up-down 3-RDM from rdm1a and a
  !> singlet rdm2ab, ket-first: every element (antisymmetric in (i1,i2), and in
  !> (j1,j2) since rdm2ab(i1,i2,j1,j2) = rdm2ab(i2,i1,j2,j1)); i1 innermost.
  subroutine RdmCumulants_FillCumulFreePartM3aab(rdm3aabCumulFree, rdm1a, rdm2ab)
    use M_Utils_RdmSingletRelations

    complex(R64), intent(out), allocatable :: rdm3aabCumulFree(:, :, :, :, :, :)
    complex(R64), intent(in), contiguous   :: rdm1a(:, :)
    complex(R64), intent(in), contiguous   :: rdm2ab(:, :, :, :)

    complex(R64), allocatable :: rdm2aa(:, :, :, :)
    integer(I32) :: nOa, i1, i2, i3, j1, j2, j3

    nOa = size(rdm2ab, 1)
    allocate (rdm3aabCumulFree(nOa, nOa, nOa, nOa, nOa, nOa))
    call RdmSingletRelations_FillRdm2aaFromRdm2ab(rdm2aa, rdm2ab)

    !$omp parallel do collapse(2) private(j3, j2, j1, i3, i2, i1)
    do j3 = 1, nOa
      do j2 = 1, nOa
        do j1 = 1, nOa
          do i3 = 1, nOa
            do i2 = 1, nOa
              do i1 = 1, nOa

                rdm3aabCumulFree(i1, i2, i3, j1, j2, j3) = rdm2aa(i1, i2, j1, j2) * rdm1a(i3, j3) + &
                                                           rdm2ab(i2, i3, j2, j3) * rdm1a(i1, j1) + &
                                                           rdm2ab(i1, i3, j1, j3) * rdm1a(i2, j2) - &
                                                           rdm2ab(i1, i3, j2, j3) * rdm1a(i2, j1) - &
                                                           rdm2ab(i2, i3, j1, j3) * rdm1a(i1, j2) - &
                                                           2.0_R64 * ( &
                                                           rdm1a(i1, j1) * rdm1a(i2, j2) - &
                                                           rdm1a(i2, j1) * rdm1a(i1, j2) &
                                                           ) * &
                                                           rdm1a(i3, j3)

              end do
            end do
          end do
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Adds the Nakatsuji-Yasuda 3-cumulant of a singlet rdm2ab to the
  !> up-up-down 3-RDM (ket-first): products of two 2-cumulants contracted
  !> through R = 2 P_occ - 1, P_occ the projector on the nBodiesUp most
  !> occupied natural orbitals of rdm1a.
  subroutine RdmCumulants_AddNyCumulantM3aab(rdm3aab, nBodiesUp, rdm1a, rdm2ab)
    use M_Utils_RdmDiagonalize
    use M_Utils_RdmSingletRelations
    use M_Utils_BlasLib

    complex(R64), intent(inout), contiguous :: rdm3aab(:, :, :, :, :, :)
    integer(I32), intent(in)                :: nBodiesUp
    complex(R64), intent(in), contiguous    :: rdm1a(:, :)
    complex(R64), intent(in), contiguous    :: rdm2ab(:, :, :, :)

    complex(R64), allocatable :: cumul2ab(:, :, :, :), cumul2aa(:, :, :, :), cumul2abPerm(:, :, :, :)
    complex(R64), allocatable :: Rxcumul2ab(:, :, :, :), Rxcumul2aa(:, :, :, :), R(:, :), natorbs(:, :)
    complex(R64), allocatable :: p1t(:, :, :, :, :), p3t(:, :, :, :, :), p5t(:, :, :, :, :), b3(:, :, :)
    real(R64), allocatable    :: natocc(:)
    integer(I32) :: nOa, i1, i2, i3, j1, j2, j3, k

    nOa = size(rdm1a, 1)

    call RdmCumulants_FillCumul2ab(cumul2ab, rdm1a, rdm1a, rdm2ab)
    call RdmSingletRelations_FillRdm2aaFromRdm2ab(cumul2aa, cumul2ab)

    call RdmDiagonalize_Rdm1(natocc, natorbs, rdm1a)
    natocc(:nBodiesUp) = 1.0_R64
    natocc(nBodiesUp + 1:) = -1.0_R64
    R = matmul(natorbs * spread(natocc, 1, nOa), conjg(transpose(natorbs)))

    ! R cumul2: one matrix product each, R (n x n) times cumul2 as (n x n^3)
    allocate (Rxcumul2ab, mold=cumul2ab)
    allocate (Rxcumul2aa, mold=cumul2ab)
    call ZGEMM('N', 'N', nOa, nOa**3, nOa, (1.0_R64, 0.0_R64), R, nOa, cumul2ab, nOa, (0.0_R64, 0.0_R64), Rxcumul2ab, nOa)
    call ZGEMM('N', 'N', nOa, nOa**3, nOa, (1.0_R64, 0.0_R64), R, nOa, cumul2aa, nOa, (0.0_R64, 0.0_R64), Rxcumul2aa, nOa)

    ! The eight products are permutations of three sums over k,
    !   p1(i1,i2,j1, i3,j2,j3) = sum_k cumul2aa(i1,i2,j1,k) Rxcumul2ab(k,i3,j2,j3)
    !   p3(i1,i3,j1, i2,j3,j2) = sum_k cumul2ab(i1,i3,j1,k) Rxcumul2ab(k,i2,j3,j2)
    !   p5(i1,i3,j3, i2,j1,j2) = sum_k cumul2ab(i1,i3,k,j3) Rxcumul2aa(k,i2,j1,j2).
    ! All eight terms of an element share its j3 (a column index of p1 and p3,
    ! a row index of p5): one thread forms the products of a j3 slice (n^5 per
    ! thread) with BLAS inside the parallel region (single-threaded there with
    ! the OpenMP build of OpenBLAS). Every element is computed; the sum is
    ! antisymmetric in (i1,i2) (p1 alone, p3 and p5 in pairs) and, for a
    ! singlet rdm2ab, in (j1,j2) (p5 alone, p1 and p3 in pairs).
    allocate (cumul2abPerm(nOa, nOa, nOa, nOa))
    do k = 1, nOa
      cumul2abPerm(:, :, :, k) = cumul2ab(:, :, k, :)
    end do

    !$omp parallel private(p1t, p3t, p5t, b3, j3, j2, j1, i3, i2, i1)
    allocate (p1t(nOa, nOa, nOa, nOa, nOa), p3t(nOa, nOa, nOa, nOa, nOa), p5t(nOa, nOa, nOa, nOa, nOa))
    allocate (b3(nOa, nOa, nOa))
    !$omp do schedule(dynamic)
    do j3 = 1, nOa
      ! p1t(i1,i2,j1, i3,j2) = p1(i1,i2,j1, i3,j2,j3)
      call ZGEMM('N', 'N', nOa**3, nOa**2, nOa, (1.0_R64, 0.0_R64), cumul2aa, nOa**3, &
                 Rxcumul2ab(1, 1, 1, j3), nOa, (0.0_R64, 0.0_R64), p1t, nOa**3)
      ! p3t(i1,i3,j1, i2,j2) = p3(i1,i3,j1, i2,j3,j2)
      b3(:, :, :) = Rxcumul2ab(:, :, j3, :)
      call ZGEMM('N', 'N', nOa**3, nOa**2, nOa, (1.0_R64, 0.0_R64), cumul2ab, nOa**3, &
                 b3, nOa, (0.0_R64, 0.0_R64), p3t, nOa**3)
      ! p5t(i1,i3, i2,j1,j2) = p5(i1,i3,j3, i2,j1,j2): the rows (i1,i3) of slice j3
      call ZGEMM('N', 'N', nOa**2, nOa**3, nOa, (1.0_R64, 0.0_R64), cumul2abPerm(1, 1, j3, 1), nOa**3, &
                 Rxcumul2aa, nOa, (0.0_R64, 0.0_R64), p5t, nOa**2)
      do j2 = 1, nOa
        do j1 = 1, nOa
          do i3 = 1, nOa
            do i2 = 1, nOa
              do i1 = 1, nOa
                rdm3aab(i1, i2, i3, j1, j2, j3) = rdm3aab(i1, i2, i3, j1, j2, j3) + &
                                                  p1t(i1, i2, j1, i3, j2) - p1t(i1, i2, j2, i3, j1) + &
                                                  p3t(i1, i3, j1, i2, j2) - p3t(i1, i3, j2, i2, j1) + &
                                                  p5t(i1, i3, i2, j1, j2) - &
                                                  p3t(i2, i3, j1, i1, j2) + p3t(i2, i3, j2, i1, j1) - &
                                                  p5t(i2, i3, i1, j1, j2)
              end do
            end do
          end do
        end do
      end do
    end do
    !$omp end do
    deallocate (p1t, p3t, p5t, b3)
    !$omp end parallel

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end module
