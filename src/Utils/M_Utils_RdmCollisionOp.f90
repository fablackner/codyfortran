! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Collision operators of the BBGKY hierarchy, Tr3[W, D3] for the 2-RDM and
!> Tr2[W, D2] for the 1-RDM equation, in the spin-orbital and singlet (aab)
!> layouts; sizes from the arguments.
module M_Utils_RdmCollisionOp
  use M_Utils_Types

  implicit none

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Tr3[W, D3] of the singlet 2-RDM equation, collOp2ab (n^4), from the
  !> ket-first rdm3aab (n^6) and h2 in spin orbitals (up block first, 2n per
  !> index): one ZGEMM of h2 as (n x n^3).
  subroutine RdmCollisionOp_FillCollisionOp2ab(collOp2ab, h2, rdm3aab)
    use M_Utils_BlasLib

    complex(R64), intent(out), allocatable :: collOp2ab(:, :, :, :)
    complex(R64), intent(in), contiguous  :: h2(:, :, :, :)
    complex(R64), intent(in), contiguous  :: rdm3aab(:, :, :, :, :, :)

    complex(R64), allocatable :: helpOp2ab(:, :, :, :), hMat(:, :), x(:, :, :, :, :, :)

    integer(I32) :: n
    integer(I32) :: i1, i2, j1, j2, k1, k2, k3

    n = size(rdm3aab, 1)
    if (size(h2, 1) .ne. 2 * n) error stop "RdmCollisionOp_FillCollisionOp2ab: h2 must span the up and down orbitals of rdm3aab"

    ! helpOp2ab(i1, i2, j1, j2) = 2 sum_{k1 k2 k3} h2(i1, k3 + n, k1, k2 + n)
    !   (rdm3aab(k1, k2, i2, j1, k3, j2) + rdm3aab(i2, k2, k1, j2, k3, j1)),
    ! one matrix product hMat(i1, (k1 k2 k3)) x((k1 k2 k3), (i2 j1 j2))
    allocate (hMat(n, n**3), x(n, n, n, n, n, n))

    do k3 = 1, n
      do k2 = 1, n
        do k1 = 1, n
          hMat(:, k1 + n * (k2 - 1) + n * n * (k3 - 1)) = 2.0_R64 * h2(1:n, k3 + n, k1, k2 + n)
        end do
      end do
    end do

    !$omp parallel do collapse(2) private(j2, j1, i2, k3, k2, k1)
    do j2 = 1, n
      do j1 = 1, n
        do i2 = 1, n
          do k3 = 1, n
            do k2 = 1, n
              do k1 = 1, n
                x(k1, k2, k3, i2, j1, j2) = rdm3aab(k1, k2, i2, j1, k3, j2) + rdm3aab(i2, k2, k1, j2, k3, j1)
              end do
            end do
          end do
        end do
      end do
    end do
    !$omp end parallel do

    allocate (helpOp2ab(n, n, n, n))
    call ZGEMM('N', 'N', n, n**3, n**3, (1.0_R64, 0.0_R64), hMat, n, x, n**3, (0.0_R64, 0.0_R64), helpOp2ab, n)

    if (.not. allocated(collOp2ab)) allocate (collOp2ab(n, n, n, n))

    do j2 = 1, n
      do j1 = 1, n
        do i2 = 1, n
          do i1 = 1, n

            collOp2ab(i1, i2, j1, j2) = helpOp2ab(i1, i2, j1, j2) + helpOp2ab(i2, i1, j2, j1) - &
                                        conjg(helpOp2ab(j1, j2, i1, i2) + helpOp2ab(j2, j1, i2, i1))

          end do
        end do
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmCollisionOp_FillCollisionOp2(collOp2, h2, rdm3)

    complex(R64), intent(out), allocatable :: collOp2(:, :, :, :)
    complex(R64), intent(in), contiguous  :: h2(:, :, :, :)
    complex(R64), intent(in), contiguous  :: rdm3(:, :, :, :, :, :)

    integer(I32) :: nO
    integer :: i1, i2, j1, j2, k1, k2, k3, l3

    nO = size(rdm3, 1)

    if (.not. allocated(collOp2)) allocate (collOp2(nO, nO, nO, nO))

    collOp2(:, :, :, :) = 0.0_R64

    do j2 = 1, nO
      do j1 = 1, nO

        do i2 = 1, nO
          do i1 = 1, nO

            !==========================================================================================
            ! Tr3 [W13+W23, D123]
            !==========================================================================================

            do l3 = 1, nO
              do k3 = 1, nO

                do k1 = 1, nO
                  collOp2(i1, i2, j1, j2) = collOp2(i1, i2, j1, j2) + 2.0_R64 * h2(i1, l3, k1, k3) * rdm3(k1, i2, k3, j1, j2, l3)
                  collOp2(i1, i2, j1, j2) = collOp2(i1, i2, j1, j2) - 2.0_R64 * h2(k1, k3, j1, l3) * rdm3(i1, i2, l3, k1, j2, k3)
                end do

                do k2 = 1, nO
                  collOp2(i1, i2, j1, j2) = collOp2(i1, i2, j1, j2) + 2.0_R64 * h2(i2, l3, k2, k3) * rdm3(i1, k2, k3, j1, j2, l3)
                  collOp2(i1, i2, j1, j2) = collOp2(i1, i2, j1, j2) - 2.0_R64 * h2(k2, k3, j2, l3) * rdm3(i1, i2, l3, j1, k2, k3)
                end do

              end do
            end do

          end do
        end do
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmCollisionOp_FillCollisionOp1a(collOp1a, h2, rdm2ab)

    complex(R64), intent(out), allocatable :: collOp1a(:, :)
    complex(R64), intent(in), contiguous  :: h2(:, :, :, :)
    complex(R64), intent(in), contiguous  :: rdm2ab(:, :, :, :)

    integer(I32) :: nOa, nOb
    integer :: i1, j1, k1, k2, l2

    nOa = size(rdm2ab, 1)
    nOb = size(rdm2ab, 2)

    if (.not. allocated(collOp1a)) allocate (collOp1a(nOa, nOa))

    collOp1a(:, :) = 0.0_R64

    do j1 = 1, nOa
      do i1 = 1, nOa

        !==========================================================================================
        ! Tr2 [W12, D12] for alpha block
        !==========================================================================================

        ! Contribution from alpha-beta interaction
        do l2 = 1, nOb
          do k2 = 1, nOb
            do k1 = 1, nOa
              ! h2 indices: (i1, l2+nOa, k1, k2+nOa) -> v_{i1_a, l2_b, k1_a, k2_b}
              ! rdm2ab indices: (k1, k2, j1, l2) -> D_{k1_a, k2_b, j1_a, l2_b}
              collOp1a(i1, j1) = collOp1a(i1, j1) + 2.0_R64 * h2(i1, l2 + nOa, k1, k2 + nOa) * rdm2ab(k1, k2, j1, l2)
              collOp1a(i1, j1) = collOp1a(i1, j1) - 2.0_R64 * h2(k1, k2 + nOa, j1, l2 + nOa) * rdm2ab(i1, l2, k1, k2)
            end do
          end do
        end do

        ! Contribution from alpha-alpha interaction (using singlet relation D_aa from D_ab)
        ! D_aa(k1, k2, j1, l2) = D_ab(k1, k2, j1, l2) - D_ab(k1, k2, l2, j1) ?
        ! For singlet: D_{pqrs}^{aa} = D_{pqrs}^{ab} - D_{pqsr}^{ab}
        do l2 = 1, nOa
          do k2 = 1, nOa
            do k1 = 1, nOa
              ! h2 indices: (i1, l2, k1, k2) -> v_{i1_a, l2_a, k1_a, k2_a}
              ! D_aa term constructed from rdm2ab
              collOp1a(i1, j1) = collOp1a(i1, j1) + 2.0_R64 * h2(i1, l2, k1, k2) * (rdm2ab(k1, k2, j1, l2) - rdm2ab(k1, k2, l2, j1))
              collOp1a(i1, j1) = collOp1a(i1, j1) - 2.0_R64 * h2(k1, k2, j1, l2) * (rdm2ab(i1, l2, k1, k2) - rdm2ab(i1, l2, k2, k1))
            end do
          end do
        end do

      end do
    end do

  end subroutine

end module
