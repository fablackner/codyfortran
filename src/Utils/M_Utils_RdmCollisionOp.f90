! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Trace and contraction utilities for reduced density matrices.
!>
!> Implements trace identities, consistency checks, and partial-trace
!> operations used to reduce higher-order RDMs and verify constraints.
module M_Utils_RdmCollisionOp
  use M_Utils_Types

  implicit none

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmCollisionOp_FillCollisionOp2ab(collOp2ab, h2, rdm3aab)
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    complex(R64), intent(out), allocatable :: collOp2ab(:, :, :, :)
    complex(R64), intent(in), contiguous  :: h2(:, :, :, :)
    complex(R64), intent(in), contiguous  :: rdm3aab(:, :, :, :, :, :)

    complex(R64), allocatable :: helpOp2ab(:, :, :, :)

    integer(I32) :: nOa, nOb
    integer(I32) :: i1, i2, j1, j2, k1, k2, k3

    nOa = Method_Mb_OrbBased_nOrbs(1)
    nOb = Method_Mb_OrbBased_nOrbs(2)

    ! Layout: (i1, i2, j1, j2) = (nOa, nOb, nOa, nOb)
    ! i1,j1,k1 are α indices (range nOa); i2,j2,k2,k3 are β indices (range nOb)
    if (.not. allocated(collOp2ab)) allocate (collOp2ab(nOa, nOb, nOa, nOb))
    allocate (helpOp2ab(nOa, nOb, nOa, nOb))

    helpOp2ab(:, :, :, :) = 0.0_R64

    do i1 = 1, nOa
      do i2 = 1, nOb
        do j1 = 1, nOa
          do j2 = 1, nOb

            do k1 = 1, nOa
              do k2 = 1, nOb
                do k3 = 1, nOb

                  helpOp2ab(i1, i2, j1, j2) = helpOp2ab(i1, i2, j1, j2) + &
                                              2.0_R64 * h2(i1, k3 + nOa, k1, k2 + nOa) * &
                                              (rdm3aab(k1, k2, i2, j1, k3, j2) + rdm3aab(i2, k2, k1, j2, k3, j1))

                end do
              end do
            end do

          end do
        end do
      end do
    end do

    collOp2ab(:, :, :, :) = 0.0_R64

    do j2 = 1, nOb
      do j1 = 1, nOa
        do i2 = 1, nOb
          do i1 = 1, nOa

            collOp2ab(i1, i2, j1, j2) = collOp2ab(i1, i2, j1, j2) + &
                                        helpOp2ab(i1, i2, j1, j2) + helpOp2ab(i2, i1, j2, j1) - &
                                        conjg(helpOp2ab(j1, j2, i1, i2) + helpOp2ab(j2, j1, i2, i1))

          end do
        end do
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmCollisionOp_FillCollisionOp2(collOp2, h2, rdm3)
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    complex(R64), intent(out), allocatable :: collOp2(:, :, :, :)
    complex(R64), intent(in), contiguous  :: h2(:, :, :, :)
    complex(R64), intent(in), contiguous  :: rdm3(:, :, :, :, :, :)

    integer(I32) :: nO
    integer :: i1, i2, j1, j2, k1, k2, k3, l3

    nO = Method_Mb_OrbBased_nOrbsSum

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
  subroutine RdmCollisionOp_FillCollisionOp1(collOp1, h2, rdm2)
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    complex(R64), intent(out), allocatable :: collOp1(:, :)
    complex(R64), intent(in), contiguous  :: h2(:, :, :, :)
    complex(R64), intent(in), contiguous  :: rdm2(:, :, :, :)

    integer(I32) :: nO
    integer :: i1, j1, k1, k2, l2

    nO = Method_Mb_OrbBased_nOrbsSum

    if (.not. allocated(collOp1)) allocate (collOp1(nO, nO))

    collOp1(:, :) = 0.0_R64

    do j1 = 1, nO
      do i1 = 1, nO

        !==========================================================================================
        ! Tr2 [W12, D12]
        !==========================================================================================

        do l2 = 1, nO
          do k2 = 1, nO
            do k1 = 1, nO
              collOp1(i1, j1) = collOp1(i1, j1) + 2.0_R64 * h2(i1, l2, k1, k2) * rdm2(k1, k2, j1, l2)
              collOp1(i1, j1) = collOp1(i1, j1) - 2.0_R64 * h2(k1, k2, j1, l2) * rdm2(i1, l2, k1, k2)
            end do
          end do
        end do

      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine RdmCollisionOp_FillCollisionOp1a(collOp1a, h2, rdm2ab)
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    complex(R64), intent(out), allocatable :: collOp1a(:, :)
    complex(R64), intent(in), contiguous  :: h2(:, :, :, :)
    complex(R64), intent(in), contiguous  :: rdm2ab(:, :, :, :)

    integer(I32) :: nOa, nOb
    integer :: i1, j1, k1, k2, l2

    nOa = Method_Mb_OrbBased_nOrbs(1)
    nOb = Method_Mb_OrbBased_nOrbs(2)

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
