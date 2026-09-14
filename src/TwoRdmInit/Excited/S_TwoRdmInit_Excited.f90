! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file S_TwoRdmInit_Excited.f90
!> @brief Implementation of excited-state 2-RDM initialization.
!>
!> @details
!> Builds the 2-RDM from 1-RDM products using the mean-field factorization:
!>
!> @par Same-Spin Blocks (bt1 == bt2)
!> Includes fermionic exchange (antisymmetrization):
!> @code
!> D^2(i1,i2;j1,j2) = D^1(i1,j1) * D^1(i2,j2) - D^1(i2,j1) * D^1(i1,j2)
!> @endcode
!>
!> @par Opposite-Spin Blocks (bt1 /= bt2)
!> Simple product (no exchange between different species):
!> @code
!> D^2(i1,i2;j1,j2) = D^1(i1,j1) * D^1(i2,j2)
!> @endcode
!>
!> @par Algorithm
!>   1. Build diagonal 1-RDM with occupation n_i = 1 for i <= nBodies
!>   2. Modify occupations based on creates/destroys configuration
!>   3. Compute 2-RDM as product of 1-RDMs with appropriate symmetry
submodule(M_TwoRdmInit_Excited) S_TwoRdmInit_Excited

  implicit none

  !-----------------------------------------------------------------------------
  ! Module-level configuration (populated by Fabricate)
  !-----------------------------------------------------------------------------

  !> Orbital indices to be populated (creation operators applied)
  integer(I32), allocatable :: creates(:)

  !> Orbital indices to be vacated (annihilation operators applied)
  integer(I32), allocatable :: destroys(:)

  !> Body type index for first particle species (default: 1 = alpha)
  integer(I32) :: bodyType1 = 1

  !> Body type index for second particle species (default: 2 = beta)
  integer(I32) :: bodyType2 = 2

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Read configuration and bind procedure pointers.
!>
!> @details
!> Reads JSON parameters from `"twoRdmInit.excited"` and binds
!> `TwoRdmInit_InitializeBlock` to the local `InitializeBlock` routine.
  module subroutine TwoRdmInit_Excited_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_TwoRdmInit

    call Say_Fabricate("twoRdmInit.excited")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    creates = Json_Get("twoRdmInit.excited.creates", [1])
    destroys = Json_Get("twoRdmInit.excited.destroys", [1])
    bodyType1 = Json_Get("twoRdmInit.excited.bodyType1", 1)
    bodyType2 = Json_Get("twoRdmInit.excited.bodyType2", 2)

    ! Map the local implementation to the main procedure pointer
    TwoRdmInit_InitializeBlock => InitializeBlock

    !------------------------------------
    ! branch
    !------------------------------------

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Initialize a 2-RDM block from 1-RDM products.
!>
!> @details
!> Constructs the 2-RDM block D^2_{bt1,bt2} using the mean-field
!> factorization ansatz. The 1-RDM for each body type is built as:
!>
!>   1. Start with diagonal occupations: n_i = 1 for i <= nBodies, else 0
!>   2. Apply modifications: n[destroy] -= 1, n[create] += 1
!>
!> @par Same-Spin (bt1 == bt2)
!> Fermionic antisymmetry enforced:
!> @code
!> D^2(i1,i2;j1,j2) = D^1(i1,j1)*D^1(i2,j2) - D^1(i2,j1)*D^1(i1,j2)
!> @endcode
!>
!> @par Opposite-Spin (bt1 /= bt2)
!> No exchange term (distinguishable particles):
!> @code
!> D^2(i1,i2;j1,j2) = D^1(i1,j1)*D^1(i2,j2)
!> @endcode
!>
!> @param[out] rdm2Bt  Output tensor block, shape (nOrbs1, nOrbs2, nOrbs1, nOrbs2)
!> @param[in]  bt1     Body type for indices i1, j1
!> @param[in]  bt2     Body type for indices i2, j2
!>
!> @note The local initDestroys/initCreates arrays are currently hardcoded
!>       to [1], meaning the module-level creates/destroys from JSON are
!>       read but not fully utilized. This appears to be incomplete.
!>       TODO: Wire up the module-level creates/destroys properly.
  subroutine InitializeBlock(rdm2Bt, bt1, bt2)
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    complex(R64), intent(out), contiguous  :: rdm2Bt(:, :, :, :)
    integer(I32), intent(in)               :: bt1
    integer(I32), intent(in)               :: bt2

    complex(R64), allocatable :: rdm1bt1(:, :)
    complex(R64), allocatable :: rdm1bt2(:, :)

    integer(I32) :: i, i1, j1, i2, j2
    integer(I32) :: iDestroy, iCreate
    integer(I32), allocatable :: initDestroysBt1(:)
    integer(I32), allocatable :: initCreatesBt1(:)
    integer(I32), allocatable :: initDestroysBt2(:)
    integer(I32), allocatable :: initCreatesBt2(:)

    allocate (initDestroysBt1, source=[1])
    allocate (initCreatesBt1, source=[1])
    allocate (initDestroysBt2, source=[1])
    allocate (initCreatesBt2, source=[1])

    allocate (rdm1bt1(Method_Mb_OrbBased_nOrbs(bt1), Method_Mb_OrbBased_nOrbs(bt1)))

    rdm1bt1(:, :) = 0.0_R64
    do i = 1, Method_Mb_nBodies(bt1)
      rdm1bt1(i, i) = (1.0_R64, 0.0_R64)
    end do

    do i = 1, size(initDestroysBt1)
      iDestroy = initDestroysBt1(i)
      rdm1bt1(iDestroy, iDestroy) = rdm1bt1(iDestroy, iDestroy) - 1.0_R64
    end do

    do i = 1, size(initCreatesBt1)
      iCreate = initCreatesBt1(i)
      rdm1bt1(iCreate, iCreate) = rdm1bt1(iCreate, iCreate) + 1.0_R64
    end do

    allocate (rdm1bt2(Method_Mb_OrbBased_nOrbs(bt2), Method_Mb_OrbBased_nOrbs(bt2)))

    rdm1bt2(:, :) = 0.0_R64
    do i = 1, Method_Mb_nBodies(bt2)
      rdm1bt2(i, i) = (1.0_R64, 0.0_R64)
    end do

    do i = 1, size(initDestroysBt2)
      iDestroy = initDestroysBt2(i)
      rdm1bt2(iDestroy, iDestroy) = rdm1bt2(iDestroy, iDestroy) - 1.0_R64
    end do

    do i = 1, size(initCreatesBt2)
      iCreate = initCreatesBt2(i)
      rdm1bt2(iCreate, iCreate) = rdm1bt2(iCreate, iCreate) + 1.0_R64
    end do

    if (bt1 .eq. bt2) then

      do j2 = 1, Method_Mb_OrbBased_nOrbs(bt2)
        do j1 = 1, Method_Mb_OrbBased_nOrbs(bt1)
          do i2 = 1, Method_Mb_OrbBased_nOrbs(bt2)
            do i1 = 1, Method_Mb_OrbBased_nOrbs(bt1)

              rdm2Bt(i1, i2, j1, j2) = rdm1bt1(i1, j1) * rdm1bt2(i2, j2) - rdm1bt2(i2, j1) * rdm1bt1(i1, j2)

            end do
          end do
        end do
      end do

    else

      do j2 = 1, Method_Mb_OrbBased_nOrbs(bt2)
        do j1 = 1, Method_Mb_OrbBased_nOrbs(bt1)
          do i2 = 1, Method_Mb_OrbBased_nOrbs(bt2)
            do i1 = 1, Method_Mb_OrbBased_nOrbs(bt1)

              rdm2Bt(i1, i2, j1, j2) = rdm1bt1(i1, j1) * rdm1bt2(i2, j2)

            end do
          end do
        end do
      end do

    end if

  end subroutine

end submodule
