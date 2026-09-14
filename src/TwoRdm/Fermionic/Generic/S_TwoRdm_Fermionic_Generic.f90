! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Implementation submodule for M_TwoRdm_Fermionic_Generic.
!>
!> @details
!> Provides the concrete implementations for the generic (unrestricted) fermionic
!> 2-RDM backend. All operations work on the full nOrbs⁴ tensor without exploiting
!> spin symmetry.
submodule(M_TwoRdm_Fermionic_Generic) S_TwoRdm_Fermionic_Generic

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Fabricate the generic fermionic backend.
!>
!> @details
!> Sets TwoRdm_nElements to nOrbs⁴ and wires all procedure pointers to the
!> generic implementations.
  module subroutine TwoRdm_Fermionic_Generic_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Method_Mb
    use M_Method_Mb_OrbBased
    use M_TwoRdm

    call Say_Fabricate("twoRdm.fermionic.generic")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    TwoRdm_nElements = Method_Mb_OrbBased_nOrbsSum**4

    TwoRdm_ApplyCommutatorH1 => ApplyCommutatorH1
    TwoRdm_ApplyCommutatorH2 => ApplyCommutatorH2

    TwoRdm_FillRdm1FromTwoRdmState => FillRdm1FromTwoRdmState
    TwoRdm_FillRdm2FromTwoRdmState => FillRdm2FromTwoRdmState
    TwoRdm_FillTwoRdmStateFromRdm2 => FillTwoRdmStateFromRdm2

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Accumulate [ĥ₁, D²] into the time-derivative of the 2-RDM.
!>
!> @details
!> Computes the commutator of the one-body Hamiltonian with the 2-RDM:
!>   [ĥ₁, D²]_{i₁i₂;j₁j₂} = Σ_k (h_{i₁k}D²_{ki₂;j₁j₂} + h_{i₂k}D²_{i₁k;j₁j₂}
!>                               - h_{kj₁}D²_{i₁i₂;kj₂} - h_{kj₂}D²_{i₁i₂;j₁k})
!>
!> The result is accumulated (added) to dTwoRdmState, which must be pre-initialized.
!>
!> @param[inout] dTwoRdmState  Accumulator for time derivative (pre-initialized)
!> @param[in]    twoRdmState   Input linearized 2-RDM
!> @param[in]    h1            One-body Hamiltonian matrix (nOrbs × nOrbs)
  subroutine ApplyCommutatorH1(dTwoRdmState, twoRdmState, h1)
    use M_Method_Mb
    use M_Method_Mb_OrbBased
    use M_TwoRdm

    complex(R64), intent(in), contiguous  :: twoRdmState(:)
    complex(R64), intent(inout), contiguous  :: dTwoRdmState(:)
    complex(R64), intent(in), contiguous  :: h1(:, :)

    complex(R64), allocatable :: rdm2(:, :, :, :)
    complex(R64), allocatable :: dRdm2(:, :, :, :)

    integer(I32) :: nO
    integer(I32) :: i1, i2, j1, j2, k

    call TwoRdm_FillRdm2FromTwoRdmState(rdm2, twoRdmState)
    call TwoRdm_FillRdm2FromTwoRdmState(dRdm2, dTwoRdmState)

    nO = Method_Mb_OrbBased_nOrbsSum

    do j2 = 1, nO
      do j1 = 1, nO

        do i2 = 1, nO
          do i1 = 1, nO

            ! [ĥ₁, D²]: contract h1 with each index of D²
            do k = 1, nO

              dRdm2(i1, i2, j1, j2) = dRdm2(i1, i2, j1, j2) + h1(i1, k) * rdm2(k, i2, j1, j2)
              dRdm2(i1, i2, j1, j2) = dRdm2(i1, i2, j1, j2) + h1(i2, k) * rdm2(i1, k, j1, j2)

              dRdm2(i1, i2, j1, j2) = dRdm2(i1, i2, j1, j2) - h1(k, j1) * rdm2(i1, i2, k, j2)
              dRdm2(i1, i2, j1, j2) = dRdm2(i1, i2, j1, j2) - h1(k, j2) * rdm2(i1, i2, j1, k)

            end do

          end do
        end do
      end do
    end do

    call TwoRdm_FillTwoRdmStateFromRdm2(dTwoRdmState, dRdm2)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Accumulate [Ŵ₁₂, D²] into the time-derivative of the 2-RDM.
!>
!> @details
!> Computes the commutator of the two-body interaction with the 2-RDM:
!>   [Ŵ₁₂, D²]_{i₁i₂;j₁j₂} = 2·Σ_{k₁k₂} (W_{i₁i₂;k₁k₂}D²_{k₁k₂;j₁j₂}
!>                                       - W_{k₁k₂;j₁j₂}D²_{i₁i₂;k₁k₂})
!>
!> The factor of 2 accounts for the symmetrized form of the two-body operator.
!>
!> @param[inout] dTwoRdmState  Accumulator for time derivative (pre-initialized)
!> @param[in]    twoRdmState   Input linearized 2-RDM
!> @param[in]    h2            Two-body interaction tensor (nOrbs⁴)
  subroutine ApplyCommutatorH2(dTwoRdmState, twoRdmState, h2)
    use M_Method_Mb
    use M_Method_Mb_OrbBased
    use M_TwoRdm

    complex(R64), intent(in), contiguous  :: twoRdmState(:)
    complex(R64), intent(inout), contiguous  :: dTwoRdmState(:)
    complex(R64), intent(in), contiguous  :: h2(:, :, :, :)

    complex(R64), allocatable :: rdm2(:, :, :, :)
    complex(R64), allocatable :: dRdm2(:, :, :, :)

    integer(I32) :: nO
    integer(I32) :: i1, i2, j1, j2, k1, k2

    call TwoRdm_FillRdm2FromTwoRdmState(rdm2, twoRdmState)
    call TwoRdm_FillRdm2FromTwoRdmState(dRdm2, dTwoRdmState)

    nO = Method_Mb_OrbBased_nOrbsSum

    do j2 = 1, nO
      do j1 = 1, nO

        do i2 = 1, nO
          do i1 = 1, nO

            ! [Ŵ₁₂, D²]: double contraction over k₁,k₂
            do k2 = 1, nO
              do k1 = 1, nO

                dRdm2(i1, i2, j1, j2) = dRdm2(i1, i2, j1, j2) + 2.0_R64 * h2(i1, i2, k1, k2) * rdm2(k1, k2, j1, j2)
                dRdm2(i1, i2, j1, j2) = dRdm2(i1, i2, j1, j2) - 2.0_R64 * h2(k1, k2, j1, j2) * rdm2(i1, i2, k1, k2)

              end do
            end do

          end do
        end do
      end do
    end do

    call TwoRdm_FillTwoRdmStateFromRdm2(dTwoRdmState, dRdm2)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Compute the 1-RDM from the 2-RDM via partial trace.
!>
!> @details
!> The 1-RDM is obtained by tracing out one particle from the 2-RDM:
!>   D¹(i,j) = (1/(N-1)) · Σ_k D²(i,k;j,k)
!>
!> The normalization factor ensures Tr(D¹) = N.
!>
!> @param[out] rdm1         Allocated on output: 1-RDM (nOrbs × nOrbs)
!> @param[in]  twoRdmState  Input linearized 2-RDM
  subroutine FillRdm1FromTwoRdmState(rdm1, twoRdmState)
    use M_Method_Mb
    use M_Method_Mb_OrbBased
    use M_TwoRdm

    complex(R64), intent(out), allocatable :: rdm1(:, :)
    complex(R64), intent(in), contiguous                :: twoRdmState(:)

    complex(R64), allocatable :: rdm2(:, :, :, :)

    integer(I32) :: i, j, k, nE, nO

    nO = Method_Mb_OrbBased_nOrbsSum
    nE = TwoRdm_nElements

    call TwoRdm_FillRdm2FromTwoRdmState(rdm2, twoRdmState)

    if (.not. allocated(rdm1)) allocate (rdm1(nO, nO))

    ! Partial trace: D¹(i,j) = Σ_k D²(i,k;j,k)
    rdm1(:, :) = 0.0_R64
    do j = 1, nO
      do i = 1, nO

        do k = 1, nO
          rdm1(i, j) = rdm1(i, j) + rdm2(i, k, j, k)
        end do

      end do
    end do

    ! Normalize by (N-1)
    rdm1(:, :) = rdm1(:, :) / (Method_Mb_nBodiesSum - 1)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Unpack the linearized 2-RDM into a 4-index tensor.
!>
!> @details
!> Reconstructs the full 4-index 2-RDM tensor from the linearized storage.
!> Column-major ordering: i₁ varies fastest, then i₂, j₁, j₂.
!>
!> @param[out] rdm2         Allocated on output: 2-RDM tensor (nOrbs⁴)
!> @param[in]  twoRdmState  Input linearized 2-RDM
  subroutine FillRdm2FromTwoRdmState(rdm2, twoRdmState)
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    complex(R64), intent(out), allocatable :: rdm2(:, :, :, :)
    complex(R64), intent(in), contiguous                :: twoRdmState(:)

    integer(I32) :: nO, count
    integer(I32) :: i1, i2, j1, j2

    nO = Method_Mb_OrbBased_nOrbsSum

    if (.not. allocated(rdm2)) allocate (rdm2(nO, nO, nO, nO))

    count = 0
    do j2 = 1, nO
      do j1 = 1, nO
        do i2 = 1, nO
          do i1 = 1, nO

            count = count + 1

            rdm2(i1, i2, j1, j2) = twoRdmState(count)

          end do
        end do
      end do
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Pack a 4-index 2-RDM tensor into linearized storage.
!>
!> @details
!> Serializes the 4-index tensor to a 1D array using column-major ordering.
!> The output array must be pre-allocated with size nOrbs⁴.
!>
!> @param[out] twoRdmState  Output linearized 2-RDM (pre-allocated)
!> @param[in]  rdm2         Input 2-RDM tensor (nOrbs⁴)
  subroutine FillTwoRdmStateFromRdm2(twoRdmState, rdm2)
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    complex(R64), intent(out), contiguous  :: twoRdmState(:)
    complex(R64), intent(in), contiguous   :: rdm2(:, :, :, :)

    integer(I32) :: nO, count
    integer(I32) :: i1, i2, j1, j2

    nO = Method_Mb_OrbBased_nOrbsSum
    count = 0

    do j2 = 1, nO
      do j1 = 1, nO
        do i2 = 1, nO
          do i1 = 1, nO

            count = count + 1

            twoRdmState(count) = rdm2(i1, i2, j1, j2)

          end do
        end do
      end do
    end do

  end subroutine

end submodule
