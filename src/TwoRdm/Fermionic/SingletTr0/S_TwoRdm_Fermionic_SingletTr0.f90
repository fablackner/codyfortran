! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Implementation submodule for M_TwoRdm_Fermionic_SingletTr0.
!>
!> @details
!> Provides the trace-free singlet backend that stores the 1-RDM explicitly
!> alongside the trace-free αβ 2-RDM block. This representation is particularly
!> suited for TD-2RDM propagation where maintaining trace consistency is critical.
submodule(M_TwoRdm_Fermionic_SingletTr0) S_TwoRdm_Fermionic_SingletTr0

  use M_Utils_RdmEnforceTraceM2ab

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Fabricate the singletTr0 fermionic backend.
!>
!> @details
!> Sets TwoRdm_nElements to nα² + nα²×nβ² to accommodate the explicit 1-RDM
!> storage plus the trace-free αβ block.
  module subroutine TwoRdm_Fermionic_SingletTr0_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Method_Mb
    use M_Method_Mb_OrbBased
    use M_TwoRdm

    call Say_Fabricate("twoRdm.fermionic.singletTr0")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    ! Storage: nα² (1-RDM) + nα²×nβ² (trace-free 2-RDM)
    TwoRdm_nElements = Method_Mb_OrbBased_nOrbs(1)**2 + &
                       (Method_Mb_OrbBased_nOrbs(1)**2) * (Method_Mb_OrbBased_nOrbs(2)**2)

    TwoRdm_ApplyCommutatorH1 => ApplyCommutatorH1
    TwoRdm_ApplyCommutatorH2 => ApplyCommutatorH2

    TwoRdm_FillRdm1FromTwoRdmState => FillRdm1FromTwoRdmState
    TwoRdm_FillRdm2FromTwoRdmState => FillRdm2FromTwoRdmState
    TwoRdm_FillTwoRdmStateFromRdm2 => FillTwoRdmStateFromRdm2

    TwoRdm_Fermionic_SingletTr0_FillTwoRdmStateFromRdm2ab => FillTwoRdmStateFromRdm2ab
    TwoRdm_Fermionic_SingletTr0_FillTwoRdmStateFromRdms => FillTwoRdmStateFromRdms
    TwoRdm_Fermionic_SingletTr0_FillRdmsFromTwoRdmState => FillRdmsFromTwoRdmState
    TwoRdm_Fermionic_SingletTr0_FillRdm2abFromTwoRdmState => FillRdm2abFromTwoRdmState
    TwoRdm_Fermionic_SingletTr0_FillRdm1aFromTwoRdmState => FillRdm1aFromTwoRdmState

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Accumulate [ĥ₁, D²] and [ĥ₁, D¹] for the trace-free representation.
!>
!> @details
!> This commutator has two parts:
!> 1. The αβ 2-RDM commutator (same as Singlet backend)
!> 2. The 1-RDM commutator: [ĥ₁, D¹]_{ij} = Σ_k (h_{ik}D¹_{kj} - h_{kj}D¹_{ik})
!>
!> After computing, the trace-free constraint is re-enforced on the 2-RDM part.
!>
!> @param[inout] dTwoRdmState  Accumulator for time derivative
!> @param[in]    twoRdmState   Input state (1-RDM + trace-free 2-RDM)
!> @param[in]    h1            One-body Hamiltonian (spin-orbital basis)
  subroutine ApplyCommutatorH1(dTwoRdmState, twoRdmState, h1)
    use M_Method_Mb
    use M_Method_Mb_OrbBased
    use M_Utils_RdmEnforceTraceM2ab

    complex(R64), intent(inout), contiguous  :: dTwoRdmState(:)
    complex(R64), intent(in), contiguous     :: twoRdmState(:)
    complex(R64), intent(in), contiguous     :: h1(:, :)

    complex(R64), allocatable :: rdm2ab(:, :, :, :)
    complex(R64), allocatable :: dRdm2ab(:, :, :, :)
    complex(R64), allocatable :: rdm1a(:, :), drdm1a(:, :)

    integer(I32) :: nOa, nOb
    integer(I32) :: i1, i2, j1, j2, k

    nOa = Method_Mb_OrbBased_nOrbs(1)
    nOb = Method_Mb_OrbBased_nOrbs(2)

    ! Extract 1-RDM components for explicit commutator
    call TwoRdm_Fermionic_SingletTr0_FillRdmsFromTwoRdmState(twoRdmState, rdm1a_=rdm1a)
    call TwoRdm_Fermionic_SingletTr0_FillRdmsFromTwoRdmState(dTwoRdmState, rdm1a_=drdm1a)

    ! Get full αβ 2-RDM (trace restored)
    call TwoRdm_Fermionic_SingletTr0_FillRdm2abFromTwoRdmState(rdm2ab, twoRdmState)
    call TwoRdm_Fermionic_SingletTr0_FillRdm2abFromTwoRdmState(dRdm2ab, dTwoRdmState)

    ! Compute [ĥ₁, D²_αβ]
    do j2 = 1, nOb
      do j1 = 1, nOa

        do i2 = 1, nOb
          do i1 = 1, nOa

            do k = 1, nOa

              dRdm2ab(i1, i2, j1, j2) = dRdm2ab(i1, i2, j1, j2) + h1(i1, k) * rdm2ab(k, i2, j1, j2)
              dRdm2ab(i1, i2, j1, j2) = dRdm2ab(i1, i2, j1, j2) + h1(i2 + nOa, k + nOa) * rdm2ab(i1, k, j1, j2)

              dRdm2ab(i1, i2, j1, j2) = dRdm2ab(i1, i2, j1, j2) - h1(k, j1) * rdm2ab(i1, i2, k, j2)
              dRdm2ab(i1, i2, j1, j2) = dRdm2ab(i1, i2, j1, j2) - h1(k + nOa, j2 + nOa) * rdm2ab(i1, i2, j1, k)

            end do

          end do
        end do
      end do
    end do

    ! Compute [ĥ₁, D¹_α] explicitly for the 1-RDM segment
    do j1 = 1, nOa
      do i1 = 1, nOa
        do k = 1, nOa
          drdm1a(i1, j1) = drdm1a(i1, j1) + h1(i1, k) * rdm1a(k, j1) - h1(k, j1) * rdm1a(i1, k)
        end do
      end do
    end do

    ! Enforce trace-free constraint on 2-RDM part before repacking
    call RdmEnforceTraceM2ab_Apply(dRdm2ab)
    call TwoRdm_Fermionic_SingletTr0_FillTwoRdmStateFromRdms(dTwoRdmState, rdm2abTr0_=dRdm2ab, rdm1a_=drdm1a)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Accumulate [Ŵ₁₂, D²] for the trace-free representation.
!>
!> @details
!> The two-body commutator only affects the 2-RDM part (the 1-RDM is unchanged
!> by [Ŵ, D²] in the BBGKY hierarchy - that term comes from the 3-RDM closure).
!>
!> @param[inout] dTwoRdmState  Accumulator for time derivative
!> @param[in]    twoRdmState   Input state (1-RDM + trace-free 2-RDM)
!> @param[in]    h2            Two-body interaction (spin-orbital basis)
  subroutine ApplyCommutatorH2(dTwoRdmState, twoRdmState, h2)
    use M_Method_Mb
    use M_Method_Mb_OrbBased
    use M_Utils_RdmEnforceTraceM2ab

    complex(R64), intent(inout), contiguous  :: dTwoRdmState(:)
    complex(R64), intent(in), contiguous     :: twoRdmState(:)
    complex(R64), intent(in), contiguous     :: h2(:, :, :, :)

    complex(R64), allocatable :: rdm2ab(:, :, :, :)
    complex(R64), allocatable :: dRdm2ab(:, :, :, :)

    integer(I32) :: nOa, nOb
    integer(I32) :: i1, i2, j1, j2, k1, k2

    call TwoRdm_Fermionic_SingletTr0_FillRdm2abFromTwoRdmState(rdm2ab, twoRdmState)
    call TwoRdm_Fermionic_SingletTr0_FillRdm2abFromTwoRdmState(dRdm2ab, dTwoRdmState)

    nOa = Method_Mb_OrbBased_nOrbs(1)
    nOb = Method_Mb_OrbBased_nOrbs(2)

    do j2 = 1, nOb
      do j1 = 1, nOa

        do i2 = 1, nOb
          do i1 = 1, nOa

            ! [Ŵ₁₂, D²_αβ]
            do k2 = 1, nOb
              do k1 = 1, nOa

                dRdm2ab(i1, i2, j1, j2) = dRdm2ab(i1, i2, j1, j2) + &
                                          2.0_R64 * h2(i1, i2 + nOa, k1, k2 + nOa) * rdm2ab(k1, k2, j1, j2)
                dRdm2ab(i1, i2, j1, j2) = dRdm2ab(i1, i2, j1, j2) - &
                                          2.0_R64 * h2(k1, k2 + nOa, j1, j2 + nOa) * rdm2ab(i1, i2, k1, k2)

              end do
            end do

          end do
        end do
      end do
    end do

    ! Enforce trace-free and repack (1-RDM unchanged by [W,D²])
    call RdmEnforceTraceM2ab_Apply(dRdm2ab)
    call TwoRdm_Fermionic_SingletTr0_FillTwoRdmStateFromRdms(dTwoRdmState, rdm2abTr0_=dRdm2ab)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Extract the full 1-RDM from the stored state.
!>
!> @details
!> Returns the full spin-orbital 1-RDM by reading the explicit α block and
!> using singlet symmetry (D¹_β = D¹_α) to fill the β block.
!>
!> @param[out] rdm1         Allocated on output: full 1-RDM (nOrbs × nOrbs)
!> @param[in]  twoRdmState  Input state with explicit 1-RDM segment
  subroutine FillRdm1FromTwoRdmState(rdm1, twoRdmState)
    use M_Utils_RdmTrace
    use M_Method_Mb_OrbBased

    complex(R64), intent(out), allocatable :: rdm1(:, :)
    complex(R64), intent(in), contiguous  :: twoRdmState(:)

    integer(I32) :: nOa, nOb, nO
    integer(I32) :: i, j, offset

    nOa = Method_Mb_OrbBased_nOrbs(1)
    nOb = Method_Mb_OrbBased_nOrbs(2)
    nO = nOa + nOb

    if (.not. allocated(rdm1)) allocate (rdm1(nO, nO))
    rdm1 = 0.0_R64

    ! Read α block from first segment
    offset = 0
    do j = 1, nOa
      do i = 1, nOa
        offset = offset + 1
        rdm1(i, j) = twoRdmState(offset)
      end do
    end do

    ! Singlet symmetry: D¹_β = D¹_α
    rdm1(nOa + 1:nO, nOa + 1:nO) = rdm1(1:nOa, 1:nOa)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Reconstruct the full 2-RDM tensor from the stored state.
!>
!> @details
!> First reconstructs the full αβ block (adding trace contribution back),
!> then uses singlet relations for the same-spin blocks.
!>
!> @param[out] rdm2         Allocated on output: full 2-RDM (nOrbs⁴)
!> @param[in]  twoRdmState  Input state
  subroutine FillRdm2FromTwoRdmState(rdm2, twoRdmState)
    use M_Utils_RdmSingletRelations
    use M_Method

    complex(R64), intent(out), allocatable :: rdm2(:, :, :, :)
    complex(R64), intent(in), contiguous                :: twoRdmState(:)

    complex(R64), allocatable :: rdm2ab(:, :, :, :)

    call TwoRdm_Fermionic_SingletTr0_FillRdm2abFromTwoRdmState(rdm2ab, twoRdmState)
    call RdmSingletRelations_FillRdm2FromRdm2ab(rdm2, rdm2ab)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Convert a full 2-RDM to the trace-free singlet representation.
!>
!> @param[out] twoRdmState  Output state (pre-allocated)
!> @param[in]  rdm2         Input full 2-RDM tensor
  subroutine FillTwoRdmStateFromRdm2(twoRdmState, rdm2)
    use M_Utils_RdmTrace
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    complex(R64), intent(out), contiguous  :: twoRdmState(:)
    complex(R64), intent(in), contiguous   :: rdm2(:, :, :, :)

    complex(R64), allocatable :: rdm2ab(:, :, :, :)
    complex(R64), allocatable :: rdm1a(:, :)

    call RdmTrace_FillRdm2abFromRdm2(rdm2ab, rdm2, Method_Mb_OrbBased_nOrbs(1), Method_Mb_OrbBased_nOrbs(2))
    call RdmTrace_FillRdm1aFromRdm2ab(rdm1a, rdm2ab)
    call TwoRdm_Fermionic_SingletTr0_FillTwoRdmStateFromRdm2ab(twoRdmState, rdm2ab)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Reconstruct full αβ 2-RDM from trace-free storage + 1-RDM.
!>
!> @details
!> Adds back the trace contribution to the stored trace-free αβ block:
!>   D²_αβ(i₁,i₂;j₁,j₂) = D²_αβ,Tr0(i₁,i₂;j₁,j₂) + D¹_α(i₁,j₁)·δ_{i₂j₂}·Nβ/...
!>
!> The reconstruction uses RdmEnforceTraceM2ab_Apply with onefold traces.
!>
!> @param[out] rdm2ab       Allocated on output: full αβ tensor
!> @param[in]  twoRdmState  Input state
  subroutine FillRdm2abFromTwoRdmState(rdm2ab, twoRdmState)
    use M_Utils_RdmSingletRelations
    use M_Method_Mb
    use M_Method_Mb_OrbBased
    use M_Utils_RdmEnforceTraceM2ab

    complex(R64), intent(out), allocatable :: rdm2ab(:, :, :, :)
    complex(R64), intent(in), contiguous   :: twoRdmState(:)

    complex(R64), allocatable :: onefoldTraces(:, :, :)
    integer(I32) :: nOa, nOb
    integer(I32) :: count
    integer(I32) :: i1, i2, j1, j2
    complex(R64) :: val

    nOa = Method_Mb_OrbBased_nOrbs(1)
    nOb = Method_Mb_OrbBased_nOrbs(2)

    if (.not. allocated(rdm2ab)) allocate (rdm2ab(nOa, nOb, nOa, nOb))

    allocate (onefoldTraces(nOa, nOa, 2))

    ! Read 1-RDM segment and prepare onefold traces for reconstruction
    count = 0
    do j1 = 1, nOa
      do i1 = 1, nOa
        count = count + 1
        val = twoRdmState(count)
        ! Vertical contraction: Nβ · D¹_α
        onefoldTraces(i1, j1, 1) = val * Method_Mb_nBodies(2)
        ! Diagonal contraction: D¹_α
        onefoldTraces(i1, j1, 2) = val
      end do
    end do

    ! Read trace-free 2-RDM segment
    do j2 = 1, nOb
      do j1 = 1, nOa
        do i2 = 1, nOb
          do i1 = 1, nOa

            count = count + 1

            rdm2ab(i1, i2, j1, j2) = twoRdmState(count)

          end do
        end do
      end do
    end do

    ! Add back trace contribution
    call RdmEnforceTraceM2ab_Apply(rdm2ab, onefoldTraces)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Convert αβ 2-RDM to trace-free representation with explicit 1-RDM.
!>
!> @details
!> Extracts the 1-RDM via partial trace, removes the trace from the 2-RDM,
!> and packs both into the linearized state vector.
!>
!> Includes a sanity check that the resulting trace-free part has zero trace.
!>
!> @param[out] twoRdmState  Output state (pre-allocated)
!> @param[in]  rdm2ab       Input full αβ 2-RDM tensor
  subroutine FillTwoRdmStateFromRdm2ab(twoRdmState, rdm2ab)
    use M_Method_Mb
    use M_Method_Mb_OrbBased
    use M_Utils_RdmTrace
    use M_Utils_RdmEnforceTraceM2ab

    complex(R64), intent(out), contiguous  :: twoRdmState(:)
    complex(R64), intent(in), contiguous   :: rdm2ab(:, :, :, :)

    complex(R64), allocatable :: rdm2abTr0(:, :, :, :)
    complex(R64), allocatable :: rdm1a(:, :)

    integer(I32) :: nOa, nOb

    nOa = Method_Mb_OrbBased_nOrbs(1)
    nOb = Method_Mb_OrbBased_nOrbs(2)

    ! Extract 1-RDM via partial trace
    call RdmTrace_FillRdm1aFromRdm2ab(rdm1a, rdm2ab)

    ! Create trace-free copy
    allocate (rdm2abTr0(nOa, nOb, nOa, nOb))
    rdm2abTr0 = rdm2ab

    call RdmEnforceTraceM2ab_Apply(rdm2abTr0)

    ! Sanity check: verify trace-free property
    block
      complex(R64), allocatable :: testRdm1a(:, :)
      real(R64) :: maxDiff
      call RdmTrace_FillRdm1aFromRdm2ab(testRdm1a, rdm2abTr0)
      maxDiff = maxval(abs(testRdm1a))
      if (maxDiff > 1.0e-12_R64) then
        print *, "WARNING: FillTwoRdmStateFromRdm2ab: rdm2abTr0 is not traceless! maxDiff =", maxDiff
      end if
    end block

    call FillTwoRdmStateFromRdms(twoRdmState, rdm2abTr0_=rdm2abTr0, rdm1a_=rdm1a)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Pack 1-RDM and/or trace-free 2-RDM into the state vector.
!>
!> @details
!> Flexible packing routine that can update either or both segments.
!> If an argument is not present, the corresponding segment is skipped
!> (offset advanced but no data written).
!>
!> @param[out] twoRdmState  Output state vector
!> @param[in]  rdm2abTr0_   Optional: trace-free αβ 2-RDM to pack
!> @param[in]  rdm1a_       Optional: α 1-RDM to pack
  subroutine FillTwoRdmStateFromRdms(twoRdmState, rdm2abTr0_, rdm1a_)
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    complex(R64), intent(out), contiguous  :: twoRdmState(:)
    complex(R64), intent(in), contiguous, optional   :: rdm2abTr0_(:, :, :, :)
    complex(R64), intent(in), contiguous, optional   :: rdm1a_(:, :)

    integer(I32) :: nOa, nOb
    integer(I32) :: count
    integer(I32) :: i1, i2, j1, j2

    nOa = Method_Mb_OrbBased_nOrbs(1)
    nOb = Method_Mb_OrbBased_nOrbs(2)

    count = 0
    if (present(rdm1a_)) then
      do j1 = 1, nOa
        do i1 = 1, nOa
          count = count + 1
          twoRdmState(count) = rdm1a_(i1, j1)
        end do
      end do
    else
      count = count + nOa * nOa
    end if

    if (present(rdm2abTr0_)) then
      do j2 = 1, nOb
        do j1 = 1, nOa
          do i2 = 1, nOb
            do i1 = 1, nOa
              count = count + 1
              twoRdmState(count) = rdm2abTr0_(i1, i2, j1, j2)
            end do
          end do
        end do
      end do
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Extract 1-RDM and/or trace-free 2-RDM from the state vector.
!>
!> @details
!> Flexible extraction routine. Only requested arrays are allocated and filled.
!>
!> @param[in]  twoRdmState  Input state vector
!> @param[out] rdm2abTr0_   Optional: receives trace-free αβ 2-RDM
!> @param[out] rdm1a_       Optional: receives α 1-RDM
  subroutine FillRdmsFromTwoRdmState(twoRdmState, rdm2abTr0_, rdm1a_)
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    complex(R64), intent(in), contiguous   :: twoRdmState(:)
    complex(R64), intent(out), allocatable, optional :: rdm2abTr0_(:, :, :, :)
    complex(R64), intent(out), allocatable, optional :: rdm1a_(:, :)

    integer(I32) :: nOa, nOb
    integer(I32) :: count
    integer(I32) :: i1, i2, j1, j2

    nOa = Method_Mb_OrbBased_nOrbs(1)
    nOb = Method_Mb_OrbBased_nOrbs(2)

    count = 0

    if (present(rdm1a_)) then
      if (.not. allocated(rdm1a_)) allocate (rdm1a_(nOa, nOa))
      do j1 = 1, nOa
        do i1 = 1, nOa
          count = count + 1
          rdm1a_(i1, j1) = twoRdmState(count)
        end do
      end do
    else
      count = count + nOa * nOa
    end if

    if (present(rdm2abTr0_)) then
      if (.not. allocated(rdm2abTr0_)) allocate (rdm2abTr0_(nOa, nOb, nOa, nOb))
      do j2 = 1, nOb
        do j1 = 1, nOa
          do i2 = 1, nOb
            do i1 = 1, nOa
              count = count + 1
              rdm2abTr0_(i1, i2, j1, j2) = twoRdmState(count)
            end do
          end do
        end do
      end do
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Extract only the α 1-RDM from the state vector.
!>
!> @param[out] rdm1a        Allocated on output: α 1-RDM (nα × nα)
!> @param[in]  twoRdmState  Input state vector
  subroutine FillRdm1aFromTwoRdmState(rdm1a, twoRdmState)
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    complex(R64), intent(out), allocatable :: rdm1a(:, :)
    complex(R64), intent(in), contiguous   :: twoRdmState(:)

    integer(I32) :: nOa
    integer(I32) :: count
    integer(I32) :: i1, j1

    nOa = Method_Mb_OrbBased_nOrbs(1)

    if (.not. allocated(rdm1a)) allocate (rdm1a(nOa, nOa))

    count = 0
    do j1 = 1, nOa
      do i1 = 1, nOa
        count = count + 1
        rdm1a(i1, j1) = twoRdmState(count)
      end do
    end do

  end subroutine

end submodule
