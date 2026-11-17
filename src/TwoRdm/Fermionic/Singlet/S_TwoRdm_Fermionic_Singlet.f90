! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Implementation submodule for M_TwoRdm_Fermionic_Singlet.
!>
!> @details
!> Provides concrete implementations for the singlet-restricted fermionic backend.
!> All operations work on the αβ block only and use utility modules to handle
!> singlet symmetry relations for reconstruction of the full 2-RDM.
submodule(M_TwoRdm_Fermionic_Singlet) S_TwoRdm_Fermionic_Singlet

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Fabricate the singlet fermionic backend.
!>
!> @details
!> Sets TwoRdm_nElements to nα²×nβ² and wires procedure pointers to singlet
!> implementations. Also exposes backend-specific ab conversion procedures.
  module subroutine TwoRdm_Fermionic_Singlet_Fabricate
    use M_Utils_Say
    use M_Method_Mb
    use M_Method_Mb_OrbBased
    use M_TwoRdm

    call Say_Fabricate("twoRdm.fermionic.singlet")

    ! The singlet layouts store one (up,down) block: equal orbital and particle
    ! counts of the two spins
    if (Method_Mb_nBodyTypes .ne. 2) error stop "twoRdm.fermionic.singlet requires two body types (the spins)"
    if (Method_Mb_OrbBased_nOrbs(1) .ne. Method_Mb_OrbBased_nOrbs(2) .or. Method_Mb_nBodies(1) .ne. Method_Mb_nBodies(2)) &
      error stop "twoRdm.fermionic.singlet requires equal nOrbs and nBodies for both spins"

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    TwoRdm_nElements = (Method_Mb_OrbBased_nOrbs(1)**2) * (Method_Mb_OrbBased_nOrbs(2)**2)

    TwoRdm_ApplyCommutatorH1 => ApplyCommutatorH1
    TwoRdm_ApplyCommutatorH2 => ApplyCommutatorH2

    TwoRdm_FillRdm1FromTwoRdmState => FillRdm1FromTwoRdmState
    TwoRdm_FillRdm2FromTwoRdmState => FillRdm2FromTwoRdmState
    TwoRdm_FillTwoRdmStateFromRdm2 => FillTwoRdmStateFromRdm2

    TwoRdm_Fermionic_Singlet_FillTwoRdmStateFromRdm2ab => FillTwoRdmStateFromRdm2ab
    TwoRdm_Fermionic_Singlet_FillRdm2abFromTwoRdmState => FillRdm2abFromTwoRdmState

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Accumulate [ĥ₁, D²] for singlet states.
!>
!> @details
!> Computes the one-body commutator in the αβ representation. For singlet states,
!> α and β one-body Hamiltonians are identical in the spatial basis, but are
!> accessed at different index offsets (h1(i,j) for α, h1(i+nα, j+nα) for β).
!>
!> @param[inout] dTwoRdmState  Accumulator for αβ block time derivative
!> @param[in]    twoRdmState   Input linearized αβ 2-RDM
!> @param[in]    h1            Full one-body Hamiltonian (nOrbs × nOrbs spin-orbital basis)
  subroutine ApplyCommutatorH1(dTwoRdmState, twoRdmState, h1)
    use M_Utils_RdmAlgebra

    complex(R64), intent(inout), contiguous  :: dTwoRdmState(:)
    complex(R64), intent(in), contiguous     :: twoRdmState(:)
    complex(R64), intent(in), contiguous     :: h1(:, :)

    complex(R64), allocatable :: rdm2ab(:, :, :, :)
    complex(R64), allocatable :: dRdm2ab(:, :, :, :)

    call TwoRdm_Fermionic_Singlet_FillRdm2abFromTwoRdmState(rdm2ab, twoRdmState)
    call TwoRdm_Fermionic_Singlet_FillRdm2abFromTwoRdmState(dRdm2ab, dTwoRdmState)

    ! [ĥ₁, D²_αβ]: α indices use h1(i,k), β indices use h1(i+nOa, k+nOa)
    call RdmAlgebra_AddCommutatorH1ab(dRdm2ab, rdm2ab, h1)

    call TwoRdm_Fermionic_Singlet_FillTwoRdmStateFromRdm2ab(dTwoRdmState, dRdm2ab)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Accumulate [Ŵ₁₂, D²] for singlet states.
!>
!> @details
!> Computes the two-body commutator in the αβ representation. The interaction
!> tensor h2 is accessed with mixed α/β indices: h2(iα, iβ+nα, kα, kβ+nα).
!>
!> @param[inout] dTwoRdmState  Accumulator for αβ block time derivative
!> @param[in]    twoRdmState   Input linearized αβ 2-RDM
!> @param[in]    h2            Two-body interaction in spin-orbital basis
  subroutine ApplyCommutatorH2(dTwoRdmState, twoRdmState, h2)
    use M_Utils_RdmAlgebra

    complex(R64), intent(inout), contiguous  :: dTwoRdmState(:)
    complex(R64), intent(in), contiguous     :: twoRdmState(:)
    complex(R64), intent(in), contiguous     :: h2(:, :, :, :)

    complex(R64), allocatable :: rdm2ab(:, :, :, :)
    complex(R64), allocatable :: dRdm2ab(:, :, :, :)

    call TwoRdm_Fermionic_Singlet_FillRdm2abFromTwoRdmState(rdm2ab, twoRdmState)
    call TwoRdm_Fermionic_Singlet_FillRdm2abFromTwoRdmState(dRdm2ab, dTwoRdmState)

    call RdmAlgebra_AddCommutatorH2ab(dRdm2ab, rdm2ab, h2)

    call TwoRdm_Fermionic_Singlet_FillTwoRdmStateFromRdm2ab(dTwoRdmState, dRdm2ab)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Compute the 1-RDM from the αβ 2-RDM block.
!>
!> @details
!> Delegates to RdmTrace utility which handles the singlet-specific contraction.
!>
!> @param[out] rdm1         Allocated on output: full 1-RDM (nOrbs × nOrbs)
!> @param[in]  twoRdmState  Input linearized αβ 2-RDM
  subroutine FillRdm1FromTwoRdmState(rdm1, twoRdmState)
    use M_Utils_RdmTrace
    use M_Method_Mb

    complex(R64), intent(out), allocatable :: rdm1(:, :)
    complex(R64), intent(in), contiguous  :: twoRdmState(:)

    complex(R64), allocatable :: rdm2ab(:, :, :, :)

    call TwoRdm_Fermionic_Singlet_FillRdm2abFromTwoRdmState(rdm2ab, twoRdmState)
    call RdmTrace_FillRdm1FromRdm2ab(rdm1, rdm2ab, Method_Mb_nBodies(1), Method_Mb_nBodies(2))

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Reconstruct the full 2-RDM tensor from the αβ block.
!>
!> @details
!> Uses singlet symmetry relations to fill in D²_αα and D²_ββ blocks from D²_αβ.
!> Delegates to RdmSingletRelations utility module.
!>
!> @param[out] rdm2         Allocated on output: full 2-RDM (nOrbs⁴)
!> @param[in]  twoRdmState  Input linearized αβ 2-RDM
  subroutine FillRdm2FromTwoRdmState(rdm2, twoRdmState)
    use M_Utils_RdmSingletRelations

    complex(R64), intent(out), allocatable :: rdm2(:, :, :, :)
    complex(R64), intent(in), contiguous                :: twoRdmState(:)

    complex(R64), allocatable :: rdm2ab(:, :, :, :)

    call TwoRdm_Fermionic_Singlet_FillRdm2abFromTwoRdmState(rdm2ab, twoRdmState)
    call RdmSingletRelations_FillRdm2FromRdm2ab(rdm2, rdm2ab)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Extract and pack the αβ block from a full 2-RDM tensor.
!>
!> @details
!> Converts from full spin-orbital 2-RDM to the singlet αβ representation.
!>
!> @param[out] twoRdmState  Output linearized αβ 2-RDM (pre-allocated)
!> @param[in]  rdm2         Input full 2-RDM tensor
  subroutine FillTwoRdmStateFromRdm2(twoRdmState, rdm2)
    use M_Utils_RdmTrace
    use M_Method_Mb_OrbBased

    complex(R64), intent(out), contiguous  :: twoRdmState(:)
    complex(R64), intent(in), contiguous   :: rdm2(:, :, :, :)

    complex(R64), allocatable :: rdm2ab(:, :, :, :)

    call RdmTrace_FillRdm2abFromRdm2(rdm2ab, rdm2, Method_Mb_OrbBased_nOrbs(1), Method_Mb_OrbBased_nOrbs(2))
    call TwoRdm_Fermionic_Singlet_FillTwoRdmStateFromRdm2ab(twoRdmState, rdm2ab)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Unpack the linearized storage to αβ tensor form.
!>
!> @param[out] rdm2ab       Allocated on output: αβ tensor (nα × nβ × nα × nβ)
!> @param[in]  twoRdmState  Input linearized αβ 2-RDM
  subroutine FillRdm2abFromTwoRdmState(rdm2ab, twoRdmState)
    use M_Method_Mb_OrbBased

    complex(R64), intent(out), allocatable :: rdm2ab(:, :, :, :)
    complex(R64), intent(in), contiguous                :: twoRdmState(:)

    integer(I32) :: nOa, nOb
    integer(I32) :: count
    integer(I32) :: i1, i2, j1, j2

    nOa = Method_Mb_OrbBased_nOrbs(1)
    nOb = Method_Mb_OrbBased_nOrbs(2)

    if (.not. allocated(rdm2ab)) allocate (rdm2ab(nOa, nOb, nOa, nOb))

    count = 0
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

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Pack an αβ tensor into linearized storage.
!>
!> @param[out] twoRdmState  Output linearized αβ 2-RDM (pre-allocated)
!> @param[in]  rdm2ab       Input αβ tensor (nα × nβ × nα × nβ)
  subroutine FillTwoRdmStateFromRdm2ab(twoRdmState, rdm2ab)
    use M_Method_Mb_OrbBased

    complex(R64), intent(out), contiguous  :: twoRdmState(:)
    complex(R64), intent(in), contiguous   :: rdm2ab(:, :, :, :)

    integer(I32) :: nOa, nOb
    integer(I32) :: count
    integer(I32) :: i1, i2, j1, j2

    nOa = Method_Mb_OrbBased_nOrbs(1)
    nOb = Method_Mb_OrbBased_nOrbs(2)

    count = 0
    do j2 = 1, nOb
      do j1 = 1, nOa
        do i2 = 1, nOb
          do i1 = 1, nOa

            count = count + 1

            twoRdmState(count) = rdm2ab(i1, i2, j1, j2)

          end do
        end do
      end do
    end do

  end subroutine

end submodule
