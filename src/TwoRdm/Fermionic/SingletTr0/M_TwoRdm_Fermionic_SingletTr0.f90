! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief SingletTr0 (trace-free singlet) fermionic backend for 2‑RDM operations.
!>
!> @details
!> This backend extends the Singlet backend by explicitly storing the 1-RDM
!> alongside a trace-free αβ 2-RDM block. This parameterization has reduced
!> degrees of freedom and is advantageous for:
!>
!> 1. **Numerical stability**: The 1-RDM is propagated independently, avoiding
!>    accumulation of trace errors during long time propagations.
!> 2. **Reduced redundancy**: The trace constraint Tr_2(D²) ∝ D¹ is built into
!>    the representation, eliminating one source of N-representability violation.
!> 3. **Purification compatibility**: Works naturally with purification schemes
!>    that separately enforce D¹ and D² constraints.
!>
!> ## Storage Layout
!>
!> The linearized state is partitioned into two contiguous segments:
!> ```
!> twoRdmState = [ rdm1a (nα²) | rdm2abTr0 (nα² × nβ²) ]
!>
!> Segment 1 (indices 1 to nα²):
!>   rdm1a(i₁, j₁) at offset = i₁ + nα*(j₁-1)
!>
!> Segment 2 (indices nα²+1 to nα² + nα²×nβ²):
!>   rdm2abTr0(i₁, i₂, j₁, j₂) at offset = nα² + i₁ + nα*(i₂-1) + nα*nβ*(j₁-1) + ...
!> ```
!>
!> ## Trace-Free Constraint
!>
!> The stored D²_αβ satisfies: Σ_k D²_αβ(i,k;j,k) = 0 for all i,j.
!> The full D²_αβ is reconstructed by adding back the 1-RDM contribution.
!>
!> ## Mathematical Background
!>
!> For singlet states, the contraction relation holds:
!>   D¹_α(i,j) = (1/Nβ) · Σ_k D²_αβ(i,k;j,k)
!>
!> We decompose: D²_αβ = D²_αβ,Tr0 + D¹_α ⊗ I_β / Nβ
!> where D²_αβ,Tr0 is the trace-free part stored in this backend.
module M_TwoRdm_Fermionic_SingletTr0
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Activate the singlet fermionic backend and wire procedure pointers.
    !>
    !> Responsibilities:
    !> - Assign `M_TwoRdm` procedure pointers to singlet‑optimized implementations.
    !> - Define the linearization layout consistent with the ab spatial tensor mapping.
    module subroutine TwoRdm_Fermionic_SingletTr0_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! local data
  !=============================================================================

  !> Pointer to the procedure converting a 4D ab spatial‑orbital 2‑RDM to linearized storage.
  procedure(I_TwoRdm_Fermionic_SingletTr0_FillTwoRdmStateFromRdm2ab), pointer &
    :: TwoRdm_Fermionic_SingletTr0_FillTwoRdmStateFromRdm2ab
  abstract interface
    !> Pack the trace-free part of the ab spatial‑orbital 2‑RDM tensor into the backend’s linearized format.
    subroutine I_TwoRdm_Fermionic_SingletTr0_FillTwoRdmStateFromRdm2ab(twoRdmState, rdm2ab)
      import :: R64
      !> Output linearized 2‑RDM (pre‑allocated; size must match backend layout).
      complex(R64), intent(out), contiguous  :: twoRdmState(:)
      !> Input ab 2‑RDM tensor with shape (n, n, n, n) over spatial orbitals.
      complex(R64), intent(in), contiguous   :: rdm2ab(:, :, :, :)
    end subroutine
  end interface

  !> Pointer to the procedure filling the linearized state from trace-free 2-RDM and 1-RDM.
  procedure(I_TwoRdm_Fermionic_SingletTr0_FillTwoRdmStateFromRdms), pointer &
    :: TwoRdm_Fermionic_SingletTr0_FillTwoRdmStateFromRdms
  abstract interface
    !> Pack the trace-free 2-RDM and 1-RDM into the backend’s linearized format.
    subroutine I_TwoRdm_Fermionic_SingletTr0_FillTwoRdmStateFromRdms(twoRdmState, rdm2abTr0_, rdm1a_)
      import :: R64
      !> Output linearized 2‑RDM (pre‑allocated; size must match backend layout).
      complex(R64), intent(out), contiguous  :: twoRdmState(:)
      !> Input trace-free ab 2‑RDM tensor with shape (n, n, n, n).
      complex(R64), intent(in), contiguous, optional   :: rdm2abTr0_(:, :, :, :)
      !> Input 1-RDM tensor with shape (n, n).
      complex(R64), intent(in), contiguous, optional   :: rdm1a_(:, :)
    end subroutine
  end interface

  !> Pointer to the procedure extracting trace-free 2-RDM and 1-RDM from linearized state.
  procedure(I_TwoRdm_Fermionic_SingletTr0_FillRdmsFromTwoRdmState), pointer :: TwoRdm_Fermionic_SingletTr0_FillRdmsFromTwoRdmState
  abstract interface
    !> Unpack the backend’s linearized 2‑RDM into trace-free 2-RDM and 1-RDM.
    subroutine I_TwoRdm_Fermionic_SingletTr0_FillRdmsFromTwoRdmState(twoRdmState, rdm2abTr0_, rdm1a_)
      import :: R64
      !> Input linearized 2‑RDM.
      complex(R64), intent(in), contiguous   :: twoRdmState(:)
      !> Output trace-free ab 2‑RDM tensor with shape (n, n, n, n).
      complex(R64), intent(out), allocatable, optional :: rdm2abTr0_(:, :, :, :)
      !> Output 1-RDM tensor with shape (n, n).
      complex(R64), intent(out), allocatable, optional :: rdm1a_(:, :)
    end subroutine
  end interface

  !> Pointer to the procedure converting the linearized 2‑RDM to a 4D ab spatial tensor.
  procedure(I_TwoRdm_Fermionic_SingletTr0_FillRdm2abFromTwoRdmState), pointer &
    :: TwoRdm_Fermionic_SingletTr0_FillRdm2abFromTwoRdmState
  abstract interface
    !> Unpack the backend’s linearized 2‑RDM into an ab 4‑index spatial‑orbital tensor.
    subroutine I_TwoRdm_Fermionic_SingletTr0_FillRdm2abFromTwoRdmState(rdm2ab, twoRdmState)
      import :: R64
      !> Allocated on output: ab tensor with shape (n, n, n, n) over spatial orbitals.
      complex(R64), intent(out), allocatable :: rdm2ab(:, :, :, :)
      !> Input linearized 2‑RDM.
      complex(R64), intent(in), contiguous   :: twoRdmState(:)
    end subroutine
  end interface

  !> Pointer to the procedure extracting the 1-RDM from the linearized 2-RDM.
  procedure(I_TwoRdm_Fermionic_SingletTr0_FillRdm1aFromTwoRdmState), pointer &
    :: TwoRdm_Fermionic_SingletTr0_FillRdm1aFromTwoRdmState
  abstract interface
    !> Extract the 1-RDM from the backend’s linearized 2‑RDM.
    subroutine I_TwoRdm_Fermionic_SingletTr0_FillRdm1aFromTwoRdmState(rdm1a, twoRdmState)
      import :: R64
      !> Allocated on output: 1-RDM with shape (n, n).
      complex(R64), intent(out), allocatable :: rdm1a(:, :)
      !> Input linearized 2‑RDM.
      complex(R64), intent(in), contiguous   :: twoRdmState(:)
    end subroutine
  end interface

end module
