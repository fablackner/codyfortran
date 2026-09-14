! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Public runtime API for the two‑body reduced density matrix (2‑RDM) subsystem.
!>
!> @details
!> The two-body reduced density matrix (2-RDM) encapsulates two-particle correlations
!> in a many-body quantum system. Given an N-particle wavefunction |Ψ⟩, the 2-RDM is
!> defined as:
!>
!>   D²(i₁,i₂;j₁,j₂) = ⟨Ψ| a†_{j₁} a†_{j₂} a_{i₂} a_{i₁} |Ψ⟩
!>
!> where a†/a are fermionic creation/annihilation operators.
!>
!> This module does not contain concrete implementations. Instead, it exposes
!> shared state and a set of procedure pointers that are wired up at runtime by
!> `TwoRdm_Fabricate` to a concrete backend (e.g. fermionic generic, singlet,
!> or singletTr0).
!>
!> ## Architecture
!>
!> The module follows the Interface/Implementation separation pattern:
!> - **M_TwoRdm** (this file): Declares procedure pointer interfaces and shared state
!> - **S_TwoRdm** (submodule): Implements fabrication logic
!> - **Fermionic backends**: Provide concrete implementations for different symmetry classes
!>
!> ## Key Points for Consumers
!>
!> - `TwoRdm_Setup` allocates/associates the 2‑RDM portion of the global state
!>   and initializes module bookkeeping. It is a no‑op until `TwoRdm_Fabricate`
!>   assigns a concrete implementation.
!> - `TwoRdm_twoRdmState(:)` is the linearized storage used by all backends; its
!>   serialization layout is backend‑defined but consistent across calls once
!>   fabricated.
!> - The commutator routines (`ApplyCommutatorH1`, `ApplyCommutatorH2`) accumulate
!>   results into the provided output array and never allocate.
!> - Conversion routines (`FillRdm1FromTwoRdmState`, `FillRdm2FromTwoRdmState`)
!>   allocate their outputs.
!>
!> ## Memory Layout
!>
!> The linearized 2-RDM storage `twoRdmState(:)` packing depends on the backend:
!> - **Generic**: Full nOrbs⁴ tensor in column-major order (i₁ fastest)
!> - **Singlet**: αβ block only: nα² × nβ² elements
!> - **SingletTr0**: 1-RDM (nα²) + trace-free αβ block (nα² × nβ²)
!>
!> ## Typical Usage
!>
!> ```fortran
!> call TwoRdm_Fabricate()    ! Wire procedure pointers from JSON config
!> call TwoRdm_Setup()        ! Allocate storage
!> ! ... populate TwoRdm_twoRdmState ...
!> call TwoRdm_ApplyCommutatorH1(dState, state, h1)  ! Time evolution
!> call TwoRdm_FillRdm1FromTwoRdmState(rdm1, state)  ! Extract 1-RDM
!> ```
module M_TwoRdm
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Selects and wires a concrete 2‑RDM backend at runtime.
    !>
    !> Responsibilities:
    !> - Choose a concrete representation (e.g. fermionic generic/singlet).
    !> - Assign all procedure pointers in this module to backend procedures.
    !> - Optionally perform early validation of configuration relevant to 2‑RDMs.
    module subroutine TwoRdm_Fabricate()
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Total number of elements in the linearized 2‑RDM representation currently in use.
  !> Set during setup/fabrication and consistent with `TwoRdm_twoRdmState`.
  integer(I32) :: TwoRdm_nElements

  !> Pointer to the 2‑RDM portion of the global/state vector in linearized form.
  !> - Association and storage are backend‑controlled (usually during `TwoRdm_Setup`).
  !> - When associated, its size equals `TwoRdm_nElements`.
  !> - The exact packing order is backend‑defined but stable for a given backend.
  complex(R64), contiguous, pointer :: TwoRdm_twoRdmState(:)

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> Pointer to the setup procedure for initializing the 2‑RDM subsystem.
  procedure(I_TwoRdm_Setup), pointer :: TwoRdm_Setup => NoOpProcedures_Setup
  abstract interface
    !> Initialize and/or associate 2‑RDM storage and bookkeeping.
    !>
    !> Expected effects:
    !> - Set `TwoRdm_nElements`.
    !> - Allocate/associate `TwoRdm_twoRdmState(:)` if needed.
    !> - Prepare any backend‑specific constants/indices used by other procedures.
    subroutine I_TwoRdm_Setup
    end subroutine
  end interface

  !> Pointer to the procedure that accumulates the [H1, D(2)] contribution.
  procedure(I_TwoRdm_ApplyCommutatorH1), pointer :: TwoRdm_ApplyCommutatorH1
  abstract interface
    !> Accumulate the one‑body Hamiltonian commutator contribution in linearized form.
    !>
    !> Contract a one‑body operator H1 with the current 2‑RDM and add the result to
    !> `dTwoRdmState(:)`. The concrete index algebra and serialization are backend‑defined.
    subroutine I_TwoRdm_ApplyCommutatorH1(dTwoRdmState, twoRdmState, h1)
      import :: R64
      !> In/out accumulator for d/dt of the linearized 2‑RDM (pre‑initialized by caller).
      complex(R64), intent(inout), contiguous :: dTwoRdmState(:)
      !> Input linearized 2‑RDM.
      complex(R64), intent(in), contiguous :: twoRdmState(:)
      !> One‑body operator matrix H1 with shape (nOrb, nOrb).
      complex(R64), intent(in), contiguous :: h1(:, :)
    end subroutine
  end interface

  !> Pointer to the procedure that accumulates the [H2, D(2)] contribution.
  procedure(I_TwoRdm_ApplyCommutatorH2), pointer :: TwoRdm_ApplyCommutatorH2
  abstract interface
    !> Accumulate the two‑body Hamiltonian commutator contribution in linearized form.
    !>
    !> Contract a two‑body operator H2 with the current 2‑RDM and add the result to
    !> `dTwoRdmState(:)`. The concrete index algebra and serialization are backend‑defined.
    subroutine I_TwoRdm_ApplyCommutatorH2(dTwoRdmState, twoRdmState, h2)
      import :: R64
      !> In/out accumulator for d/dt of the linearized 2‑RDM (pre‑initialized by caller).
      complex(R64), intent(inout), contiguous :: dTwoRdmState(:)
      !> Input linearized 2‑RDM.
      complex(R64), intent(in), contiguous :: twoRdmState(:)
      !> Two‑body operator tensor H2 with shape (nOrb, nOrb, nOrb, nOrb).
      complex(R64), intent(in), contiguous :: h2(:, :, :, :)
    end subroutine
  end interface

  !> Pointer to the procedure that computes the 1‑RDM from the linearized 2‑RDM.
  procedure(I_TwoRdm_FillRdm1FromTwoRdmState), pointer :: TwoRdm_FillRdm1FromTwoRdmState
  abstract interface
    !> Build the one‑body reduced density matrix from the current linearized 2‑RDM.
    !>
    !> The contraction and normalization follow the backend’s conventions.
    subroutine I_TwoRdm_FillRdm1FromTwoRdmState(rdm1, twoRdmState)
      import :: R64
      !> Allocated on output: 1‑RDM with shape (nOrb, nOrb).
      complex(R64), intent(out), allocatable :: rdm1(:, :)
      !> Input linearized 2‑RDM.
      complex(R64), intent(in), contiguous :: twoRdmState(:)
    end subroutine
  end interface

  !> Pointer to the procedure that unpacks the linearized 2‑RDM to a 4‑index tensor.
  procedure(I_TwoRdm_FillRdm2FromTwoRdmState), pointer :: TwoRdm_FillRdm2FromTwoRdmState
  abstract interface
    !> Reconstruct the full 2‑RDM tensor from the linearized storage used by the backend.
    !>
    !> Index/order conventions are backend‑defined but documented by the selected backend.
    subroutine I_TwoRdm_FillRdm2FromTwoRdmState(rdm2, twoRdmState)
      import :: R64
      !> Allocated on output: 2‑RDM with shape (nOrb, nOrb, nOrb, nOrb).
      complex(R64), intent(out), allocatable :: rdm2(:, :, :, :)
      !> Input linearized 2‑RDM.
      complex(R64), intent(in), contiguous :: twoRdmState(:)
    end subroutine
  end interface

  !> Pointer to the procedure that packs a 4‑index 2‑RDM tensor into linearized storage.
  procedure(I_TwoRdm_FillTwoRdmStateFromRdm2), pointer :: TwoRdm_FillTwoRdmStateFromRdm2
  abstract interface
    !> Convert a 4‑index 2‑RDM tensor into the backend’s linearized representation.
    !>
    !> The caller must ensure the provided tensor follows the backend’s symmetry
    !> and ordering conventions; no allocation is performed here.
    subroutine I_TwoRdm_FillTwoRdmStateFromRdm2(twoRdmState, rdm2)
      import :: R64
      !> Output linearized 2‑RDM (pre‑allocated, size = `TwoRdm_nElements`).
      complex(R64), intent(out), contiguous :: twoRdmState(:)
      !> Input full 2‑RDM tensor with shape (nOrb, nOrb, nOrb, nOrb).
      complex(R64), intent(in), contiguous :: rdm2(:, :, :, :)
    end subroutine
  end interface

  !> Pointer to the procedure that enforces backend‑defined normalization/consistency.
  procedure(I_TwoRdm_Normalize), pointer :: TwoRdm_Normalize
  abstract interface
    !> Bring the linearized 2‑RDM into a normalized/consistent form required by the backend.
    subroutine I_TwoRdm_Normalize(twoRdmState)
      import :: R64
      !> In/out linearized 2‑RDM. Modified in place.
      complex(R64), intent(inout), contiguous :: twoRdmState(:)
    end subroutine
  end interface

  !> Pointer to the procedure that persists the linearized 2‑RDM to disk or a stream.
  procedure(I_TwoRdm_SaveTwoRdmState), pointer :: TwoRdm_SaveTwoRdmState
  abstract interface
    !> Persist the 2‑RDM in its linearized form. The concrete format/location is backend‑ or
    !> application‑defined.
    subroutine I_TwoRdm_SaveTwoRdmState(twoRdmState)
      import :: R64
      !> Linearized 2‑RDM to be saved.
      complex(R64), intent(in), contiguous :: twoRdmState(:)
    end subroutine
  end interface

end module

