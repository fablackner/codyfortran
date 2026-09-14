! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file M_TwoRdmInit.f90
!> @brief High-level entry point for initializing two-body reduced density
!>        matrices (2-RDM).
!>
!> @details
!> The two-body reduced density matrix D²(i,j;k,l) encodes all two-particle
!> correlations in a quantum many-body system. It is defined as:
!>
!>     D²(i,j;k,l) = ⟨Ψ| a†_i a†_j a_l a_k |Ψ⟩
!>
!> This module provides a **factory pattern** interface for initializing
!> 2-RDMs from various sources:
!>
!>   - **Excited**: Build D² from 1-RDM products (mean-field approximation
!>                  with fermionic exchange for same-spin blocks)
!>   - **Load**: Read D² from binary files (restart/continuation runs)
!>
!> After fabrication, procedure pointers are bound to concrete implementations
!> allowing runtime polymorphism without modifying client code.
!>
!> @par Architecture
!> Follows the Interface Module + Implementation Submodule pattern:
!>   - M_TwoRdmInit.f90       : Public API (procedure pointers)
!>   - S_TwoRdmInit.f90       : Factory dispatcher
!>   - Excited/, Load/        : Concrete implementations
!>
!> @par Usage Sequence
!> @code{.f90}
!> call TwoRdmInit_Fabricate      ! Reads JSON, binds implementation
!> call TwoRdmInit_Setup          ! One-time setup (currently no-op)
!> call TwoRdmInit_Initialize(twoRdmState)  ! Produce full 2-RDM state
!> @endcode
!>
!> @see AGENTS.md for onboarding documentation and extension guide.
!> @see M_TwoRdm for 2-RDM state vector layout and operations.
module M_TwoRdmInit
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> @brief Bind procedure pointers to a concrete 2-RDM initializer
    !>        based on JSON configuration.
    !>
    !> @details
    !> Reads the `"twoRdmInit"` block from the active JSON configuration
    !> and dispatches to one of:
    !>   - TwoRdmInit_Load_Fabricate    if `"twoRdmInit.load"` exists
    !>   - TwoRdmInit_Excited_Fabricate if `"twoRdmInit.excited"` exists
    !>
    !> @par Expected Effects
    !>   - Binds `TwoRdmInit_Setup` (currently no-op for all implementations)
    !>   - Binds `TwoRdmInit_Initialize` to the unified block-assembly routine
    !>   - Binds `TwoRdmInit_InitializeBlock` to implementation-specific logic
    !>   - Caches configuration parameters (e.g., creates/destroys indices)
    !>
    !> @par JSON Configuration
    !> @code{.json}
    !> { "twoRdmInit": { "excited": { ... } } }
    !> { "twoRdmInit": { "load": { } } }
    !> @endcode
    !>
    !> @note Exactly one of `"load"` or `"excited"` must be specified.
    !>       The routine will `error stop` otherwise.
    module subroutine TwoRdmInit_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> @brief Pointer to the setup procedure for the selected 2-RDM initializer.
  !>
  !> @details
  !> Call this once after fabrication and before any `TwoRdmInit_Initialize`
  !> calls. Implementations may:
  !>   - Validate array dimensions against system configuration
  !>   - Precompute lookup tables or index mappings
  !>   - Allocate persistent scratch buffers
  !>
  !> Currently all implementations use the default no-op; setup logic is
  !> embedded in `TwoRdmInit_Fabricate`. This pointer exists for consistency
  !> with the framework's Fabricate -> Setup -> Use lifecycle.
  !>
  !> @note Default: `NoOpProcedures_Setup` (does nothing).
  procedure(I_TwoRdmInit_Setup), pointer :: TwoRdmInit_Setup => NoOpProcedures_Setup
  abstract interface
    !> @brief Setup hook interface for 2-RDM initializer implementations.
    subroutine I_TwoRdmInit_Setup
    end subroutine
  end interface

  !> @brief Pointer to the procedure that initializes a full 2-RDM state vector.
  !>
  !> @details
  !> This is the main entry point for clients needing a complete 2-RDM.
  !> The routine internally:
  !>   1. Iterates over all body-type pairs (bt1, bt2) in {1,2} x {1,2}
  !>   2. Calls `TwoRdmInit_InitializeBlock` for each block
  !>   3. Assembles blocks into a full 4-index tensor
  !>   4. Packs the tensor into the linearized state vector format
  !>
  !> @par Contract
  !>   - `twoRdmState(:)` must have length `TwoRdm_nTwoRdmState` (from M_TwoRdm)
  !>   - The array is fully overwritten on return
  !>
  !> @par Body-Type Block Structure (2-species fermionic system)
  !>   - (1,1): alpha-alpha block (antisymmetric under exchange)
  !>   - (1,2): alpha-beta block (no exchange symmetry)
  !>   - (2,1): beta-alpha block (no exchange symmetry)
  !>   - (2,2): beta-beta block (antisymmetric under exchange)
  !>
  !> @note Not associated until `TwoRdmInit_Fabricate` is called.
  procedure(I_TwoRdmInit_Initialize), pointer :: TwoRdmInit_Initialize
  abstract interface
    !> @brief Produce a fully-initialized 2-RDM state vector.
    !> @param[out] twoRdmState  Linearized 2-RDM state (length TwoRdm_nTwoRdmState)
    subroutine I_TwoRdmInit_Initialize(twoRdmState)
      import R64
      complex(R64), intent(out), contiguous :: twoRdmState(:)
    end subroutine
  end interface

  !> @brief Pointer to the procedure that initializes a single 2-RDM block.
  !>
  !> @details
  !> Produces the 2-RDM tensor block D^2_{bt1,bt2}(i1,i2;j1,j2) where:
  !>   - i1, j1 index orbitals of body type bt1
  !>   - i2, j2 index orbitals of body type bt2
  !>
  !> @par Tensor Layout
  !> @code
  !> rdm2Bt(i1, i2, j1, j2) = <Psi| a^dag_{i1,bt1} a^dag_{i2,bt2} a_{j2,bt2} a_{j1,bt1} |Psi>
  !> @endcode
  !>
  !> @par Symmetry Properties
  !>   - **Same-spin (bt1 == bt2)**: Antisymmetric under (i1<->i2, j1<->j2)
  !>     due to fermionic exchange: D^2(i,j;k,l) = -D^2(j,i;l,k)
  !>   - **Opposite-spin (bt1 /= bt2)**: No exchange symmetry requirement
  !>
  !> @par Contract
  !>   - `rdm2Bt` shape: (nOrbs_bt1, nOrbs_bt2, nOrbs_bt1, nOrbs_bt2)
  !>   - Array is fully overwritten on return
  !>   - Other blocks remain unchanged
  !>
  !> @note Not associated until `TwoRdmInit_Fabricate` is called.
  procedure(I_TwoRdmInit_InitializeBlock), pointer :: TwoRdmInit_InitializeBlock
  abstract interface
    !> @brief Produce a single 2-RDM tensor block for a body-type pair.
    !> @param[out] rdm2Bt  4-index tensor block (i1, i2, j1, j2)
    !> @param[in]  bt1     Body type index for first particle (1 or 2)
    !> @param[in]  bt2     Body type index for second particle (1 or 2)
    subroutine I_TwoRdmInit_InitializeBlock(rdm2Bt, bt1, bt2)
      import I32, R64
      complex(R64), intent(out), contiguous :: rdm2Bt(:, :, :, :)
      integer(I32), intent(in)               :: bt1
      integer(I32), intent(in)               :: bt2
    end subroutine
  end interface

end module
