! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Root interface module for the 2-RDM purification subsystem.
!>
!> @details
!> In time-dependent 2-RDM (TD-2RDM) propagation, numerical errors and
!> approximate BBGKY closures can cause the density matrix to violate
!> N-representability conditions—specifically, the two-particle RDM (D²) or
!> its hole counterpart (Q-RDM) may acquire negative eigenvalues.
!>
!> This module provides the runtime-dispatched entry point for purification
!> algorithms that iteratively restore positivity constraints:
!>
!>   D² ≥ 0   (2-RDM positive semidefiniteness)
!>   Q² ≥ 0   (2-hole-RDM positive semidefiniteness)
!>
!> The concrete implementation is selected at runtime via `Purification_Fabricate`,
!> which reads the JSON configuration tree under `purification.*` and binds the
!> appropriate backend. If no purification key is present, the pointers remain
!> bound to no-op procedures and purification is effectively bypassed.
!>
!> @par Available Backends
!> | JSON Key                        | Backend Module                       | Description                                     |
!> |---------------------------------|--------------------------------------|-------------------------------------------------|
!> | `purification.singlet.reduce`   | `M_Purification_Singlet_Reduce`      | Singlet spin-symmetry with iterative reduction  |
!> | `purification.singletTr0.reduce`| `M_Purification_SingletTr0_Reduce`   | SingletTr0 (trace-free) with iterative reduction|
!> | `purification.generic.alcoba`   | `M_Purification_Generic_Alcoba`      | Spin-agnostic Alcoba/Valdemoro-style (stub)     |
!>
!> @par Algorithm Overview (Reduce Method)
!> The "reduce" backends implement the following iterative scheme:
!> 1. Unpack the 2-RDM from the compressed state vector.
!> 2. Split into spin-adapted blocks (singlet/triplet).
!> 3. Compute defective components via eigendecomposition (negative eigenspace).
!> 4. Subtract trace-free projections of defects from both D² and Q².
!> 5. Optionally project onto the energy-conserving subspace.
!> 6. Repeat until min(eigenvalues) > threshold or max iterations reached.
!> 7. Repack into the state vector.
!>
!> @par Usage Contract
!> @code{.f90}
!> call Purification_Fabricate     ! Once during initialization
!> call Purification_Setup         ! Optional backend-specific setup
!> ...
!> call Purification_PurifyTwoRdmState(twoRdmState, h2)  ! During propagation
!> @endcode
!>
!> @par Thread Safety
!> The procedure pointers are module-level and thus not thread-safe.
!>
!> @see M_TwoRdm for the 2-RDM state container
!> @see M_Reconstruction for the BBGKY closure approximations
!> @see M_Utils_RdmHole for Q-RDM (hole-RDM) computation
module M_Purification
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> @brief Dynamically bind purification procedures based on JSON configuration.
    !>
    !> @details
    !> Inspects the configuration tree under `purification.*` and delegates to
    !> the appropriate branch fabricator:
    !>
    !> | JSON Path               | Delegates To                         |
    !> |-------------------------|--------------------------------------|
    !> | `purification.singlet`  | `Purification_Singlet_Fabricate`     |
    !> | `purification.singletTr0`| `Purification_SingletTr0_Fabricate` |
    !> | `purification.generic`  | `Purification_Generic_Fabricate`     |
    !>
    !> If none of these keys exist, purification is disabled and a diagnostic
    !> message is emitted. The procedure pointers remain bound to no-ops.
    !>
    !> @note Must be called after JSON configuration is loaded but before
    !>       any call to `Purification_PurifyTwoRdmState`.
    module subroutine Purification_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  ! (No module-level data; all state is held by submodules.)

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> @brief Setup procedure pointer for the selected purification backend.
  !> @details
  !> Called once after fabrication to perform any backend-specific initialization
  !> (e.g., allocating workspace, precomputing invariants). Most backends leave
  !> this as a no-op. Default: `NoOpProcedures_Setup`.
  procedure(I_Purification_Setup), pointer :: Purification_Setup => NoOpProcedures_Setup
  abstract interface
    subroutine I_Purification_Setup
    end subroutine
  end interface

  !> @brief Main entry point for 2-RDM purification.
  !> @details
  !> Applies the configured purification algorithm in-place to enforce
  !> N-representability conditions on the packed 2-RDM state vector.
  !> Default: `NoOpProcedures_PurifyTwoRdmState` (identity operation).
  !>
  !> @par Physical Background
  !> The 2-RDM must satisfy ensemble N-representability conditions, among which
  !> the D- and Q-conditions require:
  !>   - D²(ijkl) ≥ 0  (particle-particle RDM is positive semidefinite)
  !>   - Q²(ijkl) ≥ 0  (hole-hole RDM is positive semidefinite)
  !>
  !> When these are violated due to numerical error or approximate closures,
  !> the purification algorithm iteratively projects the defective (negative)
  !> eigenspace out of both D² and Q².
  procedure(I_Purification_PurifyTwoRdmState), pointer :: Purification_PurifyTwoRdmState => NoOpProcedures_PurifyTwoRdmState
  abstract interface
    !> @brief Apply purification to a packed 2-RDM state vector.
    !>
    !> @param[inout] twoRdmState  Packed 2-RDM buffer (layout depends on backend:
    !>                            Singlet uses α-β block, SingletTr0 adds trace-free
    !>                            constraint). Modified in-place to restore positivity.
    !> @param[in]    h2           Two-body Hamiltonian tensor h²(p,q,r,s). Used only
    !>                            when `conserveEnergyQ=.TRUE.` to project corrections
    !>                            onto the energy-conserving subspace (orthogonal to h²).
    subroutine I_Purification_PurifyTwoRdmState(twoRdmState, h2)
      import :: R64
      complex(R64), intent(inout), contiguous :: twoRdmState(:)
      complex(R64), intent(in), contiguous    :: h2(:, :, :, :)
    end subroutine
  end interface

end module

