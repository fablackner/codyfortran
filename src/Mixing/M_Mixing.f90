! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Iteration-mixing (convergence acceleration) interface module.
!>
!> @details
!> Abstracts self-consistency mixers: given the currently accepted iterate x
!> and a new candidate xNew produced by one fixed-point step (e.g., one SCF
!> cycle), a mixer computes the next accepted iterate. The mixer is agnostic
!> to what the vector represents — consumers decide whether to feed orbitals,
!> a density/potential, or any other flattened complex quantity.
!>
!> ## Architecture
!>
!> The module follows the CodyFortran fabrication pattern:
!>   - **Interface Module** (`M_Mixing.f90`): Declares the procedure pointers
!>   - **Dispatch Submodule** (`S_Mixing.f90`): Routes JSON to backends
!>   - **Implementation Submodules**: Provide concrete algorithms
!>
!> There is one mixer per run: its state (dimension, history) lives in the
!> implementation submodule. `Mixing_Fabricate` starts it empty;
!> `Mixing_Reset` empties it when the fixed-point map changes (e.g. a new
!> functional or another system) while the iteration continues.
!>
!> ## Usage
!>
!> 1. Call `Mixing_Fabricate()` to bind the JSON-selected backend
!>    (called by the program alongside the other module fabrications)
!> 2. Per iteration: `call Mixing_Mix(x, xNew)`
!> 3. Before iterating another fixed-point map: `call Mixing_Reset()`
!>
!> ## Available Backends
!>
!> - **linear**: x ← (1−α)·x + α·xNew (stateless; the default when the JSON
!>   has no `mixing` block)
!> - **diis**: Pulay/Anderson mixing — extrapolates over a history of
!>   (iterate, residual) pairs and takes a damped step of size α along the
!>   extrapolated residual; with an empty history it reduces to linear mixing
!>
!> @note Mixers operate on raw complex vectors with the Euclidean inner
!>       product; they carry no Grid dependency (the grid metric does not
!>       speed up DIIS on orbitals: A_Mixing.md). Consumers must remove any
!>       gauge freedom from xNew before mixing (see `Orbs_AlignOnReference`).
!>
!> @see M_GroundSolver for the main consumer
module M_Mixing
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface Mixing_Fabricate
    !> @brief Bind the JSON-selected mixing backend to the procedure pointers.
    !>
    !> @details
    !> Reads the JSON configuration key `mixing.*` and dispatches to the
    !> appropriate backend's `_Fabricate` routine. If no `mixing` block is
    !> present, the linear backend is selected (backward-compatible default).
    !> The mixer starts empty (also when fabricated again).
    !>
    !> @pre JSON configuration loaded via `Json_LoadJsonFile`
    !> @post `Mixing_Mix` and `Mixing_Reset` are bound
    module subroutine Mixing_Fabricate()
    end subroutine
  end interface

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> @brief Pointer to the mixing step of the selected backend.
  procedure(I_Mixing_Mix), pointer :: Mixing_Mix => null()
  abstract interface
    !> @brief Compute the next accepted iterate from (x, xNew).
    !>
    !> @details
    !> Contract: xMixed holds the currently accepted iterate on entry and the next
    !> accepted iterate on exit; xRaw is the candidate produced by one
    !> fixed-point step applied to xMixed. The damping parameter α is stored
    !> internally in the backend.
    !>
    !> If the vector dimension changes between calls, the backend
    !> (re)allocates its internal state automatically.
    !>
    !> @param[inout] xMixed  The quantity being iterated (e.g., orbitals, potential).
    !>                       On input: the current accepted iteration.
    !>                       On output: the new mixed (damped/extrapolated) step.
    !> @param[in]    xRaw    The raw, unmixed new iteration produced by the fixed-point map.
    subroutine I_Mixing_Mix(xMixed, xRaw)
      import :: R64
      !> The quantity being iterated, updated in place.
      complex(R64), intent(inout), contiguous :: xMixed(:)
      !> The raw new iteration produced by the fixed-point map.
      complex(R64), intent(in), contiguous :: xRaw(:)
    end subroutine
  end interface

  !> @brief Pointer to the reset of the selected backend.
  procedure(I_Mixing_Reset), pointer :: Mixing_Reset => null()
  abstract interface
    !> @brief Drop the mixer history (DIIS: stored pairs and start gate).
    !>
    !> @details
    !> Call it when the fixed-point map changes but the iteration continues
    !> from the current iterate (another functional, another system): the
    !> stored (iterate, residual) pairs belong to the old map and would
    !> misguide the extrapolation. No-op for the stateless linear backend.
    subroutine I_Mixing_Reset()
    end subroutine
  end interface

end module
