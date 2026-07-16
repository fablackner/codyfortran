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
!> The module follows the CodyFortranRDM fabrication pattern:
!>   - **Interface Module** (`M_Mixing.f90`): Declares the procedure pointer
!>   - **Dispatch Submodule** (`S_Mixing.f90`): Routes JSON to backends
!>   - **Implementation Submodules**: Provide concrete algorithms
!>
!> Mixer state (dimension, history, damping parameter α) is held as module
!> variables inside the implementation submodules — there is one global mixer
!> active at a time.
!>
!> ## Usage
!>
!> 1. Call `Mixing_Fabricate()` to bind the JSON-selected backend
!>    (called by the program alongside the other module fabrications)
!> 2. Per iteration: `call Mixing_Mix(x, xNew)`
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
!>       product; they carry no Grid dependency. Consumers must remove any
!>       gauge freedom from xNew before mixing (see
!>       `Orbs_AlignOnReference`).
!>
!> @see M_GroundSolver for the main consumer
module M_Mixing
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface Mixing_Fabricate
    !> @brief Bind the JSON-selected mixing backend to the procedure pointer.
    !>
    !> @details
    !> Reads the JSON configuration key `mixing.*` and dispatches to the
    !> appropriate backend's `_Fabricate` routine. If no `mixing` block is
    !> present, the linear backend is selected (backward-compatible default).
    !>
    !> @pre JSON configuration loaded via `Json_LoadJsonFile`
    !> @post `Mixing_Mix` is bound
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

end module
