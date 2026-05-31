! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Public absorber interface and procedure pointers.
!>
!> @details
!> Absorbers implement boundary-damping schemes for quantum mechanical
!> simulations. When propagating wave packets on finite numerical grids,
!> reflections from hard boundaries lead to spurious interference. An absorber
!> suppresses these artifacts by applying a smooth mask function M(x) that
!> attenuates the wavefunction near domain edges:
!>
!>   ψ(x) ← M(x) · ψ(x)
!>
!> where M(x) = 1 in the interior (physical region) and M(x) → 0 near the
!> boundaries.
!>
!> This module defines the abstract absorber API and holds the module-level
!> procedure pointers that a concrete implementation installs at runtime via
!> `Absorber_Fabricate`. It does not implement an absorber itself; instead it
!> provides the stable entry points used by the rest of the code base.
!>
!> @par Available Implementations
!> - **Linear/Cosinus**: Cosine-profile mask for 1D linear grids
!>   (see `M_Absorber_Linear_Cosinus`)
!>
!> @par Typical Usage
!> @code{.f90}
!> call Absorber_Fabricate   ! Binds implementation from JSON config
!> call Absorber_Setup       ! Precomputes mask arrays
!> ! ... inside time-stepping loop ...
!> call Absorber_ApplyAbsorber(orbs)  ! Damp orbitals each timestep
!> @endcode
!>
!> @par JSON Configuration
!> If no absorber section is present, a no-op absorber is installed
!> (wavefunctions pass through unchanged).
!>
!> @see M_Absorber_Linear, M_Absorber_Linear_Cosinus
module M_Absorber
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> @brief Selects and installs a concrete absorber implementation.
    !>
    !> @details
    !> This routine reads the JSON configuration (e.g., "absorber.linear.cosinus")
    !> and assigns the module-level procedure pointers `Absorber_Setup` and
    !> `Absorber_ApplyAbsorber` to the selected implementation. If no absorber
    !> block is present in the JSON, a no-op implementation is installed that
    !> leaves wavefunctions unchanged.
    !>
    !> @pre  The JSON configuration must be loaded (via M_Utils_Json).
    !> @post Procedure pointers are bound; call `Absorber_Setup` next.
    module subroutine Absorber_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !---------------------------------------------------------------------------
  !> @brief Pointer to the setup procedure for initializing the absorber.
  !>
  !> @details
  !> After fabrication, call this once to allocate and precompute any
  !> implementation-specific data structures (such as mask arrays or index
  !> ranges). It performs no time stepping; it only prepares state so that
  !> subsequent calls to `Absorber_ApplyAbsorber` are efficient.
  !>
  !> Default: `NoOpProcedures_Setup` (does nothing).
  procedure(I_Absorber_Setup), pointer :: Absorber_Setup => NoOpProcedures_Setup
  abstract interface
    !> @brief Initializes the absorber implementation.
    !>
    !> Typical responsibilities include computing mask functions, validating
    !> configuration, and caching parameters needed by the application phase.
    subroutine I_Absorber_Setup
    end subroutine
  end interface

  !---------------------------------------------------------------------------
  !> @brief Pointer to the procedure for applying the absorber to orbitals.
  !>
  !> @details
  !> Applies the absorption mask function in place to a set of orbitals.
  !> The mask is unity in the interior and smoothly decays to zero near the
  !> boundaries. Conceptually:
  !>
  !>   ψ(x) ← M(x) · ψ(x)
  !>
  !> where M(x) is the implementation-defined mask function.
  !>
  !> @warning This pointer is **not** initialized to a default. It must be
  !>          assigned by `Absorber_Fabricate` before use (either to a real
  !>          implementation or to the no-op fallback).
  procedure(I_Absorber_ApplyAbsorber), pointer :: Absorber_ApplyAbsorber
  abstract interface
    !> @brief Applies the absorption mask function in place.
    !>
    !> @param[in,out] orbs  Orbital matrix with shape (nGrid, nOrbitals).
    !>                      On exit, each orbital has been multiplied
    !>                      element-wise by the absorption mask M(x).
    subroutine I_Absorber_ApplyAbsorber(orbs)
      import :: R64
      complex(R64), intent(inout), contiguous :: orbs(:, :)
    end subroutine
  end interface

end module
