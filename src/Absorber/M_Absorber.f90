! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Public absorber interface and procedure pointers.
!>
!> This module defines the abstract absorber API and holds the module-level
!> procedure pointers that a concrete implementation installs at runtime via
!> `Absorber_Fabricate`. It does not implement an absorber itself; instead it
!> provides the stable entry points used by the rest of the code base.
module M_Absorber
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Selects and installs a concrete absorber implementation and reads input.
    !>
    !> This routine assigns the module-level procedure pointers below to a
    !> specific implementation (for example a linear-grid absorber with a
    !> cosine profile) and ingests user input from the project configuration
    !> (e.g., JSON). Call this once during initialization before using
    !> `Absorber_Setup` or `Absorber_ApplyAbsorber`.
    module subroutine Absorber_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> Pointer to the setup procedure for initializing the absorber system.
  !>
  !> After fabrication, call this once to allocate and/or precompute any
  !> implementation-specific data structures (such as mask arrays or index
  !> ranges). It performs no time stepping; it only prepares state so that
  !> subsequent calls to `Absorber_ApplyAbsorber` are fast.
  procedure(I_Absorber_Setup), pointer :: Absorber_Setup => NoOpProcedures_Setup
  abstract interface
    !> Initializes the absorber implementation.
    !>
    !> Typical responsibilities include computing mask functions, validating
    !> configuration, and caching parameters needed by the application phase.
    subroutine I_Absorber_Setup
    end subroutine
  end interface

  !> Pointer to the procedure for applying the absorber to orbitals.
  procedure(I_Absorber_ApplyAbsorber), pointer :: Absorber_ApplyAbsorber
  abstract interface
    !> Applies the absorption mask function in place to a set of orbitals.
    !>
    !> The mask is generally unity in the interior and smoothly decays to zero
    !> near the boundaries. Conceptually: $\psi(x) \leftarrow M(x)\,\psi(x)$,
    !> where $M(x)$ is the implementation-defined mask function.
    subroutine I_Absorber_ApplyAbsorber(orbs)
      import :: R64
      !> Input/output orbitals with shape (nGrid, nOrbitals). On exit, the
      !> array has been multiplied element-wise by the absorption mask.
      complex(R64), intent(inout), contiguous :: orbs(:, :)
    end subroutine
  end interface

end module
