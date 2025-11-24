! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Core configuration abstraction used across the project.
!>
!> This module defines the abstract base type that all configuration list
!> implementations derive from (e.g. all-active bosons/fermions). It
!> establishes:
!> - a uniform interface for setup, fabrication (construction from input), and
!>   applying excitations to configurations
!> - common storage for excitation connectivity (singles/doubles)
!> - a thin container type to build heterogeneous arrays of configuration
!>   providers.
!>
!> The concrete specializations (statistics, truncation scheme, etc.) live in
!> separate modules and extend `T_ConfigList_E`.
module M_ConfigList
  use M_Utils_Types
  use M_Utils_NoOpProcedures, only: NoOpProcedures_Setup

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Populate the global configuration list.
    !>
    !> This routine is responsible for creating/allocating the concrete
    !> configuration-list elements (bosonic/fermionic, truncation scheme, ...)
    !> based on the parsed input and registering them in `configList`.
    !> The actual implementation resides in a higher layer of the build
    !> pipeline; this module only declares the interface.
    module subroutine ConfigList_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module types
  !=============================================================================

  type :: T_ConfigList_Excitation
    !> Per-configuration connectivity for a given excitation rank.
    !>
    !> The fields are structured to support fast matrix construction and
    !> application of excitation operators:
    !> - `nConnected(iC)` gives the number of connected configurations for
    !>   configuration index `iC` at this excitation rank
    !> - `excitedC(:, iC)` stores the target configuration indices
    !> - `factor(:, iC)` holds the scalar factors to apply when connecting
    !>   `iC -> excitedC(:, iC)` (includes statistics-dependent phase for
    !>   fermions); the element type is complex to support complex-valued
    !>   Hamiltonians
    !> - `orbCode(:, iC)` encodes the orbital indices that define the
    !>   excitation for that edge
    integer(I32), allocatable :: nConnected(:)
    !> Target configurations for each connection edge.
    integer(I32), allocatable :: excitedC(:, :)
    !> Factors for each connection (includes fermionic sign when applicable).
    complex(R64), allocatable :: factor(:, :)
    !> Encoded orbital indices describing the excitation pattern.
    integer(I32), allocatable :: orbCode(:, :)
  end type

  !> Abstract base class for all configuration list implementations.
  !>
  !> Concrete specializations provide statistics-specific behavior and
  !> implement the deferred bindings below. The base defines minimal shared
  !> state and the standard call surface used by the rest of the code base.
  type, abstract :: T_ConfigList_E

    !> Which body this element targets (index into a global list of bodies).
    integer(I32)     :: bodyTarget

    !> Total number of configurations represented by this element.
    integer(I32)     :: nConfigurations

    !> Optional mapping from linear configuration index to an internal code.
    !> Concrete implementations define the coding scheme (bit patterns,
    !> occupation vectors, etc.).
    integer(I64), allocatable :: codeFromConfig(:)

    !> Precomputed single and double excitation connectivity. Implementations
    !> decide when/how these are filled (during fabrication or lazy at setup).
    type(T_ConfigList_Excitation) :: singles
    type(T_ConfigList_Excitation) :: doubles

    !> Hierarchical path in the input (e.g., JSON/YAML) this element was
    !> constructed from. Primarily for traceability and diagnostics.
    character(len=:), allocatable :: path

  contains

    !> Perform runtime setup. Called after fabrication when all elements exist
    !> and global context is available. Typical tasks: allocate working buffers,
    !> build caches, finalize connectivity.
    procedure(I_Setup), deferred :: Setup

    !> Construct element data from the parsed input and global settings. Called
    !> during the build phase before `Setup`.
    procedure(I_Fabricate), deferred :: Fabricate

    !> Apply a list of creation/annihilation operators to configuration `iC` and
    !> return the resulting configuration and scalar factor.
    procedure(I_ExciteConfiguration), deferred :: ExciteConfiguration

  end type

  abstract interface
    !> Runtime setup hook executed after all elements have been fabricated and
    !> registered. Implementations may allocate resources and compute caches.
    subroutine I_Setup(this)
      import :: T_ConfigList_E
      !> The configuration element to set up.
      class(T_ConfigList_E), intent(inout) :: this
    end subroutine

    !> Fabrication hook to create and initialize element data from input. This
    !> is the place to read parameters, size arrays, and build initial mapping
    !> tables.
    subroutine I_Fabricate(this)
      import :: T_ConfigList_E
      !> The configuration element to construct.
      class(T_ConfigList_E), intent(inout) :: this
    end subroutine

    !> Apply creation and annihilation operators to transform a configuration.
    !>
    !> Returns the new configuration index and associated scalar factor. For
    !> fermions the factor includes the parity/sign due to anti-commutation; for
    !> bosons it includes any ladder prefactors but no sign. A result of
    !> `iCNew = 0` indicates the excitation is not allowed (e.g., Pauli
    !> violation or leaving the active space).
    subroutine I_ExciteConfiguration(this, iCNew, factor, creates, destroys, iC)
      import :: I32, I64, R64, T_ConfigList_E
      !> Configuration element instance (statistics-specific implementation).
      class(T_ConfigList_E), intent(in) :: this
      !> Resulting configuration index (0 if excitation not possible).
      integer(I32), intent(out)        :: iCNew
      !> Scalar factor to apply to the state coefficient.
      real(R64), intent(out)           :: factor
      !> Orbital indices for creation operators (a^\dagger).
      integer(I32), intent(in), contiguous :: creates(:)
      !> Orbital indices for annihilation operators (a).
      integer(I32), intent(in), contiguous :: destroys(:)
      !> Input configuration index the excitation is applied to.
      integer(I32), intent(in)         :: iC
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Simple container enabling polymorphic arrays of configuration elements.
  type :: T_ConfigList_Container
    !> Polymorphic configuration element (allocated to a concrete implementation).
    class(T_ConfigList_E), allocatable :: e
  end type

  !> Global list of configuration providers composing the many-body basis.
  type(T_ConfigList_Container), allocatable :: configList(:)

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> Optional hook executed during the overall setup phase of the diagonalization
  !> driver. Defaults to a no-op and can be redirected by higher-level modules
  !> to perform additional global initialization once all elements are ready.
  procedure(I_ConfigList_Setup), pointer :: ConfigList_Setup => NoOpProcedures_Setup
  abstract interface
    !> Global setup stage for configuration-related infrastructure.
    subroutine I_ConfigList_Setup
    end subroutine
  end interface

end module
