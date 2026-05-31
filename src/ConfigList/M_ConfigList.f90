! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file M_ConfigList.f90
!> @brief Core Fock-basis enumeration and excitation machinery.
!>
!> @details
!> ## Purpose
!>
!> This module provides the abstract foundation for enumerating many-body
!> configurations (Fock states) and applying second-quantized excitation
!> operators. It is central to configuration-interaction (CI), MCTDH, and
!> TD-2RDM methods where Hamiltonian matrix elements are computed via
!> `<I| a†_i a_j |J>` (singles) or `<I| a†_i a†_j a_l a_k |J>` (doubles).
!>
!> ## Responsibilities
!>
!> 1. **Abstract base class (`T_ConfigList_E`)** defines the contract all
!>    configuration implementations must satisfy: fabricate from JSON, runtime
!>    setup, and excitation application.
!>
!> 2. **Excitation connectivity (`T_ConfigList_Excitation`)** precomputes and
!>    stores sparse graphs linking configurations via single and double
!>    excitations for O(1) Hamiltonian construction.
!>
!> 3. **Global container (`configList`)** collects polymorphic configuration
!>    elements, enabling multi-species systems (mixed bosons/fermions).
!>
!> ## Encoding Schemes
!>
!> Concrete implementations encode configurations as integers:
!>
!> - **Fermions:** Bit-patterns where bit `i` is set if orbital `i` is occupied.
!>   Excitations include ±1 phase from anti-commutation.
!>
!> - **Bosons:** Base-(N+1) encoding where digit `i` gives occupation of orbital
!>   `i`. Excitations produce √n ladder factors.
!>
!> ## Usage Flow
!>
!> ```
!> ConfigList_Fabricate   → parses JSON, allocates concrete elements
!>     ↓
!> ConfigList_Setup       → computes codeFromConfig, builds singles/doubles
!>     ↓
!> ExciteConfiguration    → applies a†/a operators at runtime
!> ```
!>
!> ## Key Data Structures
!>
!> | Name                      | Type                         | Description                        |
!> |---------------------------|------------------------------|------------------------------------|
!> | `configList(:)`           | `T_ConfigList_Container[]`   | Global array of config elements    |
!> | `codeFromConfig(:)`       | `integer(I64)`               | Map iC → bit/base encoding         |
!> | `singles`, `doubles`      | `T_ConfigList_Excitation`    | Precomputed connectivity graphs    |
!>
!> @see M_ConfigList_AllActive for all-active truncation
!> @see M_Coeffs for CI coefficient storage consuming these configurations
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

  !---------------------------------------------------------------------------
  !> @brief Sparse connectivity graph for excitations of a given rank.
  !>
  !> @details
  !> Stores the precomputed mapping from each configuration iC to all
  !> configurations reachable via single (rank=1) or double (rank=2)
  !> excitations. This enables O(1) lookup during Hamiltonian construction
  !> and RDM evaluation.
  !>
  !> **Memory layout (column-major for cache efficiency):**
  !> ```
  !>   nConnected(iC)       → number of non-zero connections from config iC
  !>   excitedC(k, iC)      → target config index for k-th connection

  !>   factor(k, iC)        → scalar factor (includes fermionic sign or bosonic √n)
  !>   orbCode(k, iC)       → packed orbital indices defining the excitation
  !> ```
  !>
  !> **orbCode encoding (for nOrbs total orbitals):**
  !> - Singles: `orbCode = i * nOrbs + j` encodes a†_i a_j
  !> - Doubles: `orbCode = ((i1 * nOrbs + j1) * nOrbs + i2) * nOrbs + j2`
  !>   encodes a†_i1 a†_i2 a_j2 a_j1 (note operator ordering convention)
  !>
  !> The complex type for `factor` supports future complex Hamiltonians; real
  !> factors are stored with zero imaginary part.
  !---------------------------------------------------------------------------
  type :: T_ConfigList_Excitation
    !> Number of configurations connected to each source configuration.
    !> Shape: (nConfigurations)
    integer(I32), allocatable :: nConnected(:)

    !> Target configuration indices for each connection edge.
    !> Shape: (maxConnections, nConfigurations)
    integer(I32), allocatable :: excitedC(:, :)

    !> Scalar factors including statistics-dependent phase/ladder prefactors.
    !> Shape: (maxConnections, nConfigurations)
    complex(R64), allocatable :: factor(:, :)

    !> Packed orbital indices defining the excitation pattern per edge.
    !> Shape: (maxConnections, nConfigurations)
    integer(I32), allocatable :: orbCode(:, :)
  end type

  !---------------------------------------------------------------------------
  !> @brief Abstract base class for all Fock-basis configuration implementations.
  !>
  !> @details
  !> This is the polymorphic root of the configuration hierarchy. Concrete
  !> implementations (fermionic, bosonic, spin-adapted, ...) extend this type
  !> and provide statistics-specific excitation logic.
  !>
  !> ## Lifecycle
  !>
  !> 1. **Allocation:** Factory routine allocates concrete type based on JSON.
  !> 2. **Fabricate():** Reads parameters, computes nConfigurations.
  !> 3. **Setup():** Builds codeFromConfig mapping, precomputes singles/doubles.
  !> 4. **ExciteConfiguration():** Called at runtime to apply operators.
  !>
  !> ## Thread Safety
  !>
  !> After Setup completes, all data is read-only. ExciteConfiguration is
  !> thread-safe for concurrent calls with different iC values.
  !>
  !> ## Extension Points
  !>
  !> To add a new statistics type:
  !> 1. Create `T_ConfigList_E_YourType` extending this base
  !> 2. Implement `Setup`, `Fabricate`, `ExciteConfiguration`
  !> 3. Register in `ConfigList_Fabricate` branching logic
  !---------------------------------------------------------------------------
  type, abstract :: T_ConfigList_E

    !> Index of the particle species this element handles (1-based).
    !> Multi-species systems have one element per bodyTarget.
    integer(I32) :: bodyTarget

    !> Total number of configurations (Fock states) in this basis.
    !> For N particles in M orbitals: C(M,N) fermions, C(M+N-1,N) bosons.
    integer(I32) :: nConfigurations

    !> Maps linear config index iC ∈ [1, nConfigurations] to internal encoding.
    !>
    !> - Fermions: bit i set ⟺ orbital i occupied (integer as bit-field)
    !> - Bosons: base-(N+1) encoding where digit i = occupation of orbital i
    !>
    !> This enables O(1) occupation lookup and excitation validation.
    integer(I64), allocatable :: codeFromConfig(:)

    !> Precomputed single-excitation connectivity graph.
    !> Populated during Setup by iterating all (i,j) pairs.
    type(T_ConfigList_Excitation) :: singles

    !> Precomputed double-excitation connectivity graph.
    !> Populated during Setup by iterating all (i1,i2,j1,j2) tuples.
    type(T_ConfigList_Excitation) :: doubles

    !> JSON path from which this element was constructed (for diagnostics).
    !> Example: "configList.allActive1.fermionic"
    character(len=:), allocatable :: path

  contains

    !> Perform runtime setup after all elements are fabricated.
    !> Builds codeFromConfig and populates singles/doubles connectivity.
    procedure(I_Setup), deferred :: Setup

    !> Construct element from JSON input. Called before Setup.
    !> Reads bodyTarget, nExcitations; computes nConfigurations.
    procedure(I_Fabricate), deferred :: Fabricate

    !> Apply second-quantized operators to transform configuration.
    !>
    !> Given |iC⟩, computes |iCNew⟩ = a†_{creates} a_{destroys} |iC⟩
    !> and returns the scalar factor (sign for fermions, √n for bosons).
    !>
    !> Returns iCNew=0 if excitation is forbidden (Pauli, empty orbital, ...).
    procedure(I_ExciteConfiguration), deferred :: ExciteConfiguration

  end type

  abstract interface
    !---------------------------------------------------------------------------
    !> @brief Runtime setup hook for configuration elements.
    !>
    !> @details
    !> Called after all elements have been fabricated and the global context
    !> (orbital counts, particle numbers) is available. Typical tasks:
    !> - Allocate and populate `codeFromConfig(:)`
    !> - Build singles/doubles connectivity graphs
    !> - Verify configuration count matches expectations
    !---------------------------------------------------------------------------
    subroutine I_Setup(this)
      import :: T_ConfigList_E
      class(T_ConfigList_E), intent(inout) :: this
    end subroutine

    !---------------------------------------------------------------------------
    !> @brief Fabrication hook to initialize element from JSON input.
    !>
    !> @details
    !> Called during the build phase before `Setup`. This is where to:
    !> - Read `bodyTarget` and `nExcitations` from JSON
    !> - Compute `nConfigurations` using combinatorial formulas
    !> - Validate that bodyStatistics matches element type
    !---------------------------------------------------------------------------
    subroutine I_Fabricate(this)
      import :: T_ConfigList_E
      class(T_ConfigList_E), intent(inout) :: this
    end subroutine

    !---------------------------------------------------------------------------
    !> @brief Apply second-quantized excitation operators to a configuration.
    !>
    !> @details
    !> Computes the action of a†_{creates} a_{destroys} on configuration |iC⟩.
    !>
    !> ## Operator Ordering Convention
    !>
    !> The operators are applied right-to-left (physics convention):
    !> ```
    !>   |iCNew⟩ = a†_{c_n} ... a†_{c_1} a_{d_1} ... a_{d_m} |iC⟩
    !> ```
    !> where `creates = [c_1, ..., c_n]` and `destroys = [d_1, ..., d_m]`.
    !>
    !> ## Return Values
    !>
    !> - `iCNew > 0`: Valid excitation; `factor` contains the scalar prefactor
    !>   - Fermions: factor = ±1 (parity from anti-commutation)
    !>   - Bosons: factor = √(product of occupation factors)
    !>
    !> - `iCNew = 0`: Excitation forbidden (Pauli exclusion, annihilating empty
    !>   orbital, or leaving the active space)
    !>
    !> ## Thread Safety
    !>
    !> This routine is thread-safe for concurrent calls with different iC.
    !---------------------------------------------------------------------------
    subroutine I_ExciteConfiguration(this, iCNew, factor, creates, destroys, iC)
      import :: I32, I64, R64, T_ConfigList_E
      class(T_ConfigList_E), intent(in)      :: this
      integer(I32), intent(out)              :: iCNew
      real(R64), intent(out)                 :: factor
      integer(I32), intent(in), contiguous   :: creates(:)
      integer(I32), intent(in), contiguous   :: destroys(:)
      integer(I32), intent(in)               :: iC
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !---------------------------------------------------------------------------
  !> @brief Polymorphic container for configuration elements.
  !>
  !> @details
  !> Fortran requires a wrapper type to store polymorphic objects in arrays.
  !> Each container holds one concrete configuration element (fermionic,
  !> bosonic, etc.). The global `configList(:)` array collects all species.
  !---------------------------------------------------------------------------
  type :: T_ConfigList_Container
    class(T_ConfigList_E), allocatable :: e
  end type

  !---------------------------------------------------------------------------
  !> @brief Global registry of configuration providers.
  !>
  !> @details
  !> Populated by `ConfigList_Fabricate` based on JSON input. Each element
  !> handles one particle species. Multi-species CI uses the tensor product:
  !> ```
  !>   |I⟩ = |i_1⟩ ⊗ |i_2⟩ ⊗ ... (one config per species)
  !> ```
  !> Total configuration count = product of individual nConfigurations.
  !---------------------------------------------------------------------------
  type(T_ConfigList_Container), allocatable :: configList(:)

  !=============================================================================
  ! module procedure pointers
  !=============================================================================

  !---------------------------------------------------------------------------
  !> @brief Global setup hook for configuration infrastructure.
  !>
  !> @details
  !> Called during the overall setup phase after fabrication. Defaults to a
  !> no-op; the submodule redirects this to iterate over all elements and
  !> call their individual Setup methods.
  !---------------------------------------------------------------------------
  procedure(I_ConfigList_Setup), pointer :: ConfigList_Setup => NoOpProcedures_Setup

  abstract interface
    subroutine I_ConfigList_Setup
    end subroutine
  end interface

end module
