! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file M_DiagonalizerList.f90
!> @brief Abstract interface for eigenvalue problem solvers (diagonalizers).
!>
!> @details
!> This module defines the public API and abstract type for eigensolvers used
!> throughout CodyFortran. It implements a **strategy pattern** where the
!> concrete eigensolver algorithm (LAPACK, ARPACK, etc.) is selected at runtime
!> via JSON configuration and bound through polymorphism.
!>
!> **Key Design Principles:**
!> - Matrix-free operation: Operators are accessed only via callback `ApplyMatOnVec`
!> - Time-dependent support: All interfaces accept a time parameter for propagators
!> - Partial spectrum: Solvers compute a requested subset of eigenpairs
!>
!> **Workflow:**
!> 1. Client creates `T_DiagonalizerList_FabricateInput` with callback + dimension
!> 2. `DiagonalizerList_Fabricate(input)` allocates backends per JSON config
!> 3. `DiagonalizerList_Setup()` initializes all backends
!> 4. Optionally `DiagonalizerList(i)%e%SetStartVectors(vecs)`: start vectors
!>    for the next call (used by backends that can start from given vectors)
!> 5. `DiagonalizerList(i)%e%Diagonalize(time, evecsQ)` computes eigenpairs
!> 6. Results in `%evals(:)`, `%evecs(:,:)`, `%nFound`
!>
!> **Available Backends:**
!> - LAPACK: Dense solver for small-medium matrices (explicit assembly)
!> - ARPACK: Iterative Arnoldi/Lanczos for large sparse/matrix-free problems
!> - LOBPCG: Preconditioned block eigensolver for the lowest eigenpairs of
!>   operators that are self-adjoint in a (possibly weighted) metric
!>
!> **Optional inputs** (ignored by LAPACK and ARPACK):
!> - `ApplyPrecOnVec` (callback): preconditioner, approximately (A − shift)⁻¹
!> - `metricWeights` (element component, set after fabrication): the diagonal
!>   metric in which A is self-adjoint (e.g. a grid's quadrature weights);
!>   Euclidean when unallocated, the convention of IntegratorList
!>
!> @see M_DiagonalizerList_Lapack, M_DiagonalizerList_Arpack,
!>      M_DiagonalizerList_Lobpcg for implementations
!> @see A_DiagonalizerList.md in this directory for detailed onboarding documentation
module M_DiagonalizerList
  use M_Utils_Types
  use M_Utils_NoOpProcedures, only: NoOpProcedures_Setup

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  !-----------------------------------------------------------------------------
  !> @brief Input specification for creating a diagonalizer instance.
  !>
  !> @details
  !> Each element of the input array to `DiagonalizerList_Fabricate` describes
  !> one eigensolver to create. The client provides:
  !> - A callback that applies the linear operator (Hamiltonian) to a vector
  !> - The dimension of the vector space (Hilbert space size)
  !>
  !> **Example usage:**
  !> @code{.f90}
  !> type(T_DiagonalizerList_FabricateInput) :: inputs(1)
  !> inputs(1)%ApplyMatOnVec => MyHamiltonianCallback
  !> inputs(1)%dim = hilbert_space_size
  !> ! optional, used by preconditioned solvers (lobpcg) only:
  !> inputs(1)%ApplyPrecOnVec => MyPreconditioner
  !> call DiagonalizerList_Fabricate(inputs)
  !> @endcode
  !-----------------------------------------------------------------------------
  type :: T_DiagonalizerList_FabricateInput
    !> Procedure that applies the linear operator to a vector: y = A(t) · x.
    !> For stationary (time-independent) problems, the time argument is ignored.
    !> Must conform to interface `I_ApplyMatOnVec`.
    procedure(I_ApplyMatOnVec), pointer, nopass :: ApplyMatOnVec => null()
    !> Dimension of the operator matrix (size of the Hilbert space).
    integer(I32) :: dim = 0
    !> Optional preconditioner: y ≈ (A(t) − shift)⁻¹ · x (interface
    !> `I_ApplyPrecOnVec`). Unassociated means no preconditioning (identity).
    !> Only used by preconditioned backends (lobpcg).
    procedure(I_ApplyPrecOnVec), pointer, nopass :: ApplyPrecOnVec => null()
  end type

  interface
    !---------------------------------------------------------------------------
    !> @brief Constructs and initializes the global list of diagonalizer backends.
    !>
    !> @details
    !> Factory subroutine that reads the JSON configuration at path
    !> `"diagonalizerList"` and creates one eigensolver instance per JSON child.
    !> The backend type (LAPACK vs ARPACK) is determined by the child name.
    !>
    !> **JSON Structure:**
    !> @code{.json}
    !> {
    !>   "diagonalizerList": {
    !>     "arpack_ground": { "nEvals": 5, "which": "SR", ... },
    !>     "lapack_full":   { "nEvals": -1, ... }
    !>   }
    !> }
    !> @endcode
    !>
    !> **Name convention:** Child name containing "lapack" → LAPACK backend,
    !> "arpack" → ARPACK backend, "lobpcg" → LOBPCG backend.
    !>
    !> @param[in] input  Array of fabrication inputs (one per JSON child).
    !>                   `size(input)` must equal number of JSON children.
    !>
    !> @pre JSON file loaded via `Json_LoadJsonFile`
    !> @post `DiagonalizerList(:)` allocated and backends initialized
    !---------------------------------------------------------------------------
    module subroutine DiagonalizerList_Fabricate(input)
      type(T_DiagonalizerList_FabricateInput), intent(in) :: input(:)
    end subroutine
  end interface

  !=============================================================================
  ! type definition
  !=============================================================================

  !-----------------------------------------------------------------------------
  !> @brief Abstract base type for all eigensolver backends.
  !>
  !> @details
  !> Defines the common interface that all diagonalizer implementations must
  !> provide. The design supports:
  !> - **Matrix-free operation:** Only requires `ApplyMatOnVec` callback
  !> - **Partial spectrum computation:** Request `nEvals` eigenpairs
  !> - **Time-dependent operators:** `time` parameter in `Diagonalize`
  !>
  !> **Inheritance hierarchy:**
  !> @code
  !> T_DiagonalizerList_E (abstract)
  !>   ├── T_DiagonalizerList_E_Lapack   (dense, explicit matrix)
  !>   ├── T_DiagonalizerList_E_Arpack   (iterative, Arnoldi/Lanczos)
  !>   └── T_DiagonalizerList_E_Lobpcg   (iterative, preconditioned block)
  !> @endcode
  !>
  !> **Output convention:**
  !> After `Diagonalize`, eigenvalues are in `evals(1:nFound)` sorted in
  !> ascending order (for real parts). Eigenvectors (if requested) are in
  !> `evecs(:, 1:nFound)` with `evecs(:,j)` corresponding to `evals(j)`.
  !-----------------------------------------------------------------------------
  type, abstract :: T_DiagonalizerList_E
    !> Callback that performs matrix–vector product: y = A(t) · x.
    !> Used by iterative solvers directly or to assemble explicit matrices.
    procedure(I_ApplyMatOnVec), pointer, nopass :: ApplyMatOnVec => null()

    !> Optional preconditioner y ≈ (A(t) − shift)⁻¹ · x (copied from the
    !> fabrication input; unassociated = identity). Ignored by LAPACK/ARPACK.
    procedure(I_ApplyPrecOnVec), pointer, nopass :: ApplyPrecOnVec => null()

    !> Optional metric weights (size dim) defining the inner product in which
    !> the operator is self-adjoint, for metric-aware backends (lobpcg): when
    !> allocated, inner products are sum(conjg(a) * metricWeights * b); when
    !> unallocated, the Euclidean product is used. The same contract as
    !> IntegratorList's `metricWeights`: required on grids whose quadrature
    !> weights are not absorbed into the vectors (e.g. FEDVR, finite
    !> elements: A = W⁻¹H is self-adjoint only in the weighted metric).
    !> Assign directly on the registry element after fabrication and before
    !> the first `Diagonalize` (grid weights typically exist only after
    !> `Grid_Setup`). Ignored by LAPACK/ARPACK.
    real(R64), allocatable :: metricWeights(:)

    !> Computed eigenvalues (real, sorted ascending by default).
    !> Allocated to size `nEvals`, but only `nFound` entries are valid.
    real(R64), allocatable :: evals(:)

    !> Computed eigenvectors (column-major: evecs(:,j) is j-th eigenvector).
    !> Shape (dim, nEvals), with `nFound` valid columns after diagonalization.
    !> Only allocated/filled if `evecsQ=.true.` in `Diagonalize`. An output:
    !> clients pass start vectors through `SetStartVectors`.
    complex(R64), allocatable :: evecs(:, :)

    !> Number of eigenpairs requested. Set to `-1` for all eigenvalues.
    integer(I32) :: nEvals = 0

    !> Number of converged eigenpairs actually found (≤ nEvals).
    integer(I32) :: nFound = 0

    !> Whether the last `Diagonalize` call converged all `nEvals` wanted
    !> eigenpairs. Backends that can return unconverged approximations set it
    !> on every call (lobpcg: `.false.` at maxIter or on stagnation, with
    !> nFound < nEvals); LAPACK is direct and ARPACK stops instead, so they
    !> keep the default `.true.`. Clients that may continue with
    !> approximations (e.g. the SCF ground solvers) check it after every call.
    logical :: convergedQ = .true.

    !> Dimension of the operator matrix (Hilbert space size).
    integer(I32) :: dim = 0

    !> Verbosity level: 0 = silent, 1 = progress, 2+ = debug.
    integer(I32) :: printLevel = 0

    !> JSON configuration path for this instance (e.g., "diagonalizerList.arpack1").
    character(len=:), allocatable :: path
  contains
    !> Read configuration from JSON (at `this%path`) and set parameters.
    procedure(I_Fabricate), deferred :: Fabricate
    !> Allocate working arrays and prepare for diagonalization.
    procedure(I_Setup), deferred :: Setup
    !> Compute eigenvalues and optionally eigenvectors.
    procedure(I_Diagonalize), deferred :: Diagonalize
    !> Start vectors for the next `Diagonalize` (default: ignored).
    procedure :: SetStartVectors => DiagonalizerList_E_SetStartVectors
  end type

  interface
    !---------------------------------------------------------------------------
    !> @brief Start vectors for the next `Diagonalize` (binding
    !>        `SetStartVectors`).
    !>
    !> @details
    !> Contract: the columns of `vecs` (dim rows, any number of columns) are
    !> offered as start vectors for the next `Diagonalize` call, e.g. the
    !> current orbitals for the first solve of an SCF, before the solver holds
    !> eigenvectors of its own. A backend may ignore them: this default is a
    !> no-op (LAPACK is direct; ARPACK warm-starts only from its own converged
    !> eigenvectors). LOBPCG overrides it, but starts from the vectors only
    !> when its `warmStartQ` is true (the default), the switch that also
    !> governs the reuse of its previous eigenvectors; with `warmStartQ`
    !> false they are stored and ignored.
    !>
    !> @param[inout] this  Diagonalizer instance
    !> @param[in]    vecs  Start vectors, one per column (dim rows)
    !---------------------------------------------------------------------------
    module subroutine DiagonalizerList_E_SetStartVectors(this, vecs)
      class(T_DiagonalizerList_E), intent(inout) :: this
      complex(R64), intent(in), contiguous :: vecs(:, :)
    end subroutine
  end interface

  !-----------------------------------------------------------------------------
  ! Abstract Interfaces
  !-----------------------------------------------------------------------------

  abstract interface
    !> @brief Initialize diagonalizer from JSON configuration.
    !> @details Called by `DiagonalizerList_Fabricate` after allocation.
    !>          Implementation reads parameters from `this%path` and sets
    !>          `nEvals`, `printLevel`, and backend-specific options.
    !> @param[inout] this  Diagonalizer instance to initialize
    subroutine I_Fabricate(this)
      import :: T_DiagonalizerList_E
      class(T_DiagonalizerList_E), intent(inout) :: this
    end subroutine
  end interface

  abstract interface
    !> @brief Prepare backend for diagonalization.
    !> @details Called after fabrication. May allocate working arrays,
    !>          precompute factorizations, or initialize iterative state.
    !> @param[inout] this  Diagonalizer instance to set up
    subroutine I_Setup(this)
      import :: T_DiagonalizerList_E
      class(T_DiagonalizerList_E), intent(inout) :: this
    end subroutine
  end interface

  abstract interface
    !> @brief Compute eigenvalues and optionally eigenvectors.
    !>
    !> @details
    !> Core computational routine. After successful completion:
    !> - `this%nFound` contains number of converged eigenpairs
    !> - `this%evals(1:nFound)` contains eigenvalues (sorted ascending)
    !> - `this%evecs(:,1:nFound)` contains eigenvectors (if `evecsQ=.true.`)
    !> - `this%convergedQ` tells whether all `nEvals` wanted pairs converged
    !>   (lobpcg keeps all nEvals best approximations in `evals`/`evecs` and
    !>   counts the converged ones in `nFound`, see M_DiagonalizerList_Lobpcg)
    !>
    !> @param[inout] this    Diagonalizer instance
    !> @param[in]    time    Time parameter for time-dependent operators.
    !>                       Ignored for stationary Hamiltonians.
    !> @param[in]    evecsQ  If `.true.`, compute and store eigenvectors
    subroutine I_Diagonalize(this, time, evecsQ)
      import :: R64, T_DiagonalizerList_E
      class(T_DiagonalizerList_E), intent(inout) :: this
      real(R64), intent(in) :: time
      logical, intent(in) :: evecsQ
    end subroutine
  end interface

  abstract interface
    !> @brief Matrix–vector multiplication callback: y = A(t) · x.
    !>
    !> @details
    !> The fundamental operation for matrix-free eigensolvers. Implementations
    !> apply the linear operator (typically a Hamiltonian) to an input vector.
    !>
    !> **Contract:**
    !> - `dState` is overwritten with result (no aliasing with `state`)
    !> - `state` is read-only input vector
    !> - For time-independent operators, `time` may be ignored
    !>
    !> @param[out] dState  Result vector y = A(t) · x
    !> @param[in]  state   Input vector x
    !> @param[in]  time    Time at which operator is evaluated
    subroutine I_ApplyMatOnVec(dState, state, time)
      import :: R64
      complex(R64), intent(out), contiguous, target :: dState(:)
      complex(R64), intent(in), contiguous, target  :: state(:)
      real(R64), intent(in)                         :: time
    end subroutine
  end interface

  abstract interface
    !> @brief Preconditioner callback: y ≈ (A(t) − shift)⁻¹ · x.
    !>
    !> @details
    !> Applied by preconditioned eigensolvers (LOBPCG) to the residual of a
    !> Ritz pair, with `shift` the current Ritz value. Any approximation that
    !> is Hermitian positive definite in the solver metric for shifts below
    !> the wanted part of the spectrum works (e.g. (T + σ)⁻¹ for a Fock
    !> operator with kinetic part T); only the direction of the result
    !> matters. Same aliasing contract as `I_ApplyMatOnVec`.
    !>
    !> @param[out] dState  Preconditioned vector y
    !> @param[in]  state   Input vector x (a residual)
    !> @param[in]  shift   Current Ritz value the preconditioner may adapt to
    !> @param[in]  time    Time at which the operator is evaluated
    subroutine I_ApplyPrecOnVec(dState, state, shift, time)
      import :: R64
      complex(R64), intent(out), contiguous, target :: dState(:)
      complex(R64), intent(in), contiguous, target  :: state(:)
      real(R64), intent(in)                         :: shift
      real(R64), intent(in)                         :: time
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !-----------------------------------------------------------------------------
  !> @brief Container for polymorphic diagonalizer storage.
  !> @details Fortran requires this wrapper to create arrays of polymorphic
  !>          objects with different dynamic types.
  !-----------------------------------------------------------------------------
  type :: T_DiagonalizerList_Container
    !> Polymorphic eigensolver instance (LAPACK, ARPACK, or other backend).
    class(T_DiagonalizerList_E), allocatable :: e
  end type

  !-----------------------------------------------------------------------------
  !> @brief Global array of diagonalizer instances.
  !>
  !> @details
  !> Populated by `DiagonalizerList_Fabricate`. Access pattern:
  !> @code{.f90}
  !> call DiagonalizerList(1)%e%Diagonalize(0.0_R64, .true.)
  !> evals => DiagonalizerList(1)%e%evals
  !> @endcode
  !>
  !> Array index corresponds to order in JSON `"diagonalizerList"` block.
  !-----------------------------------------------------------------------------
  type(T_DiagonalizerList_Container), allocatable :: DiagonalizerList(:)

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !-----------------------------------------------------------------------------
  !> @brief Global setup procedure pointer.
  !>
  !> @details
  !> Calls `Setup` on all diagonalizer instances in `DiagonalizerList(:)`.
  !> Bound during fabrication; defaults to no-op if fabrication not called.
  !>
  !> **Typical call sequence:**
  !> @code{.f90}
  !> call DiagonalizerList_Fabricate(inputs)  ! Binds this pointer
  !> call DiagonalizerList_Setup              ! Calls all backends' Setup
  !> @endcode
  !-----------------------------------------------------------------------------
  procedure(I_DiagonalizerList_Setup), pointer :: DiagonalizerList_Setup => NoOpProcedures_Setup

  abstract interface
    !> @brief Interface for global setup procedure.
    subroutine I_DiagonalizerList_Setup
    end subroutine
  end interface

end module

