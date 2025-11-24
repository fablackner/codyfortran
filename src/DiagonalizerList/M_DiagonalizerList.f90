! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Module M_DiagonalizerList defines the abstract interface and shared data for
!> eigenvalue problem solvers ("diagonalizers").
!> Concrete backends (e.g., ARPACK, LAPACK) implement the deferred procedures to
!> compute selected eigenvalues and optionally eigenvectors of a linear operator
!> that is provided through a matrix–vector callback.
module M_DiagonalizerList
  use M_Utils_Types
  use M_Utils_NoOpProcedures, only: NoOpProcedures_Setup

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  type :: T_DiagonalizerList_FabricateInput
    !> Procedure that applies the linear operator to a vector, y = A(t) x.
    !> For stationary problems the time argument can be ignored by the callback.
    procedure(I_ApplyMatOnVec), pointer, nopass :: ApplyMatOnVec
    !> Dimension of the operator matrix being diagonalized.
    integer(I32) :: dim
  end type

  interface
    !> Constructs and initializes the global list of diagonalizer backends.
    !> One input element corresponds to one diagonalizer instance. This is
    !> typically called during program initialization to allocate backends,
    !> record their dimensions, and store the operator callbacks.
    module subroutine DiagonalizerList_Fabricate(input)
      !> Inputs describing each backend to be created (operator callback and
      !> matrix dimension).
      type(T_DiagonalizerList_FabricateInput), intent(in) :: input(:)
    end subroutine
  end interface

  !=============================================================================
  ! type definition
  !=============================================================================

  !> Abstract type defining the interface for eigensolvers used to
  !> diagonalize a linear operator. Concrete implementations must provide
  !> methods for initialization, setup, and the actual diagonalization.
  type, abstract :: T_DiagonalizerList_E
    !> Procedure that performs the matrix–vector product y = A(t) x used by
    !> iterative solvers or for assembling explicit matrices.
    procedure(I_ApplyMatOnVec), pointer, nopass :: ApplyMatOnVec
    !> Eigenvalues computed by "Diagonalize" (length nFound ≤ nEvals).
    !> evals(i) is the i-th eigenvalue.
    real(R64), allocatable :: evals(:)
    !> Eigenvectors computed by "Diagonalize" when requested.
    !> evecs(i,j) is the i-th component of the j-th eigenvector.
    complex(R64), allocatable :: evecs(:, :)
    !> Number of eigenvalues/eigenvectors requested for computation.
    integer(I32) :: nEvals
    !> Number of converged eigenpairs found after diagonalization.
    integer(I32) :: nFound
    !> Dimension of the operator matrix being diagonalized.
    integer(I32) :: dim
    !> Level of verbosity for printing information during setup/solve.
    integer(I32)  :: printLevel
    !> Path to this object in configuration file.
    character(len=:), allocatable :: path
  contains
    !> Initializes the diagonalizer with configuration parameters
    !> (e.g., reads from this%path) and sets core fields such as dim,
    !> nEvals, printLevel, and ApplyMatOnVec.
    procedure(I_Fabricate), deferred :: Fabricate
    !> Performs backend-specific setup after initialization (e.g., allocate
    !> work arrays, precompute shifts/preconditioners).
    procedure(I_Setup), deferred :: Setup
    !> Computes eigenvalues and optionally eigenvectors according to the
    !> requested criteria, filling evals, evecs (if requested), and nFound.
    procedure(I_Diagonalize), deferred :: Diagonalize
  end type

  abstract interface
    !> Initializes a diagonalizer instance from configuration (referenced by
    !> this%path if applicable). Implementations should bind ApplyMatOnVec and
    !> set dim, nEvals, printLevel, and any backend parameters.
    subroutine I_Fabricate(this)
      import :: T_DiagonalizerList_E
      !> The diagonalizer instance to initialize
      class(T_DiagonalizerList_E), intent(inout) :: this
    end subroutine
  end interface

  abstract interface
    !> Performs setup operations for the diagonalizer after initialization.
    !> May allocate resources, working memory, or precompute factors.
    subroutine I_Setup(this)
      import :: T_DiagonalizerList_E
      !> The diagonalizer instance to set up
      class(T_DiagonalizerList_E), intent(inout) :: this
    end subroutine
  end interface

  abstract interface
    subroutine I_Diagonalize(this, time, evecsQ)
      import :: R64, T_DiagonalizerList_E
      !> The diagonalizer instance
      class(T_DiagonalizerList_E), intent(inout) :: this
      !> Time at which to evaluate the operator if it is time dependent.
      !> For time-independent problems, implementations may ignore this value.
      real(R64), intent(in) :: time
      !> Flag indicating whether to compute eigenvectors in addition to eigenvalues.
      logical, intent(in) :: evecsQ
    end subroutine
  end interface

  abstract interface
    !> Applies the linear operator to a vector.
    !> Given a state vector x, compute y = A(t) x, where A may be
    !> time dependent. Implementations may ignore the time argument for
    !> stationary operators.
    subroutine I_ApplyMatOnVec(dState, state, time)
      import :: R64
      !> Output vector y = A(t) x (overwritten on entry).
      complex(R64), intent(out), contiguous, target :: dState(:)
      !> Input vector x.
      complex(R64), intent(in), contiguous, target  :: state(:)
      !> Time at which the operator is evaluated.
      real(R64), intent(in)             :: time
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Container type for polymorphic diagonalizer objects.
  !> Enables arrays of different concrete backends.
  type :: T_DiagonalizerList_Container
    !> Polymorphic instance of a T_DiagonalizerList_E implementation.
    class(T_DiagonalizerList_E), allocatable :: e
  end type

  !> Global registry of diagonalizer instances used by the program.
  type(T_DiagonalizerList_Container), allocatable :: DiagonalizerList(:)

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> Pointer to a procedure that performs global setup for the diagonalization
  !> subsystem. Defaults to a no-op and can be redirected by backends.
  procedure(I_DiagonalizerList_Setup), pointer :: DiagonalizerList_Setup => NoOpProcedures_Setup
  abstract interface
    !> Performs global initialization for the diagonalization subsystem.
    subroutine I_DiagonalizerList_Setup
    end subroutine
  end interface

end module

