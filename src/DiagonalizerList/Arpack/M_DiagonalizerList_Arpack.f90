! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Module M_DiagonalizerList_Arpack provides an ARPACK-based iterative
!> eigensolver backend. It implements the abstract
!> T_DiagonalizerList_E interface using Arnoldi/Lanczos methods to
!> compute a selected subset of eigenvalues and (optionally) eigenvectors
!> of a linear operator accessed through a matrix–vector callback.
module M_DiagonalizerList_Arpack
  use M_Utils_Types
  use M_DiagonalizerList, only: T_DiagonalizerList_E

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Allocates a new ARPACK diagonalizer instance and associates it with a
    !> configuration path. Intended to be called by the list fabrication logic.
    module subroutine DiagonalizerList_Arpack_Allocate(e, path)
      class(T_DiagonalizerList_E), allocatable, intent(out) :: e
      character(len=*), intent(in) :: path
    end subroutine
  end interface

  !=============================================================================
  ! module types
  !=============================================================================

  !> ARPACK diagonalizer.
  !> Wraps an implicitly restarted Arnoldi/Lanczos iteration to compute the
  !> requested number of eigenpairs of large sparse or matrix-free operators.
  type, extends(T_DiagonalizerList_E) :: T_DiagonalizerList_E_Arpack
    !> Which Ritz values to target (e.g., 'LM','SM','LR','SR','LI','SI').
    character(2)  :: which
    !> Matrix type flag for ARPACK: 'I' for standard problem, 'G' for generalized.
    character(1)  :: bmat
    !> Dimension of the Krylov subspace (number of Arnoldi vectors).
    integer(I32)  :: nKry
    !> If true, perform explicit convergence checks and optional retries.
    logical       :: checkConvergenceQ
  contains
    !> Backend-specific setup (allocate work arrays, set ARPACK parameters).
    procedure :: Setup
    !> Read/derive ARPACK parameters from configuration and bind callbacks.
    procedure :: Fabricate
    !> Run the ARPACK iteration to compute eigenvalues (and optionally eigenvectors).
    procedure :: Diagonalize
  end type

  interface
    !> Initializes the ARPACK diagonalizer with configuration parameters
    !> (e.g., which, bmat, nKry, print level) and binds ApplyMatOnVec.
    module subroutine Fabricate(this)
      !> The ARPACK diagonalizer instance to initialize
      class(T_DiagonalizerList_E_Arpack), intent(inout) :: this
    end subroutine
  end interface

  interface
    !> Performs setup operations for the ARPACK diagonalizer after initialization
    !> (e.g., size workspaces based on dim, nEvals, and nKry).
    module subroutine Setup(this)
      !> The ARPACK diagonalizer instance to set up
      class(T_DiagonalizerList_E_Arpack), intent(inout) :: this
    end subroutine
  end interface

  interface
    module subroutine Diagonalize(this, time, evecsQ)
      !> The ARPACK diagonalizer instance
      class(T_DiagonalizerList_E_Arpack), intent(inout) :: this
      !> Time at which to evaluate the operator if time dependent.
      real(R64), intent(in) :: time
      !> Flag indicating whether to compute eigenvectors
      logical, intent(in) :: evecsQ
    end subroutine
  end interface

end module
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
