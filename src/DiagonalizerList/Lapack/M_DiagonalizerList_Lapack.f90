! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Module M_DiagonalizerList_Lapack provides a LAPACK-based dense
!> eigensolver backend. It implements the abstract
!> T_DiagonalizerList_E interface using LAPACK routines to compute
!> eigenvalues and (optionally) eigenvectors of explicitly formed matrices.
!> The dense matrix may be assembled internally via the matrix–vector callback
!> if it is not already available in explicit form.
module M_DiagonalizerList_Lapack
  use M_Utils_Types
  use M_DiagonalizerList, only: T_DiagonalizerList_E

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Allocates a new LAPACK diagonalizer instance and associates it with a
    !> configuration path. Intended to be called by the list fabrication logic.
    module subroutine DiagonalizerList_Lapack_Allocate(e, path)
      class(T_DiagonalizerList_E), allocatable, intent(out) :: e
      character(len=*), intent(in) :: path
    end subroutine
  end interface

  !=============================================================================
  ! module types
  !=============================================================================

  !> LAPACK diagonalizer.
  !> Wraps dense eigensolvers suitable for small-to-medium sized problems where
  !> explicit matrices are available or can be assembled at reasonable cost.
  type, extends(T_DiagonalizerList_E) :: T_DiagonalizerList_E_Lapack
  contains
    !> Backend-specific setup (allocate work arrays, choose LAPACK routines).
    procedure :: Setup
    !> Read/derive LAPACK parameters from configuration and bind callbacks.
    procedure :: Fabricate
    !> Execute the LAPACK call(s) to compute eigenvalues (and optionally eigenvectors).
    procedure :: Diagonalize
  end type

  interface
    !> Initializes the LAPACK diagonalizer with configuration parameters and
    !> binds ApplyMatOnVec (which may be used to assemble the dense matrix).
    module subroutine Fabricate(this)
      !> The LAPACK diagonalizer instance to initialize
      class(T_DiagonalizerList_E_Lapack), intent(inout) :: this
    end subroutine
  end interface

  interface
    !> Performs setup operations for the LAPACK diagonalizer after initialization
    !> (e.g., allocate workspaces sized by dim and the requested problem type).
    module subroutine Setup(this)
      !> The LAPACK diagonalizer instance to set up
      class(T_DiagonalizerList_E_Lapack), intent(inout) :: this
    end subroutine
  end interface

  interface
    module subroutine Diagonalize(this, time, evecsQ)
      !> The LAPACK diagonalizer instance
      class(T_DiagonalizerList_E_Lapack), intent(inout) :: this
      !> Time (unused for time-independent dense problems; may be used if the
      !> dense matrix is assembled from a time-dependent operator).
      real(R64), intent(in) :: time
      !> Flag indicating whether to compute eigenvectors
      logical, intent(in) :: evecsQ
    end subroutine
  end interface

end module
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
