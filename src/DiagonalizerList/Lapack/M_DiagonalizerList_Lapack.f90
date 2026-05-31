! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file M_DiagonalizerList_Lapack.f90
!> @brief LAPACK-based dense eigensolver backend.
!>
!> @details
!> This module provides a dense matrix eigensolver using LAPACK routines
!> (specifically `ZHEEVR` for Hermitian matrices). It is suitable for
!> small-to-medium sized problems where explicit matrix construction is feasible.
!>
!> **Algorithm:**
!> 1. Assemble explicit matrix by applying `ApplyMatOnVec` to unit vectors
!> 2. Call LAPACK symmetric eigensolver with index range [1, nEvals]
!> 3. Return sorted eigenvalues and (optionally) eigenvectors
!>
!> **Complexity:** O(dim³) for full diagonalization, O(dim² · nEvals) for partial.
!>
!> **JSON Configuration:**
!> @code{.json}
!> "diagonalizerList": {
!>   "lapack_solver": {
!>     "nEvals": 10,
!>     "printLevel": 0
!>   }
!> }
!> @endcode
!>
!> | Parameter    | Type | Default | Description                        |
!> |--------------|------|---------|------------------------------------|
!> | `nEvals`     | int  | 1       | Number of eigenvalues (-1 = all)   |
!> | `printLevel` | int  | 0       | Verbosity (0=silent, 1=progress)   |
!>
!> @see M_Utils_LapackLib for underlying LAPACK wrappers
module M_DiagonalizerList_Lapack
  use M_Utils_Types
  use M_DiagonalizerList, only: T_DiagonalizerList_E

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> @brief Allocate a LAPACK diagonalizer instance.
    !> @param[out] e     Polymorphic pointer to allocated instance
    !> @param[in]  path  JSON configuration path
    module subroutine DiagonalizerList_Lapack_Allocate(e, path)
      class(T_DiagonalizerList_E), allocatable, intent(out) :: e
      character(len=*), intent(in) :: path
    end subroutine
  end interface

  !=============================================================================
  ! module types
  !=============================================================================

  !-----------------------------------------------------------------------------
  !> @brief LAPACK dense eigensolver implementation.
  !>
  !> @details
  !> Extends the abstract `T_DiagonalizerList_E` to provide dense matrix
  !> diagonalization. The matrix is assembled on-the-fly by probing with
  !> unit vectors through the `ApplyMatOnVec` callback.
  !>
  !> **Best suited for:**
  !> - Problems with dim ≲ 1000
  !> - Full spectrum calculations
  !> - Benchmarking/validation against iterative methods
  !-----------------------------------------------------------------------------
  type, extends(T_DiagonalizerList_E) :: T_DiagonalizerList_E_Lapack
  contains
    procedure :: Setup      !< Prepare for diagonalization (minimal for LAPACK)
    procedure :: Fabricate  !< Read JSON parameters
    procedure :: Diagonalize !< Build matrix and call LAPACK
  end type

  interface
    !> @brief Read LAPACK parameters from JSON configuration.
    module subroutine Fabricate(this)
      class(T_DiagonalizerList_E_Lapack), intent(inout) :: this
    end subroutine
  end interface

  interface
    !> @brief Perform setup (currently minimal for LAPACK backend).
    module subroutine Setup(this)
      class(T_DiagonalizerList_E_Lapack), intent(inout) :: this
    end subroutine
  end interface

  interface
    !> @brief Diagonalize the operator using LAPACK.
    !>
    !> @details
    !> Assembles the full dim×dim matrix by applying `ApplyMatOnVec` to
    !> each canonical basis vector, then calls `LapackLib_DiagonalizeSym`
    !> to compute the lowest `nEvals` eigenpairs.
    !>
    !> @note Matrix assembly has O(dim²) cost in memory and O(dim² × T_matvec)
    !>       time, where T_matvec is the cost of one matrix–vector product.
    !>
    !> @param[inout] this    Diagonalizer instance (results stored in %evals, %evecs)
    !> @param[in]    time    Time for time-dependent operators
    !> @param[in]    evecsQ  Whether to compute eigenvectors
    module subroutine Diagonalize(this, time, evecsQ)
      class(T_DiagonalizerList_E_Lapack), intent(inout) :: this
      real(R64), intent(in) :: time
      logical, intent(in) :: evecsQ
    end subroutine
  end interface

end module
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
