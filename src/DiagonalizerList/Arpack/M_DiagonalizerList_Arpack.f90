! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file M_DiagonalizerList_Arpack.f90
!> @brief ARPACK-based iterative eigensolver backend.
!>
!> @details
!> This module provides an iterative eigensolver using ARPACK's Implicitly
!> Restarted Arnoldi Method (IRAM). It is suitable for large sparse or
!> matrix-free problems where only a few eigenpairs are needed.
!>
!> **Algorithm:**
!> The ARPACK backend uses the reverse-communication interface:
!> 1. Initialize Arnoldi iteration with random starting vector
!> 2. Build Krylov subspace of dimension `nKry` via repeated matvecs
!> 3. Extract Ritz values/vectors and restart until convergence
!> 4. Post-process to obtain final eigenpairs
!>
!> **Complexity:** O(nKry × dim × T_matvec + nKry³) per restart, typically
!> converges in O(nEvals/nKry) restarts for well-separated eigenvalues.
!>
!> **JSON Configuration:**
!> @code{.json}
!> "diagonalizerList": {
!>   "arpack_ground": {
!>     "nEvals": 5,
!>     "which": "SR",
!>     "bmat": "I",
!>     "nKry": 20,
!>     "checkConvergenceQ": true,
!>     "printLevel": 0
!>   }
!> }
!> @endcode
!>
!> | Parameter           | Type   | Default        | Description                         |
!> |---------------------|--------|----------------|-------------------------------------|
!> | `nEvals`            | int    | 1              | Number of eigenvalues (-1 = all)    |
!> | `which`             | string | "SR"           | Spectrum region (see below)         |
!> | `bmat`              | string | "I"            | "I" = standard, "G" = generalized   |
!> | `nKry`              | int    | 2*nEvals+1     | Krylov subspace dimension           |
!> | `tol`               | real   | 0 (machine eps)| Ritz value convergence tolerance    |
!> | `checkConvergenceQ` | bool   | true           | Verify residuals after solve        |
!> | `printLevel`        | int    | 0              | Verbosity level                     |
!>
!> **Spectrum selection (`which`):**
!> - "LM" / "SM" : Largest / Smallest Magnitude
!> - "LR" / "SR" : Largest / Smallest Real part (ground states)
!> - "LI" / "SI" : Largest / Smallest Imaginary part
!>
!> @see M_Utils_ArpackLib for underlying ARPACK wrappers
module M_DiagonalizerList_Arpack
  use M_Utils_Types
  use M_DiagonalizerList, only: T_DiagonalizerList_E

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> @brief Allocate an ARPACK diagonalizer instance.
    !> @param[out] e     Polymorphic pointer to allocated instance
    !> @param[in]  path  JSON configuration path
    module subroutine DiagonalizerList_Arpack_Allocate(e, path)
      class(T_DiagonalizerList_E), allocatable, intent(out) :: e
      character(len=*), intent(in) :: path
    end subroutine
  end interface

  !=============================================================================
  ! module types
  !=============================================================================

  !-----------------------------------------------------------------------------
  !> @brief ARPACK iterative eigensolver implementation.
  !>
  !> @details
  !> Extends the abstract `T_DiagonalizerList_E` to provide iterative
  !> Arnoldi/Lanczos diagonalization. Only requires matrix–vector products,
  !> never forms the explicit matrix.
  !>
  !> **Best suited for:**
  !> - Large sparse matrices (dim > 1000)
  !> - Matrix-free operators
  !> - Partial spectrum (nEvals << dim)
  !> - Ground state / low-lying excitations
  !-----------------------------------------------------------------------------
  type, extends(T_DiagonalizerList_E) :: T_DiagonalizerList_E_Arpack
    !> Spectrum region selector: 'LM','SM','LR','SR','LI','SI'
    character(2)  :: which = "SR"
    !> Problem type: 'I' = standard Ax=λx, 'G' = generalized Ax=λBx
    character(1)  :: bmat = "I"
    !> Krylov subspace dimension (number of Arnoldi vectors).
    !> Must satisfy: nEvals + 1 ≤ nKry ≤ dim.
    !> Larger values improve convergence but increase memory and compute cost.
    integer(I32)  :: nKry = 0
    !> Convergence tolerance for the Ritz values (0 => machine precision).
    !> Looser tolerances (e.g., 1e-12) can substantially reduce the number of
    !> restarts, especially for (nearly) degenerate eigenvalues.
    real(R64)     :: tol = 0.0_R64
    !> If true, compute and print residual norms after diagonalization.
    logical       :: checkConvergenceQ = .true.
  contains
    procedure :: Setup      !< Prepare for diagonalization
    procedure :: Fabricate  !< Read JSON parameters
    procedure :: Diagonalize !< Run ARPACK iteration
  end type

  interface
    !> @brief Read ARPACK parameters from JSON configuration.
    module subroutine Fabricate(this)
      class(T_DiagonalizerList_E_Arpack), intent(inout) :: this
    end subroutine
  end interface

  interface
    !> @brief Perform setup for ARPACK backend.
    module subroutine Setup(this)
      class(T_DiagonalizerList_E_Arpack), intent(inout) :: this
    end subroutine
  end interface

  interface
    !> @brief Diagonalize the operator using ARPACK.
    !>
    !> @details
    !> Runs the implicitly restarted Arnoldi method to compute `nEvals`
    !> eigenpairs in the spectral region specified by `which`. Uses a
    !> contained procedure to wrap `ApplyMatOnVec` for ARPACK's interface.
    !>
    !> @param[inout] this    Diagonalizer instance (results in %evals, %evecs)
    !> @param[in]    time    Time for time-dependent operators
    !> @param[in]    evecsQ  Whether to compute eigenvectors
    module subroutine Diagonalize(this, time, evecsQ)
      class(T_DiagonalizerList_E_Arpack), intent(inout) :: this
      real(R64), intent(in) :: time
      logical, intent(in) :: evecsQ
    end subroutine
  end interface

end module
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
