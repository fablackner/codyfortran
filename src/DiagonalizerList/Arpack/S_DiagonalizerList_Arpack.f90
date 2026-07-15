! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file S_DiagonalizerList_Arpack.f90
!> @brief Implementation of ARPACK iterative eigensolver backend.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_DiagonalizerList_Arpack) S_DiagonalizerList_Arpack

  implicit none

  !> @brief Module-level counter for matrix–vector products (debugging).
  !> @details Reset to 0 at start of each Diagonalize call when printLevel > 0.
  integer(I32), save :: count = 0

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Allocate a new ARPACK diagonalizer instance.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine DiagonalizerList_Arpack_Allocate(e, path)
    use M_Utils_UnusedVariables
    class(T_DiagonalizerList_E), allocatable, intent(out) :: e
    character(len=*), intent(in) :: path

    if (.false.) call UnusedVariables_Mark(path)

    allocate (T_DiagonalizerList_E_Arpack :: e)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Read configuration parameters for ARPACK backend.
!>
!> @details
!> Reads from JSON at `this%path`:
!> - `nEvals`: Number of eigenvalues (-1 = all, but defeats purpose of ARPACK)
!> - `which`: Spectrum region selector
!> - `bmat`: Problem type ('I' or 'G')
!> - `nKry`: Krylov dimension (default: min(2*nEvals+1, dim))
!> - `tol`: Ritz value convergence tolerance (default: 0 => machine precision)
!> - `checkConvergenceQ`: Post-solve residual verification
!> - `printLevel`: Verbosity level
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Fabricate(this)
    use M_Utils_Json
    use M_Utils_Say

    class(T_DiagonalizerList_E_Arpack), intent(inout) :: this

    call Say_Fabricate(this % path)

    !------------------------------------
    ! Read configuration parameters
    !------------------------------------

    this % nEvals = Json_Get("nEvals", 1, path_=this % path)
    if (this % nEvals < 0) then
      this % nEvals = this % dim  ! -1 means all (not recommended for ARPACK)
    end if

    this % which = Json_Get("which", "SR", path_=this % path)
    this % bmat = Json_Get("bmat", "I", path_=this % path)
    this % nKry = Json_Get("nKry", min(2 * this % nEvals + 1, this % dim), path_=this % path)
    this % tol = Json_Get("tol", 0.0_R64, path_=this % path)
    this % checkConvergenceQ = Json_Get("checkConvergenceQ", .true., path_=this % path)
    this % printLevel = Json_Get("printLevel", 0, path_=this % path)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Setup for ARPACK backend (currently minimal).
!>
!> @details
!> ARPACK's internal workspace is managed by ArpackLib. This routine logs
!> the setup call but does not pre-allocate resources.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Setup(this)
    use M_Utils_Say
    use M_Utils_UnusedVariables

    class(T_DiagonalizerList_E_Arpack), intent(inout) :: this

    call Say_Setup(this % path)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Run ARPACK iteration to compute eigenpairs.
!>
!> @details
!> **Algorithm:**
!> 1. Call `ArpackLib_Diagonalize` with wrapper callback
!> 2. Verify all requested eigenvalues converged
!> 3. Optionally verify residuals: ||Ax - λx||₂ for each eigenpair
!> 4. Print summary if printLevel > 0
!>
!> The wrapper callback `ApplyMatOnVec` (internal procedure) adapts the
!> 3-argument interface of `T_DiagonalizerList_E` to ARPACK's 2-argument
!> interface by capturing `time` from the enclosing scope.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Diagonalize(this, time, evecsQ)
    use M_Utils_ArpackLib
    use M_Utils_BlasLib

    class(T_DiagonalizerList_E_Arpack), intent(inout) :: this
    real(R64), intent(in) :: time
    logical, intent(in) :: evecsQ

    integer(I32) :: j
    real(R64) :: resnorm
    complex(R64), allocatable :: dState(:)

    if (this % printLevel > 0) then
      write (*, *)
      write (*, *) "ARPACK Diagonalization started."
      count = 0
    end if

    ! Call ARPACK library wrapper
    call ArpackLib_Diagonalize(this % evals, &
                               this % evecs, &
                               this % nFound, &
                               ApplyMatOnVec, &
                               this % dim, &
                               evecsQ, &
                               this % nEvals, &
                               this % which, &
                               this % bmat, &
                               this % nKry, &
                               tol_=this % tol)

    ! Verify convergence
    if (this % nFound .ne. this % nEvals) then
      write (*, *) "Warning not all eigenvalues are converged"
      error stop
    end if

    ! Optional residual check: ||H·ψ - E·ψ||
    if (evecsQ .and. this % checkConvergenceQ) then
      write (*, *)
      write (*, *) "checkConvergenceQ:"
      allocate (dState(this % dim))
      do j = 1, this % nFound
        call this % ApplyMatOnVec(dState, this % evecs(:, j), time)
        dState = dState - this % evals(j) * this % evecs(:, j)
        resnorm = BlasLib_CalcNorm(dState)
        write (*, *) "Eigenvector ", j, " residual 2-norm: ", resnorm
      end do
      write (*, *)
    end if

    ! Summary output
    if (this % printLevel > 0) then
      write (*, *)
      write (*, *) "Summary ARPACK:"
      write (*, *) "======================"
      write (*, *) "Size of the matrix is ", this % dim, "x", this % dim
      write (*, *) "The number of Ritz values requested is ", this % nEvals
      write (*, *) "The number of Arnoldi vectors generated (nKry) is ", this % nKry
      write (*, *) "What portion of the spectrum: ", this % which
      write (*, *) "The number of converged Ritz values is ", this % nFound
      write (*, *) "ARPACK Diagonalization finished."
      write (*, *)
    end if

  contains

    !--------------------------------------------------------------------------
    !> @brief Internal callback adapting 3-arg interface to ARPACK's 2-arg.
    !> @details Captures `time` from enclosing Diagonalize scope.
    !--------------------------------------------------------------------------
    subroutine ApplyMatOnVec(dState, state)

      complex(R64), intent(out), contiguous, target :: dState(:)
      complex(R64), intent(in), contiguous, target  :: state(:)

      if (this % printLevel > 0) then
        count = count + 1
        write (*, *) "ApplyMatOnVec count: ", count
      end if
      call this % ApplyMatOnVec(dState, state, time)

    end subroutine

  end subroutine

end submodule
