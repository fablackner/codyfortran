! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file S_DiagonalizerList_Lapack.f90
!> @brief Implementation of LAPACK dense eigensolver backend.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_DiagonalizerList_Lapack) S_DiagonalizerList_Lapack

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Allocate a new LAPACK diagonalizer instance.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine DiagonalizerList_Lapack_Allocate(e, path)
    use M_Utils_UnusedVariables
    class(T_DiagonalizerList_E), allocatable, intent(out) :: e
    character(len=*), intent(in) :: path

    if (.false.) call UnusedVariables_Mark(path)

    allocate (T_DiagonalizerList_E_Lapack :: e)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Read configuration parameters for LAPACK backend.
!>
!> @details
!> Reads from JSON at `this%path`:
!> - `nEvals`: Number of eigenvalues to compute (-1 means all)
!> - `printLevel`: Verbosity level
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Fabricate(this)
    use M_Utils_Json
    use M_Utils_Say

    class(T_DiagonalizerList_E_Lapack), intent(inout) :: this

    call Say_Fabricate(this % path)

    !------------------------------------
    ! Read configuration parameters
    !------------------------------------

    this % nEvals = Json_Get("nEvals", 1, path_=this % path)
    if (this % nEvals < 0) then
      this % nEvals = this % dim  ! -1 means compute all eigenvalues
    end if
    this % printLevel = Json_Get("printLevel", 0, path_=this % path)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Setup for LAPACK backend (currently no-op).
!>
!> @details
!> LAPACK backend allocates working arrays during Diagonalize, so Setup
!> is minimal. Future versions could pre-allocate reusable workspaces here.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Setup(this)
    use M_Utils_Say
    use M_Utils_UnusedVariables

    class(T_DiagonalizerList_E_Lapack), intent(inout) :: this

    call Say_Setup(this % path)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Diagonalize operator by explicit matrix assembly + LAPACK.
!>
!> @details
!> **Algorithm:**
!> 1. Allocate dim×dim Hermitian matrix H
!> 2. For each column i: compute H(:,i) = A · e_i (unit vector)
!> 3. Call LapackLib_DiagonalizeSym to get lowest nEvals eigenpairs
!>
!> @note The matrix is stored in row-major convention: H(i,:) = A·e_i,
!>       which for Hermitian operators is equivalent to column storage.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Diagonalize(this, time, evecsQ)
    use M_Utils_LapackLib

    class(T_DiagonalizerList_E_Lapack), intent(inout) :: this
    real(R64), intent(in) :: time
    logical, intent(in) :: evecsQ

    complex(R64), allocatable :: dState(:)    ! Result of A·x
    complex(R64), allocatable :: stateTmp(:)  ! Unit vector probe
    complex(R64), allocatable :: H(:, :)      ! Explicit matrix
    integer(I32)              :: i

    ! Allocate working arrays
    allocate (H(this % dim, this % dim))
    allocate (dState(this % dim))
    allocate (stateTmp(this % dim))

    if (this % printLevel > 0) write (*, *) "Lapack Diagonalization started."

    ! Assemble explicit matrix by probing with unit vectors
    do i = 1, this % dim
      if (this % printLevel > 0) write (*, *) "Column", i, "of", this % dim
      stateTmp(:) = 0.0_R64
      stateTmp(i) = 1.0_R64
      call this % ApplyMatOnVec(dState, stateTmp, time)
      H(i, :) = dState(:)
      stateTmp(i) = 0.0_R64
    end do

    ! Call LAPACK eigensolver for eigenvalue indices [1, nEvals]
    call LapackLib_DiagonalizeSym(this % evals, &
                                  this % evecs, &
                                  this % nFound, &
                                  H, &
                                  evecsQ, &
                                  imin_=1, &
                                  imax_=this % nEvals)

    if (this % printLevel > 0) write (*, *) "Lapack Diagonalization finished."
  end subroutine

end submodule
