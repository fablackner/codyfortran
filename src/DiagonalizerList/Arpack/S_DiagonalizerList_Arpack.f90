! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_DiagonalizerList_Arpack) S_DiagonalizerList_Arpack

  implicit none

  integer(I32), save :: count = 0

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine DiagonalizerList_Arpack_Allocate(e, path)
    use M_Utils_UnusedVariables
    class(T_DiagonalizerList_E), allocatable, intent(out) :: e
    character(len=*), intent(in) :: path

    if (.false.) call UnusedVariables_Mark(path)

    allocate (T_DiagonalizerList_E_Arpack :: e)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Fabricate(this)
    use M_Utils_Json
    use M_Utils_Say

    class(T_DiagonalizerList_E_Arpack), intent(inout) :: this

    call Say_Fabricate(this % path)

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    this % nEvals = Json_Get("nEvals", 1, path_=this % path)
    if (this % nEvals < 0) then
      this % nEvals = this % dim
    end if

    this % which = Json_Get("which", "SR", path_=this % path)
    this % bmat = Json_Get("bmat", "I", path_=this % path)
    this % nKry = Json_Get("nKry", min(2 * this % nEvals + 1, this % dim), path_=this % path)
    this % checkConvergenceQ = Json_Get("checkConvergenceQ", .true., path_=this % path)
    this % printLevel = Json_Get("printLevel", 0, path_=this % path)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Setup(this)
    use M_Utils_Say
    use M_Utils_UnusedVariables

    class(T_DiagonalizerList_E_Arpack), intent(inout) :: this

    call Say_Setup(this % path)

  end subroutine

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

    call ArpackLib_Diagonalize(this % evals, &
                               this % evecs, &
                               this % nFound, &
                               ApplyMatOnVec, &
                               this % dim, &
                               evecsQ, &
                               this % nEvals, &
                               this % which, &
                               this % bmat, &
                               this % nKry)

    if (this % nFound .ne. this % nEvals) then
      write (*, *) "Warning not all eigenvalues are converged"
      error stop
    end if

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
