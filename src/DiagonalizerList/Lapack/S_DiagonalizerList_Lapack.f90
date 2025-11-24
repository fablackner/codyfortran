! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_DiagonalizerList_Lapack) S_DiagonalizerList_Lapack

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine DiagonalizerList_Lapack_Allocate(e, path)
    use M_Utils_UnusedVariables
    class(T_DiagonalizerList_E), allocatable, intent(out) :: e
    character(len=*), intent(in) :: path

    if (.false.) call UnusedVariables_Mark(path)

    allocate (T_DiagonalizerList_E_Lapack :: e)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Fabricate(this)
    use M_Utils_Json
    use M_Utils_Say

    class(T_DiagonalizerList_E_Lapack), intent(inout) :: this

    call Say_Fabricate(this % path)

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    this % nEvals = Json_Get("nEvals", 1, path_=this % path)
    if (this % nEvals < 0) then
      this % nEvals = this % dim
    end if
    this % printLevel = Json_Get("printLevel", 0, path_=this % path)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Setup(this)
    use M_Utils_Say
    use M_Utils_UnusedVariables

    class(T_DiagonalizerList_E_Lapack), intent(inout) :: this

    call Say_Setup(this % path)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Diagonalize(this, time, evecsQ)
    use M_Utils_LapackLib

    class(T_DiagonalizerList_E_Lapack), intent(inout) :: this
    real(R64), intent(in) :: time
    logical, intent(in) :: evecsQ

    complex(R64), allocatable :: dState(:)
    complex(R64), allocatable :: stateTmp(:)
    complex(R64), allocatable :: H(:, :)
    integer(I32)              :: i

    allocate (H(this % dim, this % dim))
    allocate (dState(this % dim))
    allocate (stateTmp(this % dim))

    if (this % printLevel > 0) write (*, *) "Lapack Diagonalization started."

    do i = 1, this % dim
      if (this % printLevel > 0) write (*, *) "Column", i, "of", this % dim
      stateTmp(:) = 0.0_R64
      stateTmp(i) = 1.0_R64
      call this % ApplyMatOnVec(dState, stateTmp, time)
      H(i, :) = dState(:)
      stateTmp(i) = 0.0_R64
    end do

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
