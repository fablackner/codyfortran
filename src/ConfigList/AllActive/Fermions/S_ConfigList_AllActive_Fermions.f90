! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_ConfigList_AllActive_Fermionic) S_ConfigList_AllActive_Fermionic

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine ConfigList_E_AllActive_Fermionic_Allocate(e, path)
    use M_Utils_UnusedVariables

    class(T_ConfigList_E), allocatable, intent(out) :: e
    character(len=*), intent(in) :: path

    if (.false.) call UnusedVariables_Mark(path)

    allocate (T_ConfigList_E_AllActive_Fermionic :: e)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Fabricate(this)
    use M_Utils_Json
    use M_Utils_Say
    use M_Utils_SfGslLib
    use M_Method_Mb
    use M_Method_Mb_OrbBased
    use M_ConfigList

    class(T_ConfigList_E_AllActive_Fermionic), intent(inout) :: this

    integer(I32) :: bt, nOBt, nBBt, nE, i

    call Say_Fabricate(this % path//".fermionic")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    this % bodyTarget = Json_Get("bodyTarget", 1, path_=this % path//".fermionic")
    bt = this % bodyTarget

    nOBt = Method_Mb_OrbBased_nOrbs(bt)
    nBBt = Method_Mb_nBodies(bt)

    this % nExcitations = Json_Get("nExcitations", nBBt, path_=this % path//".fermionic")
    nE = this % nExcitations

    this % nConfigurations = 0
    do i = 0, nE
      this % nConfigurations = this % nConfigurations + SfGslLib_Binomial(nOBt - nBBt, i) * SfGslLib_Binomial(nBBt, i)
    end do

    if (Method_Mb_bodyStatistics(this % bodyTarget) .ne. 'f') error stop "bodyTarget not fermionic"

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Setup(this)
    use M_Utils_Say
    use M_Utils_Combinatorics
    use M_Utils_CombinationGslLib
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    class(T_ConfigList_E_AllActive_Fermionic), intent(inout) :: this

    integer(I32) :: iC, j, bt, nConfigurations, nBBt, nOBt
    integer(I32), allocatable :: i1(:, :)

    call Say_Setup(this % path//".fermionic")

    allocate (this % codeFromConfig(this % nConfigurations))

    bt = this % bodyTarget
    nOBt = Method_Mb_OrbBased_nOrbs(bt)
    nBBt = Method_Mb_nBodies(bt)
    nConfigurations = this % nConfigurations

    allocate (i1(nBBt, nConfigurations))

    !-------------------------------------------------------------------------------------
    ! i1(i, iC) is the orbital of the i-th body in config with index iC
    !-------------------------------------------------------------------------------------

    call CombinationGslLib_CombiNoRepeat(i1, nOBt)

    !-------------------------------------------------------------------------------------
    ! this%codeFromConfig(iC) encodes the distribution of bodies of this body type
    ! the binary digit of this%codeFromConfig(iC) is equal 1 for bodies and 0 for holes
    !-------------------------------------------------------------------------------------

    this % codeFromConfig = 0
    do j = 1, nBBt
      do iC = 1, nConfigurations
        this % codeFromConfig(iC) = ibset(this % codeFromConfig(iC), i1(j, iC) - 1)
      end do
    end do

    ! Verify correct index is recovered from config
    do iC = 1, nConfigurations
      if ((Combinatorics_indexOfCombiNoRepeat(nOBt, i1(:, iC)) - iC) .ne. 0) error stop "index calc failed"
    end do

    deallocate (i1)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine ExciteConfiguration(this, iCNew, factor, creates, destroys, iC)
    use M_Utils_Combinatorics
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    class(T_ConfigList_E_AllActive_Fermionic), intent(in) :: this
    integer(I32), intent(out)        :: iCNew
    real(R64), intent(out)           :: factor
    integer(I32), intent(in), contiguous          :: creates(:)
    integer(I32), intent(in), contiguous          :: destroys(:)
    integer(I32), intent(in)         :: iC

    integer(I32) :: nBBt, nOBt
    integer(I32) :: iDestroy, iCreate, i, j, sum

    integer(I64) :: code

    integer(I32), allocatable :: i1(:)
    integer(I32), allocatable :: occupation(:)

    nOBt = Method_Mb_OrbBased_nOrbs(this % bodyTarget)
    nBBt = Method_Mb_nBodies(this % bodyTarget)

    factor = 1.0_R64

    allocate (i1(nBBt))
    allocate (occupation(nOBt))

    code = this % codeFromConfig(iC)

    occupation = 0
    do i = 1, nOBt
      if (btest(code, i - 1)) then
        occupation(i) = 1
      end if
    end do

    do i = 1, size(destroys)
      iDestroy = destroys(i)

      if (occupation(iDestroy) .eq. 0) then
        iCNew = 0
        return
      end if

      do j = 1, iDestroy - 1
        if (occupation(j) .eq. 1) factor = (-1) * factor
      end do

      occupation(iDestroy) = occupation(iDestroy) - 1

    end do

    do i = 1, size(creates)
      iCreate = creates(size(creates) - i + 1)

      if (occupation(iCreate) .eq. 1) then
        iCNew = 0
        return
      end if

      do j = 1, iCreate - 1
        if (occupation(j) .eq. 1) factor = (-1) * factor
      end do

      occupation(iCreate) = occupation(iCreate) + 1

    end do

    i = 0
    sum = 1
    do while (sum <= nBBt)

      i = i + 1
      if (occupation(i) .eq. 0) cycle

      i1(sum) = i
      sum = sum + occupation(i)
    end do

    iCNew = Combinatorics_indexOfCombiNoRepeat(nOBt, i1(:))

  end subroutine

end submodule
