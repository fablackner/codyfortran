! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_ConfigList_AllActive_Bosonic) S_ConfigList_AllActive_Bosonic

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine ConfigList_E_AllActive_Bosonic_Allocate(e, path)
    use M_Utils_UnusedVariables
    class(T_ConfigList_E), allocatable, intent(out) :: e
    character(len=*), intent(in) :: path

    if (.false.) call UnusedVariables_Mark(path)

    allocate (T_ConfigList_E_AllActive_Bosonic :: e)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Fabricate(this)
    use M_Utils_Json
    use M_Utils_Say
    use M_Utils_SfGslLib
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    class(T_ConfigList_E_AllActive_Bosonic), intent(inout) :: this

    integer(I32) :: bt, nOBt, nBBt, nE

    call Say_Fabricate(this % path//".bosonic")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    this % bodyTarget = Json_Get("bodyTarget", 1, path_=this % path//".bosonic")
    bt = this % bodyTarget

    nOBt = Method_Mb_OrbBased_nOrbs(bt)
    nBBt = Method_Mb_nBodies(bt)

    this % nExcitations = Json_Get("nExcitations", nBBt, path_=this % path//".bosonic")
    nE = this % nExcitations

    this % nConfigurations = SfGslLib_Binomial(nOBt + nE - 1, nE)

    if (Method_Mb_bodyStatistics(this % bodyTarget) .ne. 'b') error stop "bodyTarget not bosonic"

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Setup(this)
    use M_Utils_Say
    use M_Utils_Combinatorics
    use M_Utils_MultisetGslLib
    use M_Utils_SfGslLib
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    class(T_ConfigList_E_AllActive_Bosonic), intent(inout) :: this

    integer(I32) :: iC, j, bt, nConfigurations, nOBt
    integer(I32), allocatable :: i1(:, :)

    integer(I32) :: nBBt

    call Say_Setup(this % path//".bosonic")

    allocate (this % codeFromConfig(this % nConfigurations))

    bt = this % bodyTarget
    nOBt = Method_Mb_OrbBased_nOrbs(bt)
    nBBt = Method_Mb_nBodies(bt)
    nConfigurations = this % nConfigurations

    allocate (i1(nBBt, nConfigurations))

    !-------------------------------------------------------------------------------------
    ! i1(i, iC) is the orbital of the i-th body in config with index iC
    !-------------------------------------------------------------------------------------

    call MultisetGslLib_CombiWithRepeat(i1, nOBt)

    !-------------------------------------------------------------------------------------
    ! this%codeFromConfig(iC) encodes the distribution of bodies of this body type
    ! the digit of this%codeFromConfig(iC) in base (nBBt+1) is equal to the occupation
    !-------------------------------------------------------------------------------------

    this % codeFromConfig(:) = 0
    do iC = 1, nConfigurations

      do j = 1, nBBt
        this % codeFromConfig(iC) = this % codeFromConfig(iC) + (nBBt + 1)**(i1(j, iC) - 1)
      end do

    end do

    !-------------------------------------------------------------------------------------
    ! because there are always nBBt bodies in the config the number representation
    ! is divisible by nBBt. We use this to reduce the number this%codeFromConfig.
    !-------------------------------------------------------------------------------------
    this % codeFromConfig(:) = this % codeFromConfig(:) / nBBt

    ! Verify correct index is recovered from config
    do iC = 1, nConfigurations
      if ((Combinatorics_indexOfCombiWithRepeat(nOBt, i1(:, iC)) - iC) .ne. 0) error stop "index calc failed"
    end do

    deallocate (i1)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine ExciteConfiguration(this, iCNew, factor, creates, destroys, iC)
    use M_Utils_Combinatorics
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    class(T_ConfigList_E_AllActive_Bosonic), intent(in) :: this
    integer(I32), intent(out)        :: iCNew
    real(R64), intent(out)           :: factor
    integer(I32), intent(in), contiguous          :: creates(:)
    integer(I32), intent(in), contiguous          :: destroys(:)
    integer(I32), intent(in)         :: iC

    integer(I32) :: nBBt, nOBt
    integer(I32) :: iFactor, base
    integer(I32) :: orb, i
    integer(I32) :: sum

    integer(I64) :: code
    integer(I64) :: codeTmp

    integer(I32), allocatable :: i1(:)
    integer(I32), allocatable :: occupation(:)

    nOBt = Method_Mb_OrbBased_nOrbs(this % bodyTarget)
    nBBt = Method_Mb_nBodies(this % bodyTarget)

    base = nBBt + 1
    iFactor = 1

    allocate (i1(nBBt))
    allocate (occupation(nOBt))

    code = this % codeFromConfig(iC)
    ! we multiply by nBBt to obtain the actual base (number) representations
    codeTmp = code * nBBt

    do i = 1, nOBt
      occupation(i) = int(Mod(codeTmp, int(base, kind=I64)), kind=I32)
      codeTmp = codeTmp / base
    end do

    do i = 1, size(destroys)
      orb = destroys(i)

      iFactor = iFactor * occupation(orb)

      if (iFactor .eq. 0) then
        iCNew = 0
        return
      end if

      occupation(orb) = occupation(orb) - 1

    end do

    do i = 1, size(creates)
      orb = creates(size(creates) - i + 1) ! is the order necessary

      iFactor = iFactor * (1 + occupation(orb))

      occupation(orb) = occupation(orb) + 1

    end do

    factor = sqrt(real(iFactor, R64))

    i = 0
    sum = 1
    do while (sum <= nBBt)

      i = i + 1
      if (occupation(i) .eq. 0) cycle

      i1(sum:sum + occupation(i) - 1) = i
      sum = sum + occupation(i)
    end do

    iCNew = Combinatorics_indexOfCombiWithRepeat(nOBt, i1)

  end subroutine

end submodule
