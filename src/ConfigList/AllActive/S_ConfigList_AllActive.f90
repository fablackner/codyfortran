! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_ConfigList_AllActive) S_ConfigList_AllActive

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine ConfigList_AllActive_Allocate(e, path)
    use M_Utils_Json
    use M_Utils_Say
    use M_ConfigList_AllActive_Fermionic, only: ConfigList_E_AllActive_Fermionic_Allocate
    use M_ConfigList_AllActive_Bosonic, only: ConfigList_E_AllActive_Bosonic_Allocate

    class(T_ConfigList_E), allocatable, intent(out) :: e
    character(len=*), intent(in) :: path

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence(path//".fermionic")) then
      call ConfigList_E_AllActive_Fermionic_Allocate(e, path//".fermionic")

    else if (Json_GetExistence(path//".bosonic")) then
      call ConfigList_E_AllActive_Bosonic_Allocate(e, path//".bosonic")

    else
      error stop path//". is missing one of: fermionic, bosonic"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Fabricate(this)
    use M_Utils_Say

    class(T_ConfigList_E_AllActive), intent(inout) :: this

    call Say_Fabricate(this % path)

    call this % FabricateLevel2()

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Setup(this)
    use M_Utils_Say

    class(T_ConfigList_E_AllActive), intent(inout) :: this

    call Say_Setup(this % path)

    call this % SetupLevel2()

    call SetupSinglesData(this)
    call SetupDoublesData(this)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine SetupSinglesData(this)
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    class(T_ConfigList_E_AllActive), intent(inout) :: this

    integer(I32) :: i1, j1, iC, iCNew
    integer(I32) :: nConnectedSinglesMax
    real(R64) :: factor
    integer(I32) :: k, bt, nOBt, nConfigurations, orbCode, shift, nO, nBBt

    bt = this % bodyTarget
    nOBt = Method_Mb_OrbBased_nOrbs(bt)
    nBBt = Method_Mb_nBodies(bt)
    nO = Method_Mb_OrbBased_nOrbsSum
    shift = sum(Method_Mb_OrbBased_nOrbs(1:bt - 1))
    nConfigurations = this % nConfigurations

    allocate (this % singles % nConnected(nConfigurations))

    !$omp parallel default(shared) private (iC, k, j1, i1, iCNew, factor)

    !$omp do
    do iC = 1, nConfigurations

      k = 0

      do j1 = 1, nOBt
        do i1 = 1, nOBt

          call this % ExciteConfiguration(iCNew, factor, [i1], [j1], iC)
          if (iCNew .eq. 0) cycle

          k = k + 1
        end do
      end do

      this % singles % nConnected(iC) = k

    end do
    !$omp end do

    !$omp end parallel

    nConnectedSinglesMax = maxval(this % singles % nConnected(:))

    allocate (this % singles % excitedC(nConnectedSinglesMax, nConfigurations))
    allocate (this % singles % orbCode(nConnectedSinglesMax, nConfigurations))
    allocate (this % singles % factor(nConnectedSinglesMax, nConfigurations))

    !$omp parallel default(shared) private (iC, k, j1, i1, iCNew, factor, orbCode)

    !$omp do
    do iC = 1, nConfigurations

      k = 0

      do j1 = 1, nOBt
        do i1 = 1, nOBt

          call this % ExciteConfiguration(iCNew, factor, [i1], [j1], iC)
          if (iCNew .eq. 0) cycle

          k = k + 1

          orbCode = 1
          orbCode = (orbCode - 1) * nO + (i1 + shift)
          orbCode = (orbCode - 1) * nO + (j1 + shift)

          this % singles % orbCode(k, iC) = orbCode
          this % singles % excitedC(k, iC) = iCNew
          this % singles % factor(k, iC) = factor

        end do
      end do
    end do
    !$omp end do

    !$omp end parallel

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine SetupDoublesData(this)
    use M_Method_Mb
    use M_Method_Mb_OrbBased

    class(T_ConfigList_E_AllActive), intent(inout) :: this

    integer(I32) :: i1, i2, j1, j2, iC, iCNew
    integer(I32) :: k, bt, nOBt, nConfigurations, orbCode, shift, nO, nBBt
    integer(I32) :: nConnectedDoublesMax
    real(R64) :: factor

    bt = this % bodyTarget
    nOBt = Method_Mb_OrbBased_nOrbs(bt)
    nBBt = Method_Mb_nBodies(bt)
    nO = Method_Mb_OrbBased_nOrbsSum
    shift = sum(Method_Mb_OrbBased_nOrbs(1:bt - 1))
    nConfigurations = this % nConfigurations

    allocate (this % doubles % nConnected(nConfigurations))

    !$omp parallel default(shared) private (iC, k, j2, i2, j1, i1, iCNew, factor)

    !$omp do
    do iC = 1, nConfigurations

      k = 0

      do j2 = 1, nOBt
        do j1 = 1, j2

          do i2 = 1, nOBt
            do i1 = 1, i2

              call this % ExciteConfiguration(iCNew, factor, [i1, i2], [j1, j2], iC)
              if (iCNew .eq. 0) cycle

              k = k + 1
            end do
          end do
        end do
      end do

      this % doubles % nConnected(iC) = k

    end do
    !$omp end do

    !$omp end parallel

    nConnectedDoublesMax = maxval(this % doubles % nConnected(:))

    allocate (this % doubles % excitedC(nConnectedDoublesMax, nConfigurations))
    allocate (this % doubles % orbCode(nConnectedDoublesMax, nConfigurations))
    allocate (this % doubles % factor(nConnectedDoublesMax, nConfigurations))

    !$omp parallel default(shared) private (iC, k, j2, i2, j1, i1, iCNew, factor, orbCode)

    !$omp do
    do iC = 1, nConfigurations

      k = 0

      do j2 = 1, nOBt
        do j1 = 1, j2

          do i2 = 1, nOBt
            do i1 = 1, i2

              call this % ExciteConfiguration(iCNew, factor, [i1, i2], [j1, j2], iC)
              if (iCNew .eq. 0) cycle

              k = k + 1

              orbCode = 1
              orbCode = (orbCode - 1) * nO + (i1 + shift)
              orbCode = (orbCode - 1) * nO + (j1 + shift)
              orbCode = (orbCode - 1) * nO + (i2 + shift)
              orbCode = (orbCode - 1) * nO + (j2 + shift)

              this % doubles % orbCode(k, iC) = orbCode
              this % doubles % excitedC(k, iC) = iCNew
              this % doubles % factor(k, iC) = factor

            end do
          end do
        end do
      end do
    end do
    !$omp end do

    !$omp end parallel

  end subroutine

end submodule
