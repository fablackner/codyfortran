! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> High-level combinatorics utilities building on manual/GSL generators.
!>
!> Provides convenience entry points for combinations, permutations, and
!> multisets used by configuration enumeration code.
module M_Utils_Combinatorics
  use M_Utils_Types

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function Combinatorics_IndexOfCombiWithRepeat(n, combi) result(res)
    use M_Utils_SfGslLib

    integer(I32) :: res
    integer(I32), intent(in) :: n
    integer(I32), intent(in), contiguous  :: combi(:)

    integer(I32) :: i, k, s, e

    k = size(combi, 1)

    res = 1
    i = 1
    s = 1
    e = combi(i)
    res = res + SfGslLib_Binomial(n - s + k - i + 1, k - i + 1) - SfGslLib_Binomial(n - e + k - i + 1, k - i + 1)

    do i = 2, k
      s = combi(i - 1)
      e = combi(i)
      res = res + SfGslLib_Binomial(n - s + k - i + 1, k - i + 1) - SfGslLib_Binomial(n - e + k - i + 1, k - i + 1)
    end do

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function Combinatorics_IndexOfCombiNoRepeat(n, combi) result(res)
    use M_Utils_SfGslLib

    integer(I32) :: res
    integer(I32), intent(in) :: n
    integer(I32), intent(in), contiguous  :: combi(:)

    integer(I32) :: i, k, s, e

    k = size(combi, 1)

    res = 1
    i = 1
    s = 0
    e = combi(i)
    if (s + 1 < e) res = res + SfGslLib_Binomial(n - s, k - i + 1) - SfGslLib_Binomial(n - e + 1, k - i + 1)

    do i = 2, k

      s = combi(i - 1)
      e = combi(i)
      if (s + 1 < e) res = res + SfGslLib_Binomial(n - s, k - i + 1) - SfGslLib_Binomial(n - e + 1, k - i + 1)

      ! non-integrated version of the above expression
      !do j = s + 1, combi(i) - 1
      !  res = res + SfGslLib_Binomial(n - j, k - i)
      !end do

    end do

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function Combinatorics_IndexOfCombiNoRepeatSave(n, combi) result(res)
    use M_Utils_SfGslLib

    integer(I32) :: res
    integer(I32), intent(in) :: n
    integer(I32), intent(in), contiguous  :: combi(:)

    integer(I32) :: i, k, s, j

    k = size(combi, 1)

    res = 1

    i = 1
    s = 0

    do j = s + 1, combi(i) - 1
      res = res + SfGslLib_Binomial(n - j, k - i)
    end do

    do i = 2, k
      s = combi(i - 1)

      do j = s + 1, combi(i) - 1
        res = res + SfGslLib_Binomial(n - j, k - i)
      end do

    end do

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function PermutationGslLib_PermutationsSign(perm) result(res)

    integer(I32) :: res
    integer(I32), intent(in), contiguous  :: perm(:)

    integer(I32) :: i, nInversions

    nInversions = 0
    do i = 1, size(perm) - 1
      nInversions = nInversions + count(perm(i:) < perm(i))
    end do

    res = (-1)**nInversions

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Combinatorics_SortIntegerArray(array, sec, permut_)
    ! selfmade sort
    ! alternatives: call sortqq (loc(HubSys_hoppDN(:, 2)), nImportKinBodiesDN(2), SRT$INTEGER4) ! USE ifport
    !               call M01CBF (HubSys_hoppDN(:, 1), 1, nImportKinBodiesDN(1), 'a', ifail) ! NAG
    !               call FGSL_SORT_INT(HUBSYS_HOPPDN(:, 2) ,1_FGSL_SIZE_T, NIMPORTKINBODIESDN(2)) ! fgsl

    integer(I32), intent(inout), contiguous  :: array(:)
    integer(I32), intent(in)    :: sec
    integer(I32), intent(out), allocatable, optional :: permut_(:)

    integer(I32) :: mi, sav
    integer(I32) :: k

    if (present(permut_)) permut_ = [(k, k=1, sec)]

    if (sec <= 1) return ! ¯\_\u(ツ)_/¯

    k = 1
    do while (any(array(1:sec - 1) > array(2:sec)))

      mi = minloc(array(k:sec), dim=1) + k - 1

      sav = array(k)
      array(k) = array(mi)
      array(mi) = sav

      if (present(permut_)) then
        sav = permut_(k)
        permut_(k) = permut_(mi)
        permut_(mi) = sav
      end if

      k = k + 1
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Combinatorics_SortIntegerArrayabs(array, sec, permut_)

    integer(I32), intent(inout), contiguous  :: array(:)
    integer(I32), intent(in)    :: sec
    integer(I32), intent(out), allocatable, optional :: permut_(:)

    integer(I32) :: mi, sav
    integer(I32) :: k

    if (present(permut_)) permut_ = [(k, k=1, sec)]

    if (sec <= 1) return ! ¯\_\u(ツ)_/¯

    k = 1
    do while (any(abs(array(1:sec - 1)) > abs(array(2:sec))))

      mi = minloc(abs(array(k:sec)), dim=1) + k - 1

      sav = array(k)
      array(k) = array(mi)
      array(mi) = sav

      if (present(permut_)) then
        sav = permut_(k)
        permut_(k) = permut_(mi)
        permut_(mi) = sav
      end if

      k = k + 1
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Combinatorics_SortRealArray(array, sec, permut_)

    real(R64), intent(inout), contiguous  :: array(:)
    integer(I32), intent(in) :: sec
    integer(I32), intent(out), allocatable, optional :: permut_(:)

    real(R64) :: rSave
    integer(I32) :: iSave
    integer(I32) :: mi
    integer(I32) :: k

    if (present(permut_)) permut_ = [(k, k=1, sec)]

    if (sec <= 1) return ! ¯\_\u(ツ)_/¯

    k = 1
    do while (any(array(1:sec - 1) > array(2:sec)))

      mi = minloc(array(k:sec), dim=1) + k - 1

      rSave = array(k)
      array(k) = array(mi)
      array(mi) = rSave

      if (present(permut_)) then
        iSave = permut_(k)
        permut_(k) = permut_(mi)
        permut_(mi) = iSave
      end if

      k = k + 1
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  recursive subroutine Combinatorics_SortIntegerArrayabsRecursive(sec, array)

    integer(I32), intent(out), contiguous  :: array(:)
    integer(I32), intent(in)  :: sec

    integer(I32) :: val
    integer(I32) :: i

    if (sec <= 1) return ! ¯\_\u(ツ)_/¯

    call Combinatorics_SortIntegerArrayabs(array(1:sec - 1), sec - 1)

    do i = 1, sec - 1
      val = array(sec)
      if (abs(val) < abs(array(i))) then
        array(i + 1:sec) = array(i:sec - 1)

        array(i) = val
      end if
    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function Combinatorics_EqSet(array1, array2) result(res)

    logical :: res
    integer(I32), intent(in), contiguous  :: array1(:)
    integer(I32), intent(in), contiguous  :: array2(:)
    integer(I32) :: i

    res = .true.
    do i = 1, size(array1)
      if (any(array1(i) .eq. array2(:))) cycle
      res = .false.
    end do

  end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function Combinatorics_FermiSurface(index, nX, nY, nZ) result(res)
    integer(I32)                   :: res(3)
    integer(I32), intent(in)       :: index
    integer(I32), intent(in)       :: nX
    integer(I32), intent(in)       :: nY
    integer(I32), intent(in)       :: nZ

    integer(I32) :: ikx, iky, ikz, iCount
    integer(I32) :: ene(nX * nY * nZ)
    integer(I32) :: kxCoord(nX * nY * nZ)
    integer(I32) :: kyCoord(nX * nY * nZ)
    integer(I32) :: kzCoord(nX * nY * nZ)

    integer(I32), allocatable :: permut(:)

    ene(:) = 0
    kxCoord(:) = 0
    kyCoord(:) = 0
    kzCoord(:) = 0

    iCount = 0
    do ikz = 1, nZ
      do iky = 1, nY
        do ikx = 1, nX

          iCount = iCount + 1

          ene(iCount) = ikx**2 + iky**2 + ikz**2
          kxCoord(iCount) = ikx
          kyCoord(iCount) = iky
          kzCoord(iCount) = ikz

        end do
      end do
    end do

    call Combinatorics_SortIntegerArray(ene, iCount, permut_=permut)

    res = [kxCoord(permut(index)), kyCoord(permut(index)), kzCoord(permut(index))]

  end function

end module
