!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Loading orbitals with declared m (method.mb.orbBased.mOfOrbs): the
!> initialization ends with Orbs_ProjectOnSymmetry, which removes round-off
!> admixtures of other m.
!>
!> Writes orb01_XX.in of 1s 2s 2p(-1) 2p(0) 2p(1)-like orbitals (restricted,
!> Ylm grid with lmax = mmax = 1) with an admixture of 1e-14 of the squared
!> norm in a channel of another m, loads them through Method_Setup
!> (orbsInit.load) and checks:
!> - the channels of other m are exactly zero
!> - the channels of the declared m are unchanged (bitwise)
program T_OrbsInit_Load_ProjectOnSymmetry
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_DataStorage
  use M_Grid
  use M_Grid_Ylm
  use M_Orbs
  use M_OrbsInit
  use M_ConfigList
  use M_Coeffs
  use M_CoeffsInit
  use M_Method
  use M_Method_Mb_OrbBased
  use testdrive, only: check, error_type

  implicit none

  integer(I32), parameter :: lOfOrb(5) = [0, 0, 1, 1, 1]

  type(error_type), allocatable :: error
  character(len=:), allocatable :: jsonFileName
  complex(R64), allocatable :: orbsWritten(:, :)
  character(len=256) :: filename
  integer(I32) :: i, nOtherM, nChanged

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  jsonFileName = "test/componentTests/OrbsInit/Load/T_OrbsInit_Load_ProjectOnSymmetry.json"
  call Json_LoadJsonFile(manualFileName_=jsonFileName, relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call Method_Fabricate
  call Orbs_Fabricate
  call OrbsInit_Fabricate
  call ConfigList_Fabricate
  call Coeffs_Fabricate
  call CoeffsInit_Fabricate

  call Grid_Setup
  call Orbs_Setup
  call OrbsInit_Setup
  call ConfigList_Setup
  call Coeffs_Setup
  call CoeffsInit_Setup

  !==========================================
  call Say_Section("write orbitals with admixtures of other m")
  !==========================================

  allocate (orbsWritten(Grid_nPoints, Orbs_nOrbsInState))
  do i = 1, Orbs_nOrbsInState
    call FillOrbital(orbsWritten(:, i), lOfOrb(i), Method_Mb_OrbBased_mOfOrb(i), 1.0e-7_R64)
    write (filename, '("orb01_",I2.2,".in")') i
    call SaveData(orbsWritten(:, i), trim(filename), storage_size(orbsWritten(:, i)))
  end do

  !==========================================
  call Say_Section("load")
  !==========================================

  call Method_Setup

  nOtherM = 0
  nChanged = 0
  do i = 1, Orbs_nOrbsInState
    nOtherM = nOtherM + count(Grid_Ylm_mCoord .ne. Method_Mb_OrbBased_mOfOrb(i) .and. &
                              Orbs_orbs(:, i) .ne. (0.0_R64, 0.0_R64))
    nChanged = nChanged + count(Grid_Ylm_mCoord .eq. Method_Mb_OrbBased_mOfOrb(i) .and. &
                                Orbs_orbs(:, i) .ne. orbsWritten(:, i))
  end do

  print *
  print *, "nonzero entries in channels of other m: ", nOtherM
  print *, "changed entries in channels of the declared m: ", nChanged

  call check(error, nOtherM, 0)
  if (allocated(error)) error stop "T_OrbsInit_Load_ProjectOnSymmetry failure: admixture not removed"

  call check(error, nChanged, 0)
  if (allocated(error)) error stop "T_OrbsInit_Load_ProjectOnSymmetry failure: declared channels changed"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  !> Normalized r^l exp(-r) in the channel (l, m) plus admixture * r exp(-r) in
  !> the channel (1, m') of another m.
  subroutine FillOrbital(orb, l, m, admixture)
    complex(R64), intent(out) :: orb(:)
    integer(I32), intent(in) :: l
    integer(I32), intent(in) :: m
    real(R64), intent(in) :: admixture

    integer(I32) :: iGrid, mOther

    mOther = merge(m - 1, m + 1, m .eq. 1)

    do iGrid = 1, Grid_nPoints
      orb(iGrid) = (0.0_R64, 0.0_R64)
      if (Grid_Ylm_lCoord(iGrid) .eq. l .and. Grid_Ylm_mCoord(iGrid) .eq. m) &
        orb(iGrid) = Grid_Ylm_rCoord(iGrid)**l * exp(-Grid_Ylm_rCoord(iGrid))
      if (Grid_Ylm_lCoord(iGrid) .eq. 1 .and. Grid_Ylm_mCoord(iGrid) .eq. mOther) &
        orb(iGrid) = cmplx(0.0_R64, admixture, kind=R64) * Grid_Ylm_rCoord(iGrid) * exp(-Grid_Ylm_rCoord(iGrid))
    end do
    orb(:) = orb(:) / sqrt(real(Grid_InnerProduct(orb, orb), kind=R64))

  end subroutine

end program T_OrbsInit_Load_ProjectOnSymmetry
