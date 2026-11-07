!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Partial load (orbsInit.load.nOrbsIn): the first orbitals come from files,
!> the remaining ones from the second backend, as for an MCSCF ground state
!> started from Hartree-Fock 1s 2s plus hydrogen-like 2p.
!>
!> Writes orb01_01.in and orb01_02.in of two s orbitals on the radial grid
!> only (an lmax = 0 Ylm layout), then initializes 5 restricted orbitals of
!> declared m = 0 0 -1 0 1 on an lmax = mmax = 1 grid through Method_Setup
!> with orbsInit.load (lmaxIn 0, nOrbsIn 2) and orbsInit.ylm.hydrogenLike
!> (n = 2, l = 1, m = -1 0 1, listing only the remaining orbitals). Checks:
!> - orbitals 1, 2: the written radial functions in the (0,0) channel
!>   (bitwise), zero elsewhere
!> - orbitals 3-5: the normalized hydrogen-like 2p r exp(-Z r / 2) in the
!>   channel (1, m) of their declared m, zero elsewhere
program T_OrbsInit_Load_Partial
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

  real(R64), parameter :: charge = 3.0_R64

  type(error_type), allocatable :: error
  character(len=:), allocatable :: jsonFileName
  complex(R64), allocatable :: radial(:, :), expected(:)
  character(len=256) :: filename
  integer(I32) :: i, iGrid, nWrong
  real(R64) :: dev

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  jsonFileName = "test/componentTests/OrbsInit/Load/T_OrbsInit_Load_Partial.json"
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
  call Say_Section("write two s orbitals (lmax 0 layout)")
  !==========================================

  allocate (radial(Grid_Ylm_nRadial, 2))
  radial(:, 1) = exp(-2.0_R64 * Grid_Ylm_radialPoints(:))
  radial(:, 2) = (1.0_R64 - Grid_Ylm_radialPoints(:)) * exp(-Grid_Ylm_radialPoints(:))
  do i = 1, 2
    write (filename, '("orb01_",I2.2,".in")') i
    call SaveData(radial(:, i), trim(filename), storage_size(radial(:, i)))
  end do

  !==========================================
  call Say_Section("initialize")
  !==========================================

  call Method_Setup

  ! The channel layout is l-major: the (0,0) channel is the prefix of length nRadial
  nWrong = 0
  do i = 1, 2
    nWrong = nWrong + count(abs(Orbs_orbs(1:Grid_Ylm_nRadial, i) - radial(:, i)) > 0.0_R64)
    nWrong = nWrong + count(abs(Orbs_orbs(Grid_Ylm_nRadial + 1:, i)) > 0.0_R64)
  end do
  print *, "loaded orbitals: wrong entries ", nWrong
  call check(error, nWrong, 0)
  if (allocated(error)) error stop "T_OrbsInit_Load_Partial failure: loaded orbitals"

  allocate (expected(Grid_nPoints))
  dev = 0.0_R64
  do i = 3, 5
    do iGrid = 1, Grid_nPoints
      expected(iGrid) = (0.0_R64, 0.0_R64)
      if (Grid_Ylm_lCoord(iGrid) .eq. 1 .and. Grid_Ylm_mCoord(iGrid) .eq. Method_Mb_OrbBased_mOfOrb(i)) &
        expected(iGrid) = Grid_Ylm_rCoord(iGrid) * exp(-charge * Grid_Ylm_rCoord(iGrid) / 2.0_R64)
    end do
    expected(:) = expected(:) / sqrt(real(Grid_InnerProduct(expected, expected), kind=R64))
    dev = max(dev, maxval(abs(Orbs_orbs(:, i) - expected)))
  end do
  print *, "hydrogen-like 2p: max deviation ", dev
  call check(error, dev < 1e-12_R64)
  if (allocated(error)) error stop "T_OrbsInit_Load_Partial failure: hydrogen-like orbitals"

  !==========================================
  call Say_Goodbye
  !==========================================

end program T_OrbsInit_Load_Partial
