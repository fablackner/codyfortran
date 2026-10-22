!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Loading orbitals that do not have the declared m (method.mb.orbBased.mOfOrbs
!> in another order than the orbital files): the initialization must stop
!> with an error instead of propagating wrong physics.
!>
!> Writes orb01_XX.in of 1s 2s 2p(+1) 2p(0) 2p(-1)-like orbitals while the
!> configuration declares m = 0 0 -1 0 1, and loads them through Method_Setup
!> (orbsInit.load). The test passes only if the program stops with the error
!> of Orbs_ProjectOnSymmetry:
!> ctest-pass-regex: Orbs_ProjectOnSymmetry: orbital does not have the m declared
program T_OrbsInit_Load_WrongSymmetry
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

  implicit none

  integer(I32), parameter :: lOfOrb(5) = [0, 0, 1, 1, 1]
  !> m of the written orbitals: 2p(-1) and 2p(+1) swapped against the declaration
  integer(I32), parameter :: mOfFile(5) = [0, 0, 1, 0, -1]

  character(len=:), allocatable :: jsonFileName
  complex(R64), allocatable :: orb(:)
  character(len=256) :: filename
  integer(I32) :: i, iGrid

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  jsonFileName = "test/componentTests/OrbsInit/Load/T_OrbsInit_Load_WrongSymmetry.json"
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
  call Say_Section("write orbitals in another m order")
  !==========================================

  allocate (orb(Grid_nPoints))
  do i = 1, Orbs_nOrbsInState
    do iGrid = 1, Grid_nPoints
      orb(iGrid) = (0.0_R64, 0.0_R64)
      if (Grid_Ylm_lCoord(iGrid) .eq. lOfOrb(i) .and. Grid_Ylm_mCoord(iGrid) .eq. mOfFile(i)) &
        orb(iGrid) = Grid_Ylm_rCoord(iGrid)**lOfOrb(i) * exp(-Grid_Ylm_rCoord(iGrid))
    end do
    orb(:) = orb(:) / sqrt(real(Grid_InnerProduct(orb, orb), kind=R64))
    write (filename, '("orb01_",I2.2,".in")') i
    call SaveData(orb, trim(filename), storage_size(orb))
  end do

  !==========================================
  call Say_Section("load (must stop with an error)")
  !==========================================

  call Method_Setup

  error stop "T_OrbsInit_Load_WrongSymmetry failure: orbitals of the wrong m were accepted"

end program T_OrbsInit_Load_WrongSymmetry
