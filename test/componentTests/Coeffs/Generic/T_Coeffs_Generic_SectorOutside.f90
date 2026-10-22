!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> An initial CI vector outside the declared total-M sector: the unary
!> configuration (orbitals 1 2 of body type 1, 1 2 3 of body type 2, of
!> m = 0 0 and 0 0 -1, i.e. total M = -1) with coeffs.generic.totalM = 0.
!> The initialization must stop with the error of Coeffs_ProjectOnSector
!> instead of propagating a state the sector Hamiltonian does not act on:
!> ctest-pass-regex: Coeffs_ProjectOnSector: coefficients do not lie in the sector
program T_Coeffs_Generic_SectorOutside
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Grid
  use M_Method
  use M_ConfigList
  use M_Coeffs
  use M_CoeffsInit

  implicit none

  complex(R64), allocatable :: coeffs(:)

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call Json_LoadJsonFile(manualFileName_="test/componentTests/Coeffs/Generic/T_Coeffs_Generic_SectorOutside.json", &
                         relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call Method_Fabricate
  call ConfigList_Fabricate
  call Coeffs_Fabricate
  call CoeffsInit_Fabricate

  call Grid_Setup
  call ConfigList_Setup
  call Coeffs_Setup
  call CoeffsInit_Setup

  !==========================================
  call Say_Section("initialize (must stop with an error)")
  !==========================================

  allocate (coeffs(Coeffs_nCoeffs))
  call CoeffsInit_Initialize(coeffs)

  error stop "T_Coeffs_Generic_SectorOutside failure: a unary configuration outside the sector was accepted"

end program T_Coeffs_Generic_SectorOutside
