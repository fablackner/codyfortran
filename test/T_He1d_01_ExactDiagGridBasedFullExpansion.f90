program T_He1d_01_ExactDiagGridBasedFullExpansion
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_Actions
  use M_Utils_PrinterSpectrum
  use M_Grid
  use M_SysKinetic
  use M_SysPotential
  use M_SysInteraction
  use M_Method
  use M_DiagonalizerList
  use testdrive, only: check, error_type

  implicit none

  type(error_type), allocatable :: error
  character(len=:), allocatable :: jsonFileName
  type(T_DiagonalizerList_FabricateInput) :: DiagonalizerListInput(1)

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  jsonFileName = "test/T_He1d_01_ExactDiagGridBasedFullExpansion.json"
  call Json_LoadJsonFile(manualFileName_=jsonFileName, relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call SysKinetic_Fabricate
  call SysPotential_Fabricate
  call SysInteraction_Fabricate
  call Method_Fabricate
  DiagonalizerListInput(1) % ApplyMatOnVec => Actions_ApplyHamiltonian
  DiagonalizerListInput(1) % dim = Grid_nPoints**2
  call DiagonalizerList_Fabricate(DiagonalizerListInput)

  call Grid_Setup
  call SysKinetic_Setup
  call SysPotential_Setup
  call SysInteraction_Setup
  call Method_Setup
  call DiagonalizerList_Setup

  !==========================================
  call Say_Section("perform action")
  !==========================================

  call DiagonalizerList(1) % e % Diagonalize(0.0_R64, .true.)

  !==========================================
  call Say_Section("printout results")
  !==========================================

  if (DiagonalizerList(1) % e % nFound .eq. 0) error stop "no eigenvalue found!"

  write (*, *)
  write (*, *) "groundstate energy: ", DiagonalizerList(1) % e % evals(1)
  write (*, *)

  call PrinterSpectrum_DumpSpectrum("eigenspectrum.dat", .true.)

  call check(error, DiagonalizerList(1) % e % evals(1), -2.4852443924542094_R64, thr=1e-7_R64)
  if (allocated(error)) error stop "T_He1d_ExactDiagFullExpansion failure"

  !==========================================
  call Say_Goodbye
  !==========================================

end program
