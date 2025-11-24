program T_FermiHubbard_03_ExactPropagationTdci
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_Actions
  use M_Utils_PrinterSpectrum
  use M_Grid
  use M_SysKinetic
  use M_SysPotential
  use M_SysInteraction
  use M_Orbs
  use M_OrbsInit
  use M_ConfigList
  use M_Coeffs
  use M_CoeffsInit
  use M_Method
  use M_DiagonalizerList
  use M_Propagator
  use testdrive, only: check, error_type

  implicit none

  real(R64) :: t, startTime, endTime, dt
  integer(I32) :: iStep
  real(R64) :: energy
  type(error_type), allocatable :: error
  character(len=:), allocatable :: jsonFileName
  type(T_DiagonalizerList_FabricateInput) :: DiagonalizerListInput(1)

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  jsonFileName = "test/T_FermiHubbard_03_ExactPropagationTdci.json"
  call Json_LoadJsonFile(manualFileName_=jsonFileName, relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call SysKinetic_Fabricate
  call SysInteraction_Fabricate
  call Method_Fabricate
  call Orbs_Fabricate
  call OrbsInit_Fabricate
  call ConfigList_Fabricate
  call Coeffs_Fabricate
  call CoeffsInit_Fabricate
  DiagonalizerListInput(1) % ApplyMatOnVec => Actions_ApplyHamiltonian
  DiagonalizerListInput(1) % dim = Coeffs_nCoeffs
  call DiagonalizerList_Fabricate(DiagonalizerListInput)
  call Propagator_Fabricate

  call Say_Fabricate("program")
  startTime = Json_Get("program.startTime", 0.0_R64)
  endTime = Json_Get("program.endTime", 1.0_R64)
  dt = Json_Get("program.dt", 1e-3_R64)

  call Grid_Setup
  call SysKinetic_Setup
  call SysInteraction_Setup
  call Orbs_Setup
  call OrbsInit_Setup
  call ConfigList_Setup
  call Coeffs_Setup
  call CoeffsInit_Setup
  call Method_Setup
  call DiagonalizerList_Setup
  call Propagator_Setup

  !==========================================
  call Say_Section("perform action")
  !==========================================

  call DiagonalizerList(1) % e % Diagonalize(0.0_R64, .true.)
  Method_state(:) = DiagonalizerList(1) % e % evecs(:, 1)

  t = startTime
  iStep = 0
  do while (t < endTime)
    write (*, '(A,ES10.3E1,A,I0,E20.10)') ' t = ', t, '    iStep = ', iStep

    call Propagator_Propagate(Method_state, t, t + dt)

    t = t + dt
    iStep = iStep + 1

  end do

  !==========================================
  call Say_Section("printout results")
  !==========================================

  energy = Method_GetEnergy(0.0_R64)

  write (*, *)
  write (*, *) "end energy: ", energy
  write (*, *)

  call check(error, energy, -7.29797344_R64, thr=1e-7_R64)
  if (allocated(error)) error stop "T_FermiHubbard_ExactPropagationTdci failure"

  !==========================================
  call Say_Goodbye
  !==========================================

end program
