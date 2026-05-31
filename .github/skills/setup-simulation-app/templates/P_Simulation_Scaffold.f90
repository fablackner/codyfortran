program __PROGRAM_NAME__
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_PrinterObservable
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
  use M_Method_Mb_OrbBased
  use M_IntegratorList
  use M_Propagator

  implicit none

  character(len=*), parameter :: jsonFileName = "__JSON_REL_PATH__"
  real(R64) :: t, startTime, endTime, outputStep, timeStep
  integer(I32) :: iStep, nOutputSteps, nInnerSteps, iInner
  type(T_IntegratorList_FabricateInput) :: IntegratorListInput(1)

  call Say_Hello
  call Say_Section("start")

  call Json_LoadJsonFile(manualFileName_=jsonFileName, relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call SysKinetic_Fabricate
  call SysPotential_Fabricate
  call SysInteraction_Fabricate
  call Method_Fabricate
  call Orbs_Fabricate
  call OrbsInit_Fabricate
  call ConfigList_Fabricate
  call Coeffs_Fabricate
  call CoeffsInit_Fabricate
  IntegratorListInput(1) % TimeDerivative => Method_TimeDerivative
  call IntegratorList_Fabricate(IntegratorListInput)
  call Propagator_Fabricate

  call Say_Fabricate("program")
  startTime = Json_Get("program.startTime", 0.0_R64)
  endTime = Json_Get("program.endTime", 0.1_R64)
  outputStep = Json_Get("program.outputStep", 0.01_R64)
  nInnerSteps = Json_Get("program.nInnerSteps", 1)

  call Grid_Setup
  call SysKinetic_Setup
  call SysPotential_Setup
  call SysInteraction_Setup
  call Orbs_Setup
  call OrbsInit_Setup
  call ConfigList_Setup
  call Coeffs_Setup
  call CoeffsInit_Setup
  call Method_Setup
  call IntegratorList_Setup
  call Propagator_Setup

  call Say_Section("perform action")

  t = startTime
  nOutputSteps = max(1_I32, int((endTime - startTime) / outputStep + 0.5_R64))
  timeStep = outputStep / real(max(1_I32, nInnerSteps), R64)

  do iStep = 1, nOutputSteps
    do iInner = 1, nInnerSteps
      call Propagator_Propagate(Method_state, t, t + timeStep)
      t = t + timeStep
    end do

    call PrinterObservable_DumpEnergy("energy.dat", .true., t)
    call Coeffs_ApplyH1FillRdm1(Coeffs_coeffs, rdm1_=Method_Mb_OrbBased_rdm1)
    call PrinterObservable_DumpNorm(Method_Mb_OrbBased_rdm1, Method_Mb_OrbBased_rdm2, Orbs_orbs, "", .true.)
  end do

  call Say_Section("printout results")
  call Say_Goodbye

end program
