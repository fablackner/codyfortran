program T_H3d_01_ImagTimePropagationSb
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_Constants
  use M_Utils_DataStorage
  use M_Utils_Actions
  use M_Grid
  use M_SysKinetic
  use M_SysPotential
  use M_Orbs
  use M_OrbsInit
  use M_Method
  use M_IntegratorList
  use M_Propagator
  use testdrive, only: check, error_type

  implicit none

  real(R64)    :: energyNew = 1e10_R64
  real(R64)    :: conv = 1e10_R64
  real(R64)    :: energyOld
  real(R64)    :: convThresh
  real(R64)    :: outputStep
  real(R64)    :: time = 0.0_R64
  real(R64)    :: timeStep
  real(R64)    :: norm
  integer(I32) :: iStep, nTimeSteps, innerStep
  type(error_type), allocatable :: error
  character(len=:), allocatable :: jsonFileName
  type(T_IntegratorList_FabricateInput) :: IntegratorListInput(1)

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  jsonFileName = "test/T_H3d_01_ImagTimePropagationSingleBody.json"
  call Json_LoadJsonFile(manualFileName_=jsonFileName, relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call SysKinetic_Fabricate
  call SysPotential_Fabricate
  call OrbsInit_Fabricate
  call Method_Fabricate
  IntegratorListInput(1) % TimeDerivative => Actions_ImagTimeDerivative
  call IntegratorList_Fabricate(IntegratorListInput)
  call Propagator_Fabricate()

  call Say_Fabricate("program")
  convThresh = Json_Get("program.convThresh", 1e-13_R64)
  outputStep = Json_Get("program.outputStep", 1e-1_R64)
  nTimeSteps = Json_Get("program.nTimeSteps", 10)

  call Grid_Setup
  call SysKinetic_Setup
  call SysPotential_Setup
  call OrbsInit_Setup
  call Method_Setup
  call IntegratorList_Setup
  call Propagator_Setup

  !==========================================
  call Say_Section("perform action")
  !==========================================

  energyOld = Method_GetEnergy(time)  ! Use initial time

  iStep = 0
  time = 0.0_R64  ! Initialize simulation time
  timeStep = outputStep / nTimeSteps  ! Calculate step size for each iteration

  do while (abs(conv) > convThresh)
    iStep = iStep + 1

    ! Explicit loop to execute propagation nTimeSteps times
    do innerStep = 1, nTimeSteps
      call Propagator_Propagate(Method_state, time, time + timeStep)
      time = time + timeStep  ! Update current time

      norm = sqrt(real(Grid_InnerProduct(Method_state, Method_state), kind=R64))
      Method_state = Method_state / norm  ! Normalize the wavefunction
    end do

    energyOld = energyNew
    energyNew = Method_GetEnergy(time)  ! Use current time for energy calculation

    conv = energyNew - energyOld
    write (*, *) "convergence: ", energyNew, conv, " time:", time
  end do

  !==========================================
  call Say_Section("printout results")
  !==========================================

  print *
  print *, "final energy: ", energyNew

  call check(error, energyNew, -0.5_R64, thr=1e-5_R64)
  if (allocated(error)) error stop "P_Hydrogen_ImagTimePropagationSb failure"

  !==========================================
  call Say_Goodbye
  !==========================================

end program
