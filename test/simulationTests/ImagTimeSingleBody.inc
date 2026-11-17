!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Imaginary-time ground state of the single-body method (method.sb), shared
! by the single-body tests of H2 and H3d; included after `contains`.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Loads jsonFileName and propagates Method_state with
!> Method_ImagTimeDerivative in output steps of program.outputStep
!> (program.nTimeSteps steps each) until the energy changes by less than
!> program.convThresh in one output step, normalizing the state after every
!> step. Stops unless the final energy agrees with referenceEnergy to
!> threshold.
subroutine CheckImagTimeSingleBody(jsonFileName, referenceEnergy, threshold)
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Grid
  use M_SysKinetic
  use M_SysPotential
  use M_OrbsInit
  use M_Method
  use M_IntegratorList
  use M_Propagator
  use testdrive, only: check, error_type

  character(len=*), intent(in) :: jsonFileName
  real(R64), intent(in) :: referenceEnergy
  real(R64), intent(in) :: threshold

  real(R64) :: energyNew, energyOld, conv, convThresh, outputStep, time, timeStep, norm
  integer(I32) :: iStep, nTimeSteps, innerStep
  type(error_type), allocatable :: error
  type(T_IntegratorList_FabricateInput) :: IntegratorListInput(1)

  !==========================================
  call Say_Section("start")
  !==========================================

  call Json_LoadJsonFile(manualFileName_=jsonFileName, relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call SysKinetic_Fabricate
  call SysPotential_Fabricate
  call OrbsInit_Fabricate
  call Method_Fabricate
  IntegratorListInput(1) % TimeDerivative => Method_ImagTimeDerivative
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

  iStep = 0
  time = 0.0_R64
  timeStep = outputStep / nTimeSteps
  energyNew = 1e10_R64
  conv = 1e10_R64

  do while (abs(conv) > convThresh)
    iStep = iStep + 1

    do innerStep = 1, nTimeSteps
      call Propagator_Propagate(Method_state, time, time + timeStep)
      time = time + timeStep

      norm = sqrt(real(Grid_InnerProduct(Method_state, Method_state), kind=R64))
      Method_state = Method_state / norm
    end do

    energyOld = energyNew
    energyNew = Method_GetEnergy(time)

    conv = energyNew - energyOld
    write (*, *) "convergence: ", energyNew, conv, " time:", time
  end do

  !==========================================
  call Say_Section("printout results")
  !==========================================

  print *
  print *, "final energy: ", energyNew

  call check(error, energyNew, referenceEnergy, thr=threshold)
  if (allocated(error)) error stop "CheckImagTimeSingleBody: energy off the reference"

end subroutine
