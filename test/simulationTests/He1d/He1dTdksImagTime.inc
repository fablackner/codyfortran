!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Imaginary-time propagation of method.mb.orbBased.tdks on a 1D grid, shared by
! T_He1d_11 and T_He1d_13 (included after `contains`):
!
! - RunTdksImagTime: fabricate and set up the run of a JSON file (sysXc.linear),
!   propagate in imaginary time (Actions_ImagTimeDerivative, orthonormalized
!   after every step) until the Kohn-Sham energy changes by less than
!   program.convThresh per output step
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine RunTdksImagTime(energy, jsonFileName)
  use M_Utils_Json
  use M_Utils_Actions
  use M_Grid
  use M_SysKinetic
  use M_SysPotential
  use M_SysInteraction
  use M_SysXc
  use M_Orbs
  use M_OrbsInit
  use M_ConfigList
  use M_Coeffs
  use M_CoeffsInit
  use M_Method
  use M_IntegratorList
  use M_Propagator

  real(R64), intent(out) :: energy
  character(len=*), intent(in) :: jsonFileName

  type(T_IntegratorList_FabricateInput) :: IntegratorListInput(1)
  real(R64) :: energyOld, convThresh, outputStep, time, timeStep
  integer(I32) :: nTimeSteps, innerStep, iOut

  call Json_LoadJsonFile(manualFileName_=jsonFileName, relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call SysKinetic_Fabricate
  call SysPotential_Fabricate
  call SysInteraction_Fabricate
  call SysXc_Fabricate
  call Method_Fabricate
  call Orbs_Fabricate
  call OrbsInit_Fabricate
  call ConfigList_Fabricate
  call Coeffs_Fabricate
  call CoeffsInit_Fabricate
  IntegratorListInput(1) % TimeDerivative => Actions_ImagTimeDerivative
  call IntegratorList_Fabricate(IntegratorListInput)
  call Propagator_Fabricate()
  convThresh = Json_Get("program.convThresh", 1e-12_R64)
  outputStep = Json_Get("program.outputStep", 1.0_R64)
  nTimeSteps = Json_Get("program.nTimeSteps", 200)

  call Grid_Setup
  call SysKinetic_Setup
  call SysPotential_Setup
  call SysInteraction_Setup
  call SysXc_Setup
  call Orbs_Setup
  call OrbsInit_Setup
  call ConfigList_Setup
  call Coeffs_Setup
  call CoeffsInit_Setup
  call Method_Setup
  call IntegratorList_Setup
  call Propagator_Setup
  call Orbs_Orthonormalize(Orbs_orbs)

  time = 0.0_R64
  timeStep = outputStep / nTimeSteps
  energy = Method_GetEnergy(time)
  do iOut = 1, 10000
    do innerStep = 1, nTimeSteps
      call Propagator_Propagate(Method_state, time, time + timeStep)
      time = time + timeStep
      call Orbs_Orthonormalize(Orbs_orbs)
    end do
    energyOld = energy
    energy = Method_GetEnergy(time)
    write (*, '(1X, A, F10.3, A, F20.14, A, ES10.2)') "tau ", time, ": E ", energy, ", dE ", energy - energyOld
    if (abs(energy - energyOld) < convThresh) exit
  end do
  if (.not. (abs(energy - energyOld) < convThresh)) error stop "RunTdksImagTime (He1dTdksImagTime.inc): not converged"

end subroutine
