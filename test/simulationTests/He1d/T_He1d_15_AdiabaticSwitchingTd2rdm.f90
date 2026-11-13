!> TD-2RDM ground state by adiabatic switching (method.mb.orbBased.adiabaticSwitching,
!> moving orbitals, no CI): from the determinant of the first orbitals of the
!> MCTDHF ground state (twoRdmInit.excited), the interaction switched on over
!> -20 < t < 0. Checks the energy at t = 0 against its reference value and
!> its distance from the MCTDHF ground state (the reconstruction makes the
!> TD-2RDM state differ slightly).
program T_He1d_15_AdiabaticSwitchingTd2rdm
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Grid
  use M_SysKinetic
  use M_SysPotential
  use M_SysInteraction
  use M_Orbs
  use M_OrbsInit
  use M_TwoRdm
  use M_TwoRdmInit
  use M_Reconstruction
  use M_Method
  use M_Method_Mb_OrbBased
  use M_IntegratorList
  use M_Propagator
  use testdrive, only: check, error_type

  implicit none

  real(R64), parameter :: groundEnergy = -2.4852442075051160_R64
  real(R64) :: t, timeStep, energy
  integer(I32) :: innerStep, nTimeSteps
  type(error_type), allocatable :: error
  type(T_IntegratorList_FabricateInput) :: IntegratorListInput(1)

  call Say_Hello
  call Json_LoadJsonFile(manualFileName_="test/simulationTests/He1d/T_He1d_15_AdiabaticSwitchingTd2rdm.json", &
                         relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call SysKinetic_Fabricate
  call SysPotential_Fabricate
  call SysInteraction_Fabricate
  call Method_Fabricate
  call TwoRdm_Fabricate
  call TwoRdmInit_Fabricate
  call Orbs_Fabricate
  call OrbsInit_Fabricate
  call Reconstruction_Fabricate
  IntegratorListInput(1) % TimeDerivative => Method_TimeDerivative
  call IntegratorList_Fabricate(IntegratorListInput)
  call Propagator_Fabricate()

  nTimeSteps = Json_Get("program.nTimeSteps", 10)
  timeStep = Json_Get("program.outputStep", 0.5_R64) / nTimeSteps

  call Grid_Setup
  call SysKinetic_Setup
  call SysPotential_Setup
  call SysInteraction_Setup
  call TwoRdm_Setup
  call TwoRdmInit_Setup
  call Orbs_Setup
  call OrbsInit_Setup
  call Method_Setup
  call Reconstruction_Setup
  call IntegratorList_Setup
  call Propagator_Setup

  t = -Method_Mb_OrbBased_switchTime
  do while (t < -0.5_R64 * timeStep)
    do innerStep = 1, nTimeSteps
      call Propagator_Propagate(Method_state, t, t + timeStep)
      t = t + timeStep
    end do
    write (*, '(A,F8.3,A,F8.5,A,F20.12)') " t =", t, "  lambda =", Method_Mb_OrbBased_SwitchFactor(t), &
      "  energy =", Method_GetEnergy(t)
  end do
  energy = Method_GetEnergy(0.0_R64)
  write (*, '(A,F20.14,A,ES10.2)') " t = 0:", energy, "  minus the MCTDHF ground state:", energy - groundEnergy

  call check(error, energy, -2.48524382546435_R64, thr=1e-9_R64)
  if (allocated(error)) error stop "T_He1d_15: energy at t = 0 differs from the reference"
  call check(error, energy, groundEnergy, thr=1e-5_R64)
  if (allocated(error)) error stop "T_He1d_15: the switched state is far from the MCTDHF ground state"

  call Say_Goodbye

end program
