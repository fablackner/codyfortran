!> Adiabatic switching (method.mb.orbBased.adiabaticSwitching) with MCTDHF:
!> from the determinant of the first orbitals of the MCTDHF ground state, the
!> interaction switched on over -20 < t < 0, must arrive at the imaginary-time
!> MCTDHF ground state of T_He1d_04 (the switched Hamiltonian enters through
!> the shared operators of Method: the pair potentials, the one-body
!> reference, h2 of the CI part).
program T_He1d_14_AdiabaticSwitchingMctdhx
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
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
  use testdrive, only: check, error_type

  implicit none

  real(R64), parameter :: groundEnergy = -2.4852442075051160_R64
  real(R64) :: t, timeStep, energyStart, energy
  integer(I32) :: innerStep, nTimeSteps
  type(error_type), allocatable :: error
  type(T_IntegratorList_FabricateInput) :: IntegratorListInput(1)

  call Say_Hello
  call Json_LoadJsonFile(manualFileName_="test/simulationTests/He1d/T_He1d_14_AdiabaticSwitchingMctdhx.json", &
                         relativeToProjectDirQ_=.true.)

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
  call Propagator_Fabricate()

  nTimeSteps = Json_Get("program.nTimeSteps", 10)
  timeStep = Json_Get("program.outputStep", 0.5_R64) / nTimeSteps

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

  t = -Method_Mb_OrbBased_switchTime
  energyStart = Method_GetEnergy(t)
  do while (t < -0.5_R64 * timeStep)
    do innerStep = 1, nTimeSteps
      call Propagator_Propagate(Method_state, t, t + timeStep)
      t = t + timeStep
    end do
    write (*, '(A,F8.3,A,F8.5,A,F20.12)') " t =", t, "  lambda =", Method_Mb_OrbBased_SwitchFactor(t), &
      "  energy =", Method_GetEnergy(t)
  end do
  energy = Method_GetEnergy(0.0_R64)
  write (*, '(A,F20.12,A,F20.12,A,ES10.2)') " start (lambda = 0):", energyStart, "  t = 0:", energy, &
    "  minus the MCTDHF ground state:", energy - groundEnergy

  ! lambda at the ends and in the middle of the window
  call check(error, Method_Mb_OrbBased_SwitchFactor(-2.0_R64 * Method_Mb_OrbBased_switchTime), 0.0_R64, thr=1e-15_R64)
  if (allocated(error)) error stop "T_He1d_14: switching factor before the window"
  call check(error, Method_Mb_OrbBased_SwitchFactor(-0.5_R64 * Method_Mb_OrbBased_switchTime), 0.5_R64, thr=1e-14_R64)
  if (allocated(error)) error stop "T_He1d_14: switching factor in the middle of the window"
  call check(error, Method_Mb_OrbBased_SwitchFactor(0.0_R64), 1.0_R64, thr=1e-15_R64)
  if (allocated(error)) error stop "T_He1d_14: switching factor at t = 0"
  ! adiabatic: the ground state of the full Hamiltonian
  call check(error, energy, groundEnergy, thr=1e-5_R64)
  if (allocated(error)) error stop "T_He1d_14: the switched state is not the MCTDHF ground state"

  call Say_Goodbye

end program
