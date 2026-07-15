!> Real-time MCTDHF propagation of He1d driven by the split-step propagator
!> with a Short Iterative Lanczos (SIL) integrator for the linear orbital
!> part and RK4 for the nonlinear part.
!>
!> Starting from the correlated ground state (loaded from file), the
!> field-free propagation must conserve both the total energy and the norm
!> of the state. The final energy is additionally checked against the
!> reference value of the RK4-only propagation test
!> (T_He1d_05_TimePropagationMctdhx), since both integrations solve the
!> same split equations to high accuracy.
program T_He1d_06_TimePropagationMctdhxSil
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_Actions
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
  use M_Absorber
  use testdrive, only: check, error_type

  implicit none

  real(R64) :: t, startTime, endTime, outputStep
  integer(I32) :: iStep, nTimeSteps, innerStep
  real(R64) :: startEnergy, energy, timeStep
  real(R64) :: startNorm, finalNorm
  type(error_type), allocatable :: error
  character(len=:), allocatable :: jsonFileName
  type(T_IntegratorList_FabricateInput) :: IntegratorListInput(2)
  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  jsonFileName = "test/simulationTests/He1d/T_He1d_06_TimePropagationMctdhxSil.json"
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
  IntegratorListInput(1) % TimeDerivative => Method_Mb_OrbBased_TimeDerivativeOrbsLin
  IntegratorListInput(2) % TimeDerivative => Method_Mb_OrbBased_TimeDerivativeCoeffsPlusOrbsNonLin
  call IntegratorList_Fabricate(IntegratorListInput)
  call Propagator_Fabricate()
  call Absorber_Fabricate

  call Say_Fabricate("program")
  startTime = Json_Get("program.startTime", 0.0_R64)
  endTime = Json_Get("program.endTime", 1.0_R64)
  outputStep = Json_Get("program.outputStep", 2e-5_R64)
  nTimeSteps = Json_Get("program.nTimeSteps", 2)

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
  call Absorber_Setup

  !==========================================
  call Say_Section("perform action")
  !==========================================

  startEnergy = Method_GetEnergy(startTime)
  startNorm = sqrt(sum(abs(Method_state)**2))

  write (*, *)
  write (*, *) "start energy: ", startEnergy
  write (*, *)

  t = startTime
  iStep = 0
  timeStep = outputStep / nTimeSteps  ! Calculate step size for each iteration

  do while (t < endTime)
    iStep = iStep + 1

    ! Explicit loop to execute propagation nTimeSteps times
    do innerStep = 1, nTimeSteps
      call Propagator_Propagate(Method_state, t, t + timeStep)
      t = t + timeStep  ! Update current time
    end do

    ! Print output every outputStep
    energy = Method_GetEnergy(t)

    write (*, '(A,ES10.3E1,A,I0,A,E20.10)') ' time = ', t, '    iStep = ', iStep, '    energy = ', energy
  end do

  !==========================================
  call Say_Section("printout results")
  !==========================================

  energy = Method_GetEnergy(t)
  finalNorm = sqrt(sum(abs(Method_state)**2))

  print *
  print *, "start energy: ", startEnergy
  print *, "final energy: ", energy
  print *, "start norm:   ", startNorm
  print *, "final norm:   ", finalNorm

  ! Field-free propagation must conserve the energy up to the O(dt^2)
  ! splitting error of the order-2 split-step propagator.
  call check(error, energy, startEnergy, thr=1e-6_R64)
  if (allocated(error)) error stop "T_He1d_06_TimePropagationMctdhxSil failure: energy not conserved"

  ! The SIL result must agree with the RK4 reference propagation
  ! (T_He1d_05_TimePropagationMctdhx).
  call check(error, energy, -2.4852441374569345_R64, thr=1e-7_R64)
  if (allocated(error)) error stop "T_He1d_06_TimePropagationMctdhxSil failure: final energy off reference"

  ! Unitary propagation must conserve the norm of the full state vector.
  call check(error, finalNorm, startNorm, thr=1e-10_R64)
  if (allocated(error)) error stop "T_He1d_06_TimePropagationMctdhxSil failure: norm not conserved"

  !==========================================
  call Say_Goodbye
  !==========================================

end program
