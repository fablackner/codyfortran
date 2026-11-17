!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Real-time MCTDHF propagation of He1d driven by the split-step propagator
!> with a Short Iterative Lanczos (SIL) integrator for the linear orbital
!> part and RK4 for the nonlinear part.
!>
!> Starts from the correlated ground state of T_He1d_04 (orb*.in and
!> coeffs.in) after the dipole kick exp(i kick x) of every orbital. The
!> field-free propagation must conserve the energy and the norm of the state,
!> and the dipole <x> at the end time must agree with the RK4 propagation of
!> T_He1d_05_TimePropagationMctdhx, since both integrations solve the same
!> split equations.
program T_He1d_06_TimePropagationMctdhxSil
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_Constants
  use M_Printer_Observable_Linear
  use M_Grid
  use M_Grid_Linear
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
  integer(I32) :: iStep, nTimeSteps, innerStep, nOutputSteps, iOrb
  real(R64) :: startEnergy, energy, timeStep, kick, dipole
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
  kick = Json_Get("program.kick", 0.0_R64)

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

  ! The dipole kick exp(i kick x) of every orbital makes the ground state
  ! non-stationary: a propagator that leaves the state unchanged fails the
  ! dipole check
  do iOrb = 1, Orbs_nOrbsInState
    Orbs_orbs(:, iOrb) = Orbs_orbs(:, iOrb) * exp(IU * kick * Grid_Linear_xCoord(:))
  end do

  startEnergy = Method_GetEnergy(startTime)
  startNorm = sqrt(sum(abs(Method_state)**2))
  write (*, *)
  write (*, *) "start energy: ", startEnergy
  write (*, *)

  t = startTime
  timeStep = outputStep / nTimeSteps
  nOutputSteps = nint((endTime - startTime) / outputStep)

  do iStep = 1, nOutputSteps
    do innerStep = 1, nTimeSteps
      call Propagator_Propagate(Method_state, t, t + timeStep)
      t = t + timeStep
    end do

    energy = Method_GetEnergy(t)
    call Coeffs_ApplyH1FillRdm1(Coeffs_coeffs, rdm1_=Method_Mb_OrbBased_rdm1)
    call Printer_Observable_Linear_DumpDipole(Method_Mb_OrbBased_rdm1, Orbs_orbs, "", .false., t, dipole_=dipole)
    write (*, '(A,F8.4,A,I0,A,F20.14,A,ES22.14)') ' time = ', t, '    iStep = ', iStep, '    energy = ', energy, &
      '    <x> = ', dipole
  end do

  !==========================================
  call Say_Section("printout results")
  !==========================================

  finalNorm = sqrt(sum(abs(Method_state)**2))
  print *
  print *, "start energy: ", startEnergy
  print *, "final energy: ", energy
  print *, "final dipole: ", dipole
  print *, "start norm:   ", startNorm
  print *, "final norm:   ", finalNorm

  ! Field-free propagation conserves the energy up to the splitting error of
  ! the order-2 split step (measured 3.1e-8 with dt = 0.01)
  call check(error, energy, startEnergy, thr=1e-7_R64)
  if (allocated(error)) error stop "T_He1d_06_TimePropagationMctdhxSil failure: energy not conserved"

  ! The dipole of the RK4 reference propagation (T_He1d_05_TimePropagationMctdhx),
  ! measured difference 1.2e-11
  call check(error, dipole, 0.32414036926390166_R64, thr=1e-9_R64)
  if (allocated(error)) error stop "T_He1d_06_TimePropagationMctdhxSil failure: dipole off the RK4 reference"

  ! Unitary propagation must conserve the norm of the full state vector.
  call check(error, finalNorm, startNorm, thr=1e-10_R64)
  if (allocated(error)) error stop "T_He1d_06_TimePropagationMctdhxSil failure: norm not conserved"

  !==========================================
  call Say_Goodbye
  !==========================================

end program
