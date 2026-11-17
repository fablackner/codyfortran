!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> 1D helium: field-free real-time TD-2RDM (singlet, Valdemoro closure) of the
!> MCTDHF ground state written by T_He1d_04 (orb*.in, rdm2Bt*.in) after the
!> dipole kick exp(i kick x) of every orbital. Checks energy conservation
!> and the dipole <x> at the end time against the MCTDHF propagation
!> T_He1d_05_TimePropagationMctdhx (N = 2: the closure is the only
!> difference).
program T_He1d_06_TimePropagationTd2rdm
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_Constants
  use M_Printer_Observable_Linear
  use M_Printer_Density_Linear
  use M_Printer_Observable
  use M_Printer_Rdm
  use M_Grid
  use M_Grid_Linear
  use M_SysKinetic
  use M_SysPotential
  use M_SysInteraction
  use M_Orbs
  use M_OrbsInit
  use M_TwoRdm
  use M_TwoRdmInit
  use M_Method
  use M_Method_Mb_OrbBased
  use M_IntegratorList
  use M_Propagator
  use M_Purification
  use M_Reconstruction
  use M_Absorber
  use testdrive, only: check, error_type

  implicit none

  real(R64) :: t, startTime, endTime, outputStep
  integer(I32) :: iStep, nTimeSteps, innerStep, nOutputSteps, iOrb
  real(R64) :: startEnergy, energy, timeStep, kick, dipole
  type(error_type), allocatable :: error
  character(len=:), allocatable :: jsonFileName
  type(T_IntegratorList_FabricateInput) :: IntegratorListInput(1)

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  jsonFileName = "test/simulationTests/He1d/T_He1d_06_TimePropagationTd2rdm.json"
  call Json_LoadJsonFile(manualFileName_=jsonFileName, relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call SysKinetic_Fabricate
  call SysPotential_Fabricate
  call SysInteraction_Fabricate
  call Method_Fabricate
  call TwoRdm_Fabricate
  call TwoRdmInit_Fabricate
  call Orbs_Fabricate
  call OrbsInit_Fabricate
  IntegratorListInput(1) % TimeDerivative => Method_TimeDerivative
  call IntegratorList_Fabricate(IntegratorListInput)
  call Propagator_Fabricate()
  call Purification_Fabricate
  call Reconstruction_Fabricate
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
  call TwoRdm_Setup
  call TwoRdmInit_Setup
  call Orbs_Setup
  call OrbsInit_Setup
  call Method_Setup
  call IntegratorList_Setup
  call Propagator_Setup
  call Purification_Setup
  call Reconstruction_Setup
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
      call Absorber_ApplyAbsorber(Orbs_orbs)
      call Purification_PurifyTwoRdmState(TwoRdm_twoRdmState, Method_Mb_OrbBased_h2)
    end do

    energy = Method_GetEnergy(t)
    call TwoRdm_FillRdm1FromTwoRdmState(Method_Mb_OrbBased_rdm1, TwoRdm_twoRdmState)
    call TwoRdm_FillRdm2FromTwoRdmState(Method_Mb_OrbBased_rdm2, TwoRdm_twoRdmState)
    call Printer_Rdm_DumpRdm2Symmetries(Method_Mb_OrbBased_rdm2, "rdm2Symmetries.dat", .false.)
    call Printer_Observable_DumpNorm(Method_Mb_OrbBased_rdm1, Method_Mb_OrbBased_rdm2, Orbs_orbs, "", .true.)
    call Printer_Density_Linear_DumpOneBodyDensityOnGrid(Method_Mb_OrbBased_rdm1, Orbs_orbs, "density.dat", .false.)
    call Printer_Observable_Linear_DumpDipole(Method_Mb_OrbBased_rdm1, Orbs_orbs, "", .false., t, dipole_=dipole)
    write (*, '(A,F8.4,A,I0,A,F20.14,A,ES22.14)') ' time = ', t, '    iStep = ', iStep, '    energy = ', energy, &
      '    <x> = ', dipole
  end do

  !==========================================
  call Say_Section("printout results")
  !==========================================

  print *
  print *, "start energy: ", startEnergy
  print *, "final energy: ", energy
  print *, "final dipole: ", dipole

  ! The Valdemoro closure without ccQ does not conserve the energy (drift
  ! 6.5e-6 over 1 a.u., the same with dt = 0.001): pinned final energy
  call check(error, energy, -2.4297140051271939_R64, thr=1e-8_R64)
  if (allocated(error)) error stop "T_He1d_06_TimePropagationTd2rdm failure: energy not conserved"

  ! For N = 2 the exact 3-RDM vanishes and the closure error is the only
  ! difference to MCTDHF (measured 8.6e-5)
  call check(error, dipole, 0.32414036926390166_R64, thr=2e-4_R64)
  if (allocated(error)) error stop "T_He1d_06_TimePropagationTd2rdm failure: dipole off the MCTDHF reference"

  call TwoRdm_FillRdm1FromTwoRdmState(Method_Mb_OrbBased_rdm1, TwoRdm_twoRdmState)
  call Printer_Rdm_DumpRdm1(Method_Mb_OrbBased_rdm1, "rdm1.dat", .true.)
  call TwoRdm_FillRdm2FromTwoRdmState(Method_Mb_OrbBased_rdm2, TwoRdm_twoRdmState)
  call Printer_Rdm_DumpRdm2(Method_Mb_OrbBased_rdm2, "rdm2.dat", .false.)

  !==========================================
  call Say_Goodbye
  !==========================================

end program
