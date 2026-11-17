!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> 1D helium: field-free real-time MCTDHF (split step, order 2, RK4 for both
!> parts) of the ground state written by T_He1d_04 (orb*.in, coeffs.in) after
!> the dipole kick exp(i kick x) of every orbital. Checks energy conservation
!> and the dipole <x> at the end time, which is also the reference of
!> T_He1d_06_TimePropagationMctdhxSil and T_He1d_06_TimePropagationTd2rdm.
program T_He1d_05_TimePropagationMctdhx
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_Constants
  use M_Printer_Density_Linear
  use M_Printer_Observable
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
  type(error_type), allocatable :: error
  character(len=:), allocatable :: jsonFileName
  type(T_IntegratorList_FabricateInput) :: IntegratorListInput(2)
  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  jsonFileName = "test/simulationTests/He1d/T_He1d_05_TimePropagationMctdhx.json"
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
    end do

    energy = Method_GetEnergy(t)
    call Coeffs_ApplyH1FillRdm1(Coeffs_coeffs, rdm1_=Method_Mb_OrbBased_rdm1)
    call Coeffs_ApplyH2FillRdm2(Coeffs_coeffs, rdm2_=Method_Mb_OrbBased_rdm2)
    call Printer_Observable_DumpNorm(Method_Mb_OrbBased_rdm1, Method_Mb_OrbBased_rdm2, Orbs_orbs, "", .true.)
    call Printer_Density_Linear_DumpOneBodyDensityOnGrid(Method_Mb_OrbBased_rdm1, Orbs_orbs, "density.dat", .false.)
    call Printer_Observable_Linear_DumpDipole(Method_Mb_OrbBased_rdm1, Orbs_orbs, "dipole.dat", .false., t, dipole_=dipole)
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

  ! Field-free propagation conserves the energy up to the splitting error of
  ! the order-2 split step (measured 3.1e-8 with dt = 0.01)
  call check(error, energy, startEnergy, thr=1e-7_R64)
  if (allocated(error)) error stop "T_He1d_05_TimePropagationMctdhx failure: energy not conserved"

  ! <x> after 1 a.u. (the free flight of the kick alone gives 0.36; the dt
  ! error of the split step is 1e-7)
  call check(error, dipole, 0.32414036926390166_R64, thr=1e-9_R64)
  if (allocated(error)) error stop "T_He1d_05_TimePropagationMctdhx failure: dipole off reference"

  !==========================================
  call Say_Goodbye
  !==========================================

end program
