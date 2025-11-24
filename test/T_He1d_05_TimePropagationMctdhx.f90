program T_He1d_05_TimePropagationMctdhx
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_Actions
  use M_Utils_PrinterDensityLinear
  use M_Utils_PrinterObservable
  use M_Utils_PrinterObservableLinear
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
  real(R64) :: energy, timeStep
  type(error_type), allocatable :: error
  character(len=:), allocatable :: jsonFileName
  type(T_IntegratorList_FabricateInput) :: IntegratorListInput(2)
  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  jsonFileName = "test/T_He1d_05_TimePropagationMctdhx.json"
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

  energy = Method_GetEnergy(startTime)

  write (*, *)
  write (*, *) "start energy: ", energy
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

      call Absorber_ApplyAbsorber(Orbs_orbs)
    end do

    ! Print output every outputStep
    energy = Method_GetEnergy(t)
    call Coeffs_ApplyH1FillRdm1(Coeffs_coeffs, rdm1_=Method_Mb_OrbBased_rdm1)
    call Coeffs_ApplyH2FillRdm2(Coeffs_coeffs, rdm2_=Method_Mb_OrbBased_rdm2)
    call PrinterObservable_DumpNorm(Method_Mb_OrbBased_rdm1, Method_Mb_OrbBased_rdm2, Orbs_orbs, "", .true.)
    call PrinterDensityLinear_DumpOneBodyDensityOnGrid(Method_Mb_OrbBased_rdm1, Orbs_orbs, "density.dat", .false.)
    call PrinterObservableLinear_DumpDipole(Method_Mb_OrbBased_rdm1, Orbs_orbs, "dipole.dat", .true., t)

    write (*, '(A,ES10.3E1,A,I0,A,E20.10)') ' time = ', t, '    iStep = ', iStep, '    energy = ', energy
  end do

  !==========================================
  call Say_Section("printout results")
  !==========================================

  energy = Method_GetEnergy(t)

  print *
  print *, "final energy: ", energy

  call check(error, energy, -2.4852441374569345_R64, thr=1e-7_R64)
  if (allocated(error)) error stop "T_He1d_TimePropagationMctdhx failure"

  !==========================================
  call Say_Goodbye
  !==========================================

end program
