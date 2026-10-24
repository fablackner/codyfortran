!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> MCTDHF reference for HHG in Be (Lackner et al., PRA 95, 033414 (2017),
!> Fig. 4): real-time propagation in a z-polarized velocity-gauge few-cycle
!> pulse on a Ylm grid (m = 0 orbitals).
!>
!> The initial state (orbXX_YY.in, coeffs.in) is the imaginary-time MCTDHF
!> ground state on the same grid (relax/P_ImagTimeMctdhf). Split-step (Strang)
!> propagation: the stiff linear one-body part by integrator 1, the Short
!> Iterative Lanczos in the FEDVR metric; the interaction part by RK4. With
!> method.mb.orbBased.exactFieldFreeQ (the make_inputs.py default) the
!> field-free h0 (kinetic + centrifugal + Coulomb) is propagated exactly
!> (Method_PropagateExactPart) and composed with the SIL step according to
!> propagator.splitStep.exactPartScheme (make_inputs.py sets blanesMoan6; the
!> code default is strang), so SIL only
!> carries the A.p gauge coupling; otherwise SIL carries the whole linear
!> part. A radial cos^(1/order) mask absorbs outgoing flux after every time
!> step; with absorber.keepOrthonormalQ the orbitals are then
!> re-orthonormalized and the absorbed norm moves into the CI vector
!> (Method_Mb_OrbBased_OrthonormalizeKeepingState).
!>
!> Checkpoints: with program.saveStep > 0 (a multiple of program.outputStep)
!> the orbitals and the CI vector are saved (orbXX_YY.in, coeffs.in) after
!> every saveStep. A run is continued from the last checkpoint in the same
!> directory with orbsInit.load, coeffsInit.load and program.startTime set to
!> the checkpoint time; observables.dat and natocc.dat then keep their lines
!> before startTime and are continued.
!>
!> Output (one line per program.outputStep):
!>   observables.dat : time, E(t), A(t), <z>, dipole acceleration
!>                     -Z<z/r^3> (Eq. 38 of the paper), electrons left in the
!>                     box, energy
!>   natocc.dat      : time, natural occupation numbers (spin up)
program P_LaserMctdhf
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Grid
  use M_SysKinetic
  use M_SysPotential
  use M_SysGauge
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

  implicit none

  real(R64)    :: time, startTime, endTime, outputStep, timeStep, charge, saveStep
  integer(I32) :: nTimeSteps, innerStep, nOutputs, nOutputsPerSave
  type(T_IntegratorList_FabricateInput) :: IntegratorListInput(2)

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call Json_LoadJsonFile()

  call Grid_Fabricate
  call SysKinetic_Fabricate
  call SysPotential_Fabricate
  call SysGauge_Fabricate
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
  outputStep = Json_Get("program.outputStep", 0.2_R64)
  nTimeSteps = Json_Get("program.nTimeSteps", 20)
  saveStep = Json_Get("program.saveStep", 0.0_R64)
  charge = Json_Get("sysPotential.ylm.coulomb.charge", 4.0_R64)

  call Grid_Setup
  call SysKinetic_Setup
  call SysPotential_Setup
  call SysGauge_Setup
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

  ! The one-body operator is self-adjoint only in the FEDVR-weighted metric
  integratorList(1) % e % metricWeights = MetricWeights(Coeffs_nCoeffs)

  !==========================================
  call Say_Section("propagate through laser pulse")
  !==========================================

  time = startTime
  timeStep = outputStep / nTimeSteps
  nOutputsPerSave = nint(saveStep / outputStep)
  if (saveStep > 0.0_R64 .and. abs(nOutputsPerSave * outputStep - saveStep) > 1e-9_R64 * saveStep) &
    error stop "program.saveStep must be a multiple of program.outputStep"

  if (startTime > 0.0_R64) call TruncateLaserObservables(startTime - 0.5_R64 * timeStep)
  call WriteLaserObservables(time, charge, Method_GetEnergy(time), startTime <= 0.0_R64)

  nOutputs = 0
  do while (time < endTime - 0.5_R64 * timeStep)

    do innerStep = 1, nTimeSteps
      call Propagator_Propagate(Method_state, time, time + timeStep)
      time = time + timeStep
      call Absorber_ApplyAbsorber(Orbs_orbs)
      if (Absorber_keepOrthonormalQ) call Method_Mb_OrbBased_OrthonormalizeKeepingState(Method_state)
    end do

    call WriteLaserObservables(time, charge, Method_GetEnergy(time), .false.)

    nOutputs = nOutputs + 1
    if (nOutputsPerSave > 0) then
      if (mod(nOutputs, nOutputsPerSave) .eq. 0) then
        call Orbs_SaveOrbs(Orbs_orbs)
        call Coeffs_SaveCoeffs(Coeffs_coeffs)
        write (*, '(A,F9.3)') " checkpoint saved at t =", time
      end if
    end if

  end do

  !==========================================
  call Say_Section("save data")
  !==========================================

  call Orbs_SaveOrbs(Orbs_orbs)
  call Coeffs_SaveCoeffs(Coeffs_coeffs)

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  include '../common/BeCommon.inc'
  include '../common/BeObservables.inc'

end program
