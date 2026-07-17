!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Real-time MCTDHF propagation of Ne (3d, Ylm grid) in a z-polarized
!> velocity-gauge laser pulse, computing the time-dependent dipole moment.
!>
!> Velocity-gauge companion of T_Ne3d_07_TimePropagationMctdhxLaser: the laser
!> enters through the SysGauge coupling A(t)·p + A²/2 instead of the
!> length-gauge E(t)·z potential, while the static -Z/r attraction remains an
!> ordinary SysPotential. The MCSCF ground state is converged first (A = 0 at
!> t = 0, so the gauge term vanishes and the ground state must match
!> T_Ne3d_06/07 exactly). The state is then propagated in real time through
!> the rising edge of a sin² pulse with the split-step (Strang) propagator:
!> the stiff linear one-body part (now including A·p) is integrated by the
!> Short Iterative Lanczos in the FEDVR-weighted metric, the nonlinear
!> mean-field/CI part by classical RK4.
!>
!> Checks:
!> - the MCSCF ground-state energy matches T_Ne3d_06/07 (gauge off at t = 0)
!> - the weighted-metric norm of the state is conserved by the propagation
!>   (requires A·p Hermitian in the weighted metric, see SysGauge_Ylm)
!> - the dipole moment <z> responds to the field (nonzero, regression value)
program T_Ne3d_08_TimePropagationMctdhxLaserVelocity
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_Actions
  use M_Utils_PrinterObservableYlm
  use M_Grid
  use M_Grid_Ylm
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
  use M_DiagonalizerList
  use M_Mixing
  use M_GroundSolver
  use M_GroundSolver_Mcscf
  use M_GroundSolver_Mcscf_YlmOpt
  use M_IntegratorList
  use M_Propagator
  use testdrive, only: check, error_type

  implicit none

  real(R64)    :: energyNew = 1e10_R64
  real(R64)    :: conv = 1e10_R64
  real(R64)    :: energyOld
  real(R64)    :: convThresh
  real(R64)    :: time, startTime, endTime, outputStep, timeStep
  real(R64)    :: startNorm, finalNorm
  real(R64)    :: dipole, energy
  integer(I32) :: iStep, nTimeSteps, nTimeStepsPropagation, innerStep, iOrb
  real(R64), allocatable :: metricWeights(:)
  complex(R64), allocatable :: rdm1(:, :)
  type(error_type), allocatable :: error
  character(len=:), allocatable :: jsonFileName
  type(T_DiagonalizerList_FabricateInput) :: DiagonalizerListInput(3)
  type(T_IntegratorList_FabricateInput) :: IntegratorListInput(2)

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  jsonFileName = "test/simulationTests/Ne3d/T_Ne3d_08_TimePropagationMctdhxLaserVelocity.json"
  call Json_LoadJsonFile(manualFileName_=jsonFileName, relativeToProjectDirQ_=.true.)

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
  call Mixing_Fabricate
  call GroundSolver_Fabricate
  DiagonalizerListInput(1) % ApplyMatOnVec => GroundSolver_Mcscf_HamiltonianAction
  DiagonalizerListInput(1) % dim = Coeffs_nCoeffs
  DiagonalizerListInput(2) % ApplyMatOnVec => ApplyMatOnVecL0
  DiagonalizerListInput(2) % dim = Grid_Ylm_nRadial
  DiagonalizerListInput(3) % ApplyMatOnVec => ApplyMatOnVecL1
  DiagonalizerListInput(3) % dim = Grid_Ylm_nRadial
  call DiagonalizerList_Fabricate(DiagonalizerListInput)
  IntegratorListInput(1) % TimeDerivative => Method_Mb_OrbBased_TimeDerivativeOrbsLin
  IntegratorListInput(2) % TimeDerivative => Method_Mb_OrbBased_TimeDerivativeCoeffsPlusOrbsNonLin
  call IntegratorList_Fabricate(IntegratorListInput)
  call Propagator_Fabricate()

  call Say_Fabricate("program")
  convThresh = Json_Get("program.convThresh", 1e-13_R64)
  nTimeSteps = Json_Get("program.nTimeSteps", 10)
  startTime = Json_Get("program.startTime", 0.0_R64)
  endTime = Json_Get("program.endTime", 1.0_R64)
  outputStep = Json_Get("program.outputStep", 2e-1_R64)
  nTimeStepsPropagation = Json_Get("program.nTimeStepsPropagation", 4)

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
  call GroundSolver_Setup
  call DiagonalizerList_Setup

  ! The one-body operator is self-adjoint only in the FEDVR-weighted metric;
  ! supply the metric to the SIL integrator (integrator 1): unit weights for
  ! the CI coefficients, grid weights for each orbital.
  allocate (metricWeights(Coeffs_nCoeffs + Grid_nPoints * Orbs_nOrbsInState))
  metricWeights(1:Coeffs_nCoeffs) = 1.0_R64
  do iOrb = 1, Orbs_nOrbsInState
    metricWeights(Coeffs_nCoeffs + (iOrb - 1) * Grid_nPoints + 1: &
                  Coeffs_nCoeffs + iOrb * Grid_nPoints) = Grid_Ylm_weights(:)
  end do
  integratorList(1) % e % metricWeights = metricWeights

  call IntegratorList_Setup
  call Propagator_Setup

  !==========================================
  call Say_Section("converge mcscf ground state")
  !==========================================

  iStep = 0

  do while (abs(conv) > convThresh)
    iStep = iStep + 1

    do innerStep = 1, nTimeSteps
      call GroundSolver_Approach(Method_state, time=0.0_R64)
    end do

    energyOld = energyNew
    energyNew = Method_GetEnergy(0.0_R64)

    conv = energyNew - energyOld
    write (*, *) "convergence: ", energyNew, conv, " step:", iStep
  end do

  print *
  print *, "ground-state energy: ", energyNew

  ! Gauge coupling is off at t = 0: must reproduce the T_Ne3d_06/07 mcscf energy.
  call check(error, energyNew, -127.7449896_R64, thr=1e-7_R64)
  if (allocated(error)) error stop "T_Ne3d_08_TimePropagationMctdhxLaserVelocity failure: ground state off reference"

  !==========================================
  call Say_Section("propagate through laser pulse")
  !==========================================

  startNorm = MetricNorm(Method_state)

  time = startTime
  iStep = 0
  timeStep = outputStep / nTimeStepsPropagation

  do while (time < endTime - 0.5_R64 * timeStep)
    iStep = iStep + 1

    do innerStep = 1, nTimeStepsPropagation
      call Propagator_Propagate(Method_state, time, time + timeStep)
      time = time + timeStep
    end do

    energy = Method_GetEnergy(time)
    call Method_Mb_OrbBased_FillRdm1(rdm1, Method_state)
    call PrinterObservableYlm_DumpDipole(rdm1, Orbs_orbs, "dipoleVelocity.dat", .false., time, dipole_=dipole)

    write (*, '(A,ES10.3E1,A,E20.10,A,E20.10)') ' time = ', time, '    energy = ', energy, '    <z> = ', dipole
  end do

  !==========================================
  call Say_Section("printout results")
  !==========================================

  finalNorm = MetricNorm(Method_state)

  print *
  print *, "start norm:   ", startNorm
  print *, "final norm:   ", finalNorm
  print *, "final dipole: ", dipole

  ! Unitary propagation must conserve the weighted-metric norm of the state
  ! (up to the RK4 integration error of the nonlinear part).
  call check(error, finalNorm, startNorm, thr=1e-7_R64)
  if (allocated(error)) error stop "T_Ne3d_08_TimePropagationMctdhxLaserVelocity failure: norm not conserved"

  ! The dipole moment must respond to the field (regression value). It is not
  ! comparable to T_Ne3d_07: for the 1-cycle pulse the envelope-derivative
  ! part of E = -dA/dt differs substantially from the length-gauge field.
  call check(error, dipole, 1.3142887110e-3_R64, thr=1e-8_R64)
  if (allocated(error)) error stop "T_Ne3d_08_TimePropagationMctdhxLaserVelocity failure: dipole off reference"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  subroutine ApplyMatOnVecL0(dOrbLm, orbLm, time)
    use M_GroundSolver_Mcscf_YlmOpt

    complex(R64), intent(out), contiguous, target :: dOrbLm(:)
    complex(R64), intent(in), contiguous, target :: orbLm(:)
    real(R64), intent(in) :: time

    call GroundSolver_Mcscf_YlmOpt_FockAction(dOrbLm, orbLm, l=0, time=time)

  end subroutine

  subroutine ApplyMatOnVecL1(dOrbLm, orbLm, time)
    use M_GroundSolver_Mcscf_YlmOpt

    complex(R64), intent(out), contiguous, target :: dOrbLm(:)
    complex(R64), intent(in), contiguous, target :: orbLm(:)
    real(R64), intent(in) :: time

    call GroundSolver_Mcscf_YlmOpt_FockAction(dOrbLm, orbLm, l=1, time=time)

  end subroutine

  function MetricNorm(state) result(res)
    complex(R64), intent(in), contiguous :: state(:)
    real(R64) :: res

    res = sqrt(sum(metricWeights(:) * abs(state(:))**2))

  end function

end program
