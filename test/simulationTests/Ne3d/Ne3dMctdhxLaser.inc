!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Helpers shared by the Ne MCTDHF laser tests T_Ne3d_07 .. T_Ne3d_16
! (included after `contains`).
!
! - SetupMctdhxLaser:           fabricate and set up the run of a JSON file
! - ConvergedGroundStateEnergy: MCSCF ground state (Ylm-optimized orbitals)
! - BeforeEndTimeQ, PropagateOutputStep: real-time propagation loop
! - MetricNorm, OtherMEntries:  checks of the propagated state
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Fabricates and sets up the run of jsonFileName: MCSCF ground solver with
!> one diagonalizer per l (lmax = 1), split-step propagator with the linear
!> part in integrator 1 and the non-linear part in integrator 2. The
!> one-body operator is self-adjoint only in the FEDVR-weighted metric, which
!> is supplied to integrator 1 (SIL): unit weights for the CI coefficients,
!> grid weights for each orbital.
subroutine SetupMctdhxLaser(jsonFileName)
  use M_Utils_Json
  use M_Utils_Say
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
  use M_IntegratorList
  use M_Propagator

  character(len=*), intent(in) :: jsonFileName

  type(T_DiagonalizerList_FabricateInput) :: DiagonalizerListInput(3)
  type(T_IntegratorList_FabricateInput) :: IntegratorListInput(2)
  real(R64), allocatable :: metricWeights(:)
  integer(I32) :: iOrb
  logical :: gaugeQ

  call Json_LoadJsonFile(manualFileName_=jsonFileName, relativeToProjectDirQ_=.true.)
  gaugeQ = Json_GetExistence("sysGauge")

  call Grid_Fabricate
  call SysKinetic_Fabricate
  call SysPotential_Fabricate
  if (gaugeQ) call SysGauge_Fabricate
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

  call Grid_Setup
  call SysKinetic_Setup
  call SysPotential_Setup
  if (gaugeQ) call SysGauge_Setup
  call SysInteraction_Setup
  call Orbs_Setup
  call OrbsInit_Setup
  call ConfigList_Setup
  call Coeffs_Setup
  call CoeffsInit_Setup
  call Method_Setup
  call GroundSolver_Setup
  call DiagonalizerList_Setup

  allocate (metricWeights(Coeffs_nCoeffs + Grid_nPoints * Orbs_nOrbsInState))
  metricWeights(1:Coeffs_nCoeffs) = 1.0_R64
  do iOrb = 1, Orbs_nOrbsInState
    metricWeights(Coeffs_nCoeffs + (iOrb - 1) * Grid_nPoints + 1: &
                  Coeffs_nCoeffs + iOrb * Grid_nPoints) = Grid_Ylm_weights(:)
  end do
  integratorList(1) % e % metricWeights = metricWeights

  call IntegratorList_Setup
  call Propagator_Setup

end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Converges the MCSCF ground state (program.nTimeSteps solver steps between
!> energy evaluations, until the energy changes by less than
!> program.convThresh) and returns its energy.
function ConvergedGroundStateEnergy() result(energyNew)
  use M_Utils_Json
  use M_Method
  use M_GroundSolver

  real(R64) :: energyNew

  real(R64) :: energyOld, conv, convThresh
  integer(I32) :: iStep, innerStep, nTimeSteps

  convThresh = Json_Get("program.convThresh", 1e-13_R64)
  nTimeSteps = Json_Get("program.nTimeSteps", 10)

  energyNew = 1e10_R64
  conv = 1e10_R64
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

end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Whether time is before program.endTime by more than half a time step.
function BeforeEndTimeQ(time) result(res)
  use M_Utils_Json

  real(R64), intent(in) :: time
  logical :: res

  real(R64) :: timeStep

  timeStep = Json_Get("program.outputStep", 2e-1_R64) / Json_Get("program.nTimeStepsPropagation", 4)
  res = time < Json_Get("program.endTime", 1.0_R64) - 0.5_R64 * timeStep

end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Propagates the state by one output step (program.outputStep in
!> program.nTimeStepsPropagation propagator steps), appends the dipole <z>
!> to dipoleFile and prints energy and dipole.
subroutine PropagateOutputStep(time, dipole, dipoleFile)
  use M_Utils_Json
  use M_Utils_PrinterObservableYlm
  use M_Orbs
  use M_Method
  use M_Method_Mb_OrbBased
  use M_Propagator

  real(R64), intent(inout) :: time
  real(R64), intent(out) :: dipole
  character(len=*), intent(in) :: dipoleFile

  complex(R64), allocatable :: rdm1(:, :)
  real(R64) :: timeStep, energy
  integer(I32) :: innerStep, nTimeStepsPropagation

  nTimeStepsPropagation = Json_Get("program.nTimeStepsPropagation", 4)
  timeStep = Json_Get("program.outputStep", 2e-1_R64) / nTimeStepsPropagation

  do innerStep = 1, nTimeStepsPropagation
    call Propagator_Propagate(Method_state, time, time + timeStep)
    time = time + timeStep
  end do

  energy = Method_GetEnergy(time)
  call Method_Mb_OrbBased_FillRdm1(rdm1, Method_state)
  call PrinterObservableYlm_DumpDipole(rdm1, Orbs_orbs, dipoleFile, .false., time, dipole_=dipole)

  write (*, '(A,ES10.3E1,A,E20.10,A,E20.10)') ' time = ', time, '    energy = ', energy, '    <z> = ', dipole

end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Norm of a packed state in the metric of integrator 1.
function MetricNorm(state) result(res)
  use M_IntegratorList

  complex(R64), intent(in), contiguous :: state(:)
  real(R64) :: res

  res = sqrt(sum(integratorList(1) % e % metricWeights(:) * abs(state(:))**2))

end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Number of nonzero entries of the orbitals in channels other than their
!> declared m.
function OtherMEntries(orbs) result(res)
  use M_Grid_Ylm
  use M_Method_Mb_OrbBased

  complex(R64), intent(in), contiguous :: orbs(:, :)
  integer(I32) :: res

  integer(I32) :: i

  res = 0
  do i = 1, size(orbs, 2)
    res = res + count(Grid_Ylm_mCoord .ne. Method_Mb_OrbBased_mOfOrb(i) .and. orbs(:, i) .ne. (0.0_R64, 0.0_R64))
  end do

end function

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine ApplyMatOnVecL0(dOrbLm, orbLm, time)
  use M_GroundSolver_Mcscf_YlmOpt

  complex(R64), intent(out), contiguous, target :: dOrbLm(:)
  complex(R64), intent(in), contiguous, target :: orbLm(:)
  real(R64), intent(in) :: time

  call GroundSolver_Mcscf_YlmOpt_FockAction(dOrbLm, orbLm, l=0, time=time)

end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine ApplyMatOnVecL1(dOrbLm, orbLm, time)
  use M_GroundSolver_Mcscf_YlmOpt

  complex(R64), intent(out), contiguous, target :: dOrbLm(:)
  complex(R64), intent(in), contiguous, target :: orbLm(:)
  real(R64), intent(in) :: time

  call GroundSolver_Mcscf_YlmOpt_FockAction(dOrbLm, orbLm, l=1, time=time)

end subroutine
