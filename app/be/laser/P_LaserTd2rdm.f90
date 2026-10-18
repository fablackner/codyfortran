!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> TD-2RDM for HHG in Be (Lackner et al., PRA 95, 033414 (2017), Fig. 4):
!> real-time propagation of the (up,down) block of the two-particle reduced
!> density matrix together with time-dependent orbitals (MCTDHF orbital
!> equations) in a z-polarized velocity-gauge few-cycle pulse.
!>
!> The 3-RDM in the collision operator is reconstructed with the contraction-
!> consistent Nakatsuji-Yasuda functional. Optionally (purification block in
!> the input) the D- and Q-conditions are enforced after every time step as
!> in the paper (Eqs. 32-34): the contraction-free kernel of the negative
!> parts of D and Q is subtracted, which keeps norm and 1-RDM but not the
!> energy. The initial state (orbXX_YY.in, rdm2BtXX_YY.in) is the time-averaged MCTDHF
!> ground state (td2rdmInit/P_Td2rdmTimeAverage, appendix C of the paper).
!> Grid, pulse, absorber and propagation scheme are identical to
!> P_LaserMctdhf.
!>
!> Output: observables.dat and natocc.dat as in P_LaserMctdhf.
program P_LaserTd2rdm
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Grid
  use M_Grid_Ylm
  use M_SysKinetic
  use M_SysPotential
  use M_SysGauge
  use M_SysGauge_Ylm_Velocity
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

  implicit none

  real(R64)    :: time, startTime, endTime, outputStep, timeStep, charge, energy
  integer(I32) :: nTimeSteps, innerStep
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
  call TwoRdm_Fabricate
  call TwoRdmInit_Fabricate
  call Orbs_Fabricate
  call OrbsInit_Fabricate
  IntegratorListInput(1) % TimeDerivative => Method_Mb_OrbBased_TimeDerivativeOrbsLin
  IntegratorListInput(2) % TimeDerivative => Method_Mb_OrbBased_TimeDerivativeCoeffsPlusOrbsNonLin
  call IntegratorList_Fabricate(IntegratorListInput)
  call Propagator_Fabricate()
  call Purification_Fabricate
  call Reconstruction_Fabricate
  call Absorber_Fabricate

  call Say_Fabricate("program")
  startTime = Json_Get("program.startTime", 0.0_R64)
  endTime = Json_Get("program.endTime", 1.0_R64)
  outputStep = Json_Get("program.outputStep", 0.2_R64)
  nTimeSteps = Json_Get("program.nTimeSteps", 20)
  charge = Json_Get("sysPotential.ylm.coulomb.charge", 4.0_R64)

  call Grid_Setup
  call SysKinetic_Setup
  call SysPotential_Setup
  call SysGauge_Setup
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

  ! The one-body operator is self-adjoint only in the FEDVR-weighted metric
  integratorList(1) % e % metricWeights = MetricWeights(TwoRdm_nElements)

  !==========================================
  call Say_Section("propagate through laser pulse")
  !==========================================

  time = startTime
  timeStep = outputStep / nTimeSteps

  ! also fills Method_Mb_OrbBased_h2, only used by energy-conserving purification
  energy = Method_GetEnergy(time)
  call WriteObservables(time, energy, .true.)

  do while (time < endTime - 0.5_R64 * timeStep)

    do innerStep = 1, nTimeSteps
      call Propagator_Propagate(Method_state, time, time + timeStep)
      time = time + timeStep
      call Absorber_ApplyAbsorber(Orbs_orbs)
      call Purification_PurifyTwoRdmState(TwoRdm_twoRdmState, Method_Mb_OrbBased_h2)
    end do

    energy = Method_GetEnergy(time)
    call WriteObservables(time, energy, .false.)

  end do

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  subroutine WriteObservables(t, e, newFilesQ)
    real(R64), intent(in) :: t
    real(R64), intent(in) :: e
    logical, intent(in) :: newFilesQ

    complex(R64), allocatable :: rdm1(:, :)

    call TwoRdm_FillRdm1FromTwoRdmState(rdm1, TwoRdm_twoRdmState)
    call WriteLaserObservables(rdm1, t, charge, e, newFilesQ)

  end subroutine

  include '../common/BeCommon.inc'
  include '../common/BeObservables.inc'

end program
