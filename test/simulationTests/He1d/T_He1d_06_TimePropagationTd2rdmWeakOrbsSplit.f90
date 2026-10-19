!> TD-2RDM with the active/weak natural-orbital split (nWeakOrbs = 2 of 4),
!> unrestricted orbitals, splitStep order2 (SIL linear, RK4 non-linear), same
!> kicked initial state as T_He1d_06_TimePropagationTd2rdmWeakOrbs. Checks the
!> O(dt^2) energy drift and the final natural occupations against the RK4
!> reference of the full time derivative. N = 2: no reconstruction.
program T_He1d_06_TimePropagationTd2rdmWeakOrbsSplit
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_Constants
  use M_Utils_RdmDiagonalize
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
  use M_Method_Mb_OrbBased_Td2rdm
  use M_IntegratorList
  use M_Propagator
  use M_Reconstruction
  use testdrive, only: check, error_type

  implicit none

  real(R64) :: t, startTime, endTime, outputStep, timeStep
  real(R64) :: energy, energy0, maxDrift, kick
  real(R64), allocatable :: natOcc(:)
  complex(R64), allocatable :: rdm1(:, :), natOrbs(:, :)
  integer(I32) :: nTimeSteps, innerStep, iOrb
  type(error_type), allocatable :: error
  character(len=:), allocatable :: jsonFileName
  type(T_IntegratorList_FabricateInput) :: IntegratorListInput(2)

  ! final natural occupations (body type 1) of the RK4 propagation of the full
  ! time derivative (T_He1d_06_TimePropagationTd2rdmWeakOrbs)
  real(R64), parameter :: natOccRef(4) = [ &
                          9.9225227458314180e-01_R64, &
                          7.5190775378442048e-03_R64, &
                          2.2854423643292498e-04_R64, &
                          1.0364258414043142e-07_R64]

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  jsonFileName = "test/simulationTests/He1d/T_He1d_06_TimePropagationTd2rdmWeakOrbsSplit.json"
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
  IntegratorListInput(1) % TimeDerivative => Method_Mb_OrbBased_TimeDerivativeOrbsLin
  IntegratorListInput(2) % TimeDerivative => Method_Mb_OrbBased_TimeDerivativeCoeffsPlusOrbsNonLin
  call IntegratorList_Fabricate(IntegratorListInput)
  call Propagator_Fabricate()
  Reconstruction_AddCollisionOp => AddNoCollisionOp

  call Say_Fabricate("program")
  startTime = Json_Get("program.startTime", 0.0_R64)
  endTime = Json_Get("program.endTime", 1.0_R64)
  outputStep = Json_Get("program.outputStep", 0.1_R64)
  nTimeSteps = Json_Get("program.nTimeSteps", 10)
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

  !==========================================
  call Say_Section("perform action")
  !==========================================

  if (Method_Mb_OrbBased_Td2rdm_nWeakOrbs .eq. 0) error stop "T_He1d_06_TimePropagationTd2rdmWeakOrbsSplit: split not active"

  do iOrb = 1, Orbs_nOrbsInState
    Orbs_orbs(:, iOrb) = Orbs_orbs(:, iOrb) * exp(IU * kick * Grid_Linear_xCoord(:))
  end do
  ! The kick changes the weak natural orbitals on the grid that Method_Setup
  ! stored for the linear step (without this the scheme is first order).
  call Method_Mb_OrbBased_Td2rdm_RefreshWeakOrbs

  energy0 = Method_GetEnergy(startTime)
  maxDrift = 0.0_R64

  t = startTime
  timeStep = outputStep / nTimeSteps

  do while (t < endTime - 0.5_R64 * timeStep)
    do innerStep = 1, nTimeSteps
      call Propagator_Propagate(Method_state, t, t + timeStep)
      t = t + timeStep
    end do

    energy = Method_GetEnergy(t)
    maxDrift = max(maxDrift, abs(energy - energy0))
    write (*, '(A,F8.3,A,F20.14,A,ES10.2)') ' time = ', t, '    energy = ', energy, '    drift = ', energy - energy0
  end do

  !==========================================
  call Say_Section("printout results")
  !==========================================

  call TwoRdm_FillRdm1FromTwoRdmState(rdm1, TwoRdm_twoRdmState)
  call RdmDiagonalize_Rdm1(natOcc, natOrbs, rdm1(1:4, 1:4))

  print *
  print *, "max energy drift: ", maxDrift
  write (*, '(A,*(ES24.16))') " final natural occupations (body type 1): ", natOcc
  write (*, '(A,*(ES10.2))') " deviation from the RK4 reference:         ", natOcc - natOccRef

  ! measured: drift 3.7e-8 (dt = 0.002, second order in dt); natural
  ! occupations 1.6e-8 off the reference, which is the O(eps/n^2) part of
  ! D1 M that the split step drops (4.8e-10 with regularizationParameter
  ! 1e-16). A stale weak space (no refresh after the kick) gives 5e-6.
  call check(error, maxDrift, 0.0_R64, thr=1e-7_R64)
  if (allocated(error)) error stop "T_He1d_06_TimePropagationTd2rdmWeakOrbsSplit failure: energy not conserved"

  call check(error, maxval(abs(natOcc - natOccRef)), 0.0_R64, thr=1e-7_R64)
  if (allocated(error)) error stop "T_He1d_06_TimePropagationTd2rdmWeakOrbsSplit failure: natural occupations off reference"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  !> N = 2: the 3-RDM and with it the collision operator vanish.
  subroutine AddNoCollisionOp(dTwoRdmState, twoRdmState, h2)
    complex(R64), intent(inout), contiguous :: dTwoRdmState(:)
    complex(R64), intent(in), contiguous    :: twoRdmState(:)
    complex(R64), intent(in), contiguous    :: h2(:, :, :, :)

    if (.false.) dTwoRdmState(1) = dTwoRdmState(1) + twoRdmState(1) * h2(1, 1, 1, 1)

  end subroutine

end program
