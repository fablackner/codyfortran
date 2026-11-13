!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Anchored TD-2RDM reconstruction (reconstruction.singlet.anchor) with moving
!> orbitals: field-free TD-2RDM (NY closure with ccQ, exact field-free part,
!> split step) from the Be MCSCF ground state of T_Be3d_01 (its orbitals,
!> 2-RDM and up-up-down 3-RDM, read from this directory).
!>
!> With the anchor the reconstructed 3-RDM equals the exact one at the
!> reference, so the MCSCF state is a stationary point of the TD-2RDM
!> equations, also in the split-step gauge, where the orbitals rotate with
!> exp(-i h t) and the anchor is transported with their overlap with the
!> initial orbitals. The run starts from the same state in rotated orbitals
!> (1s and 2s mixed with complex coefficients), so that the overlap is neither
!> real nor symmetric and a transposed or conjugated transport shows. Checks
!> after program.endTime = 1: the natural occupations change by less than
!> 5e-6 (anchored: 1.4e-6 at the time of writing, the split-step error of
!> MCTDHF from the same state; unanchored NY: 1.8e-4) and the energy by less
!> than 5e-5 (1.8e-5, the split-step error). The energy check cannot see the
!> anchor.
program T_Be3d_03_TimePropagationTd2rdmAnchor
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_RdmDiagonalize
  use M_Utils_RdmBasisChange
  use M_Grid
  use M_SysKinetic
  use M_SysPotential
  use M_SysInteraction
  use M_Orbs
  use M_OrbsInit
  use M_TwoRdm
  use M_TwoRdm_Fermionic_Singlet
  use M_TwoRdmInit
  use M_Reconstruction
  use M_Method
  use M_Method_Mb_OrbBased
  use M_IntegratorList
  use M_Propagator
  use testdrive, only: check, error_type

  implicit none

  real(R64)    :: time, endTime, timeStep, energy0, energy
  real(R64), allocatable :: natOcc0(:), natOcc(:)
  complex(R64), allocatable :: u(:, :), rdm2ab(:, :, :, :)
  integer(I32) :: iOrb
  type(error_type), allocatable :: error
  type(T_IntegratorList_FabricateInput) :: IntegratorListInput(2)

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call Json_LoadJsonFile(manualFileName_="test/simulationTests/Be3d/T_Be3d_03_TimePropagationTd2rdmAnchor.json", &
                         relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call SysKinetic_Fabricate
  call SysPotential_Fabricate
  call SysInteraction_Fabricate
  call Method_Fabricate
  call TwoRdm_Fabricate
  call TwoRdmInit_Fabricate
  call Orbs_Fabricate
  call OrbsInit_Fabricate
  call Reconstruction_Fabricate
  IntegratorListInput(1) % TimeDerivative => Method_Mb_OrbBased_TimeDerivativeOrbsLin
  IntegratorListInput(2) % TimeDerivative => Method_Mb_OrbBased_TimeDerivativeCoeffsPlusOrbsNonLin
  call IntegratorList_Fabricate(IntegratorListInput)
  call Propagator_Fabricate()

  call Say_Fabricate("program")
  endTime = Json_Get("program.endTime", 1.0_R64)
  timeStep = Json_Get("program.timeStep", 0.02_R64)

  call Grid_Setup
  call SysKinetic_Setup
  call SysPotential_Setup
  call SysInteraction_Setup
  call TwoRdm_Setup
  call TwoRdmInit_Setup
  call Orbs_Setup
  call OrbsInit_Setup
  call Method_Setup
  call Reconstruction_Setup
  call IntegratorList_Setup
  call Propagator_Setup

  ! the same state in the orbitals phi U, U unitary and complex in the 1s-2s block
  allocate (u(Orbs_nOrbsInState, Orbs_nOrbsInState))
  u = 0.0_R64
  do iOrb = 1, Orbs_nOrbsInState
    u(iOrb, iOrb) = 1.0_R64
  end do
  u(1:2, 1:2) = reshape([(0.8_R64, 0.0_R64), (0.0_R64, 0.6_R64), (-0.6_R64, 0.0_R64), (0.0_R64, 0.8_R64)], [2, 2])
  Orbs_orbs(:, :) = matmul(Orbs_orbs, u)
  call TwoRdm_Fermionic_Singlet_FillRdm2abFromTwoRdmState(rdm2ab, TwoRdm_twoRdmState)
  call RdmBasisChange_Transform(rdm2ab, 4, conjg(transpose(u)))
  call TwoRdm_Fermionic_Singlet_FillTwoRdmStateFromRdm2ab(TwoRdm_twoRdmState, rdm2ab)

  !==========================================
  call Say_Section("perform action")
  !==========================================

  energy0 = Method_GetEnergy(0.0_R64)
  call FillNatOcc(natOcc0)

  time = 0.0_R64
  do while (time < endTime - 0.5_R64 * timeStep)
    call Propagator_Propagate(Method_state, time, time + timeStep)
    time = time + timeStep
  end do

  energy = Method_GetEnergy(time)
  call FillNatOcc(natOcc)

  !==========================================
  call Say_Section("printout results")
  !==========================================

  write (*, '(A,F22.14,A,ES10.2)') " energy: ", energy, "  change: ", energy - energy0
  write (*, '(A,*(ES14.6))') " natural occupations (spin up): ", natOcc
  write (*, '(A,ES10.2)') " max change of the natural occupations: ", maxval(abs(natOcc - natOcc0))

  call check(error, abs(energy - energy0) < 5e-5_R64)
  if (allocated(error)) error stop "T_Be3d_03_TimePropagationTd2rdmAnchor failure: energy"

  call check(error, maxval(abs(natOcc - natOcc0)) < 5e-6_R64)
  if (allocated(error)) error stop "T_Be3d_03_TimePropagationTd2rdmAnchor failure: natural occupations"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  subroutine FillNatOcc(natOcc)
    real(R64), intent(out), allocatable :: natOcc(:)

    complex(R64), allocatable :: rdm1(:, :), natOrbs(:, :)

    call TwoRdm_FillRdm1FromTwoRdmState(rdm1, TwoRdm_twoRdmState)
    call RdmDiagonalize_Rdm1(natOcc, natOrbs, rdm1(1:Orbs_nOrbsInState, 1:Orbs_nOrbsInState))

  end subroutine

end program
