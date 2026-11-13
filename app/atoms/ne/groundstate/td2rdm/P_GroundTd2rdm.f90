!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> TD-2RDM ground state by adiabatic switching, with moving orbitals and no
!> CI anywhere: the Hartree-Fock determinant (twoRdmInit.excited) in the
!> Hartree-Fock orbitals plus correlating orbitals, the interaction switched
!> on over -T < t < 0 (method.mb.orbBased.adiabaticSwitching, the frozen
!> Hartree-Fock mean field as the reference), then program.relaxTime of
!> propagation with the physical Hamiltonian, which shows how stationary the
!> state is under the TD-2RDM equations.
!>
!> Input:  orbXX_YY.in of the Hartree-Fock orbitals (orbsInit.load, the
!>         correlating orbitals from the backend next to it, nOrbsIn)
!> Output: switching.dat (t, lambda, energy <H(t)>, min eigenvalues of D and
!>         Q of the up-down block, natural occupations), orbXX_YY.in,
!>         rdm2BtXX_YY.in of the state at t = 0
program P_GroundTd2rdm
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_RdmDiagonalize
  use M_Utils_RdmHole
  use M_Utils_RdmTrace
  use M_Grid
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

  implicit none

  real(R64)    :: time, relaxTime, outputStep, timeStep, energy0
  integer(I32) :: nTimeSteps, innerStep, io
  integer(I64) :: c0, c1, cRate
  complex(R64), allocatable :: rdm2(:, :, :, :)
  type(T_IntegratorList_FabricateInput) :: IntegratorListInput(2)

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call Json_LoadJsonFile()

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
  call Reconstruction_Fabricate

  call Say_Fabricate("program")
  relaxTime = Json_Get("program.relaxTime", 10.0_R64)
  outputStep = Json_Get("program.outputStep", 0.2_R64)
  nTimeSteps = Json_Get("program.nTimeSteps", 20)
  if (.not. (Method_Mb_OrbBased_switchTime > 0.0_R64)) &
    error stop "P_GroundTd2rdm needs method.mb.orbBased.adiabaticSwitching"

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
  call Reconstruction_Setup

  ! the correlating orbitals of orbsInit are not orthogonal to the occupied
  ! ones (Gram-Schmidt keeps the occupied orbitals, and with them the
  ! reference)
  call Orbs_Orthonormalize(Orbs_orbs)
  call Method_Mb_OrbBased_Td2rdm_RefreshWeakOrbs
  call Method_Mb_OrbBased_SetSwitchingReference

  ! the one-body operator is self-adjoint only in the FEDVR-weighted metric
  integratorList(1) % e % metricWeights = MetricWeights(TwoRdm_nElements)

  !==========================================
  call Say_Section("adiabatic switching and relaxation")
  !==========================================

  open (newunit=io, file="switching.dat", status="replace")
  timeStep = outputStep / nTimeSteps
  time = -Method_Mb_OrbBased_switchTime
  call WriteProgress(time)
  call system_clock(c0, cRate)

  do while (time < relaxTime - 0.5_R64 * timeStep)
    do innerStep = 1, nTimeSteps
      call Propagator_Propagate(Method_state, time, time + timeStep)
      time = time + timeStep
    end do
    call WriteProgress(time)

    if (abs(time) < 0.5_R64 * timeStep) then
      energy0 = Method_GetEnergy(0.0_R64)
      call TwoRdm_FillRdm2FromTwoRdmState(rdm2, TwoRdm_twoRdmState)
      call Orbs_SaveOrbs(Orbs_orbs)
      call Method_Mb_OrbBased_SaveRdm2Blocks(rdm2)
      write (*, '(A,F22.12)') " ground state (t = 0, saved): energy ", energy0
    end if
  end do
  close (io)

  call system_clock(c1)
  write (*, '(A,F22.12,A,ES10.2)') " energy after the relaxation: ", Method_GetEnergy(time), &
    ",  drift ", Method_GetEnergy(time) - energy0
  write (*, '(A,F10.2,A)') " wall time ", real(c1 - c0, R64) / cRate, " s"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  subroutine WriteProgress(t)
    real(R64), intent(in) :: t

    real(R64), allocatable :: natOcc(:)
    complex(R64), allocatable :: rdm1(:, :), natOrbs(:, :), rdm2(:, :, :, :), rdm2ab(:, :, :, :)
    real(R64) :: energy, minD, minQ
    integer(I32) :: n

    n = Orbs_nOrbsInState
    call TwoRdm_FillRdm1FromTwoRdmState(rdm1, TwoRdm_twoRdmState)
    call TwoRdm_FillRdm2FromTwoRdmState(rdm2, TwoRdm_twoRdmState)
    call RdmDiagonalize_Rdm1(natOcc, natOrbs, rdm1(1:n, 1:n))
    call RdmTrace_FillRdm2abFromRdm2(rdm2ab, rdm2, n, n)
    call RdmHole_MinEigenvaluesDQab(minD, minQ, rdm2ab)
    energy = Method_GetEnergy(t)

    write (io, '(F10.4, F8.4, F24.14, 2ES13.4, *(ES15.7))') t, Method_Mb_OrbBased_SwitchFactor(t), energy, &
      minD, minQ, natOcc
    flush (io)
    write (*, '(A,F9.3,A,F7.4,A,F22.12,A,2ES11.3)') " t =", t, "  lambda =", Method_Mb_OrbBased_SwitchFactor(t), &
      "  energy =", energy, "  min D, Q:", minD, minQ

  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> FEDVR metric weights of a packed state [many-body part (nHead entries,
  !> unit weights) | orbitals (Grid_Ylm_weights each)] for the SIL integrator.
  function MetricWeights(nHead) result(res)
    use M_Grid_Ylm

    integer(I32), intent(in) :: nHead
    real(R64), allocatable :: res(:)

    integer(I32) :: iOrb

    allocate (res(nHead + Grid_nPoints * Orbs_nOrbsInState))
    res(1:nHead) = 1.0_R64
    do iOrb = 1, Orbs_nOrbsInState
      res(nHead + (iOrb - 1) * Grid_nPoints + 1:nHead + iOrb * Grid_nPoints) = Grid_Ylm_weights(:)
    end do

  end function

end program
