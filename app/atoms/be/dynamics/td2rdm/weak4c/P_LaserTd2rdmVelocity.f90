!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> TD-2RDM for HHG in Be: real-time propagation of the (up,down) block of
!> the two-particle reduced density matrix together with time-dependent
!> orbitals (MCTDHF orbital equations) in a z-polarized velocity-gauge
!> few-cycle pulse.
!>
!> The 3-RDM in the collision operator is reconstructed with the contraction-
!> consistent Nakatsuji-Yasuda functional. Optionally (purification block in
!> the input) the D- and Q-conditions are enforced after every time step:
!> the contraction-free kernel of the negative parts of D and Q is
!> subtracted, which keeps norm and 1-RDM but not the energy (with
!> conserveEnergyQ also the energy; h2 is then rebuilt from the current
!> orbitals before every purification). The initial state
!> (orbXX_YY.in, rdm2BtXX_YY.in, read from orbsInit.load.dir and
!> twoRdmInit.load.dir) is the time-averaged MCTDHF ground state
!> (groundstate/td2rdm/P_GroundTd2rdm).
!> Grid, pulse, absorber and propagation scheme are identical to
!> P_LaserMctdhfVelocity (exactFieldFreeQ is unavailable with td2rdm.nWeakOrbs > 0).
!>
!> Output: observables.dat and natocc.dat as in P_LaserMctdhfVelocity, and at the end
!> the HHG spectrum hhg.dat (harmonic order, acceleration form, length form;
!> M_Utils_HhgSpectrum, up to program.hhgMaxOrder, default 300).
program P_LaserTd2rdmVelocity
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_HhgSpectrum
  use M_Grid
  use M_SysKinetic
  use M_SysPotential
  use M_SysGauge
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

  ! also fills Method_Mb_OrbBased_h2, only used by purification if Purification_usesH2Q
  energy = Method_GetEnergy(time)
  call WriteLaserObservables(time, charge, energy, .true.)

  do while (time < endTime - 0.5_R64 * timeStep)

    do innerStep = 1, nTimeSteps
      call Propagator_Propagate(Method_state, time, time + timeStep)
      time = time + timeStep
      call Absorber_ApplyAbsorber(Orbs_orbs)
      if (Absorber_keepOrthonormalQ) call Method_Mb_OrbBased_OrthonormalizeKeepingState(Method_state)
      ! h2 of the current orbitals (Method_GetEnergy fills it only at output steps)
      if (Purification_usesH2Q) call Method_Mb_OrbBased_FillH2(Method_Mb_OrbBased_h2, Orbs_orbs, time)
      call Purification_PurifyTwoRdmState(TwoRdm_twoRdmState, Method_Mb_OrbBased_h2)
      ! weak natural orbitals of the absorbed/purified state (nWeakOrbs > 0 only)
      call Method_Mb_OrbBased_Td2rdm_RefreshWeakOrbs
    end do

    energy = Method_GetEnergy(time)
    call WriteLaserObservables(time, charge, energy, .false.)

  end do

  !==========================================
  call Say_Section("hhg spectrum")
  !==========================================

  call HhgSpectrum_WriteFromTable("hhg.dat", "observables.dat", Json_Get("sysGauge.ylm.velocity.omega", 0.057_R64), &
                                  1, 4, colAcceleration_=5, maxOrder_=Json_Get("program.hhgMaxOrder", 300.0_R64))

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> FEDVR metric weights of a packed state [many-body part (nHead entries,
  !> unit weights) | orbitals (Grid_Ylm_weights each)] for the SIL integrator.
  function MetricWeights(nHead) result(res)
    use M_Grid
    use M_Grid_Ylm
    use M_Orbs

    integer(I32), intent(in) :: nHead
    real(R64), allocatable :: res(:)

    integer(I32) :: iOrb

    allocate (res(nHead + Grid_nPoints * Orbs_nOrbsInState))
    res(1:nHead) = 1.0_R64
    do iOrb = 1, Orbs_nOrbsInState
      res(nHead + (iOrb - 1) * Grid_nPoints + 1:nHead + iOrb * Grid_nPoints) = Grid_Ylm_weights(:)
    end do

  end function

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> The Hermitian matrix U diag(f) U^H of eigenvectors U and values f.
  function HermitianFunction(evecs, f) result(res)
    complex(R64), intent(in), contiguous :: evecs(:, :)
    real(R64), intent(in), contiguous    :: f(:)
    complex(R64), allocatable :: res(:, :)

    res = matmul(evecs, spread(f, 2, size(f)) * conjg(transpose(evecs)))

  end function

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Append one line to observables.dat and natocc.dat (replaced if newFilesQ):
  !>   observables.dat: time, E(t), A(t), <z>, -Z<z/r^3>, electrons in the box,
  !>                    energy
  !>   natocc.dat:      time, natural occupations of the spin-up 1-RDM
  !> Absorbed orbitals are no longer orthonormal: the electron number is
  !> Tr(rho S) and the occupations are the eigenvalues of S^1/2 rho S^1/2 with
  !> the orbital overlap matrix S (i.e. of the density matrix inside the box).
  !> The 1-RDM is taken from Method_state via the method-agnostic
  !> Method_Mb_OrbBased_FillRdm1 (MCTDHF and TD-2RDM alike).
  subroutine WriteLaserObservables(time, charge, energy, newFilesQ)
    use M_Utils_LapackLib
    use M_Utils_PrinterObservableYlm
    use M_Grid
    use M_Orbs
    use M_SysGauge_Ylm_Velocity
    use M_Method
    use M_Method_Mb_OrbBased

    real(R64), intent(in) :: time
    real(R64), intent(in) :: charge
    real(R64), intent(in) :: energy
    logical, intent(in) :: newFilesQ

    real(R64) :: dipole, acceleration, nElectrons
    real(R64), allocatable :: evals(:), occ(:)
    complex(R64), allocatable :: rdm1(:, :), s(:, :), evecs(:, :), sSqrt(:, :), rhoUp(:, :)
    integer(I32) :: nOS, i, j, io
    character(len=7) :: openStatus

    nOS = Orbs_nOrbsInState

    call Method_Mb_OrbBased_FillRdm1(rdm1, Method_state)

    call PrinterObservableYlm_DumpDipole(rdm1, Orbs_orbs, "", .false., time, dipole_=dipole)
    call PrinterObservableYlm_DumpDipoleAcceleration(rdm1, Orbs_orbs, charge, "", .false., time, &
                                                     acceleration_=acceleration)

    ! Overlap of the stored orbitals and spin-up block of the 1-RDM
    allocate (s(nOS, nOS))
    do j = 1, nOS
      do i = 1, nOS
        s(i, j) = Grid_InnerProduct(Orbs_orbs(:, i), Orbs_orbs(:, j))
      end do
    end do
    rhoUp = rdm1(1:nOS, 1:nOS)

    ! rdm1(i,j) = <a+_j a_i>: N = sum_ij rdm1(i,j) <phi_j|phi_i> per stored set
    nElectrons = real(sum(rhoUp * transpose(s)), kind=R64) * (size(rdm1, 1) / nOS)

    call LapackLib_DiagonalizeGeneric(evals, evecs, s, .true.)
    sSqrt = HermitianFunction(evecs, sqrt(max(evals, 0.0_R64)))
    call LapackLib_DiagonalizeGeneric(occ, evecs, matmul(conjg(transpose(sSqrt)), matmul(rhoUp, sSqrt)), .false.)

    openStatus = "unknown"
    if (newFilesQ) openStatus = "replace"

    open (newunit=io, file="observables.dat", status=trim(openStatus), position="append")
    if (newFilesQ) write (io, '(A)') "# time  E(t)  A(t)  <z>  -Z<z/r^3>  N_box  energy"
    write (io, '(7E24.14E3)') time, SysGauge_Ylm_Velocity_FieldAmplitude(time), &
      SysGauge_Ylm_Velocity_VectorPotentialAmplitude(time), dipole, acceleration, nElectrons, energy
    close (io)

    open (newunit=io, file="natocc.dat", status=trim(openStatus), position="append")
    write (io, '(*(E16.6E3))') time, occ(size(occ):1:-1)
    close (io)

    write (*, '(A,F9.3,A,ES12.4,A,ES14.6,A,ES14.6,A,F9.6,A,F18.10)') " t =", time, "  E =", &
      SysGauge_Ylm_Velocity_FieldAmplitude(time), "  <z> =", dipole, "  a =", acceleration, &
      "  N =", nElectrons, "  energy =", energy

  end subroutine

end program
