!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Frozen mean-field reference (method.mb.orbBased.exactFieldFreeMeanFieldQ) in
!> MCTDHF with more orbitals than electrons, where its CI term matters: H2
!> (R = 1.4) on a small finite-element mesh (N = 15341), restricted MCTDHX
!> with nOrbs = [2, 2] (sigma_g and sigma_u from STOs, CI coefficient 1 on
!> |sigma_g sigma_g>), every orbital kicked by exp(i kappa z) with kappa =
!> program.kappa = 0.1. The split with the reference (exact part h0 + U;
!> nonlinear part: orbitals -i Q (W - U) phi, CI -i (H2 - U^) c) is a gauge
!> transformation of the plain split (exact part h0), so both must follow one
!> trajectory. With nOrbs = nBodies (TDHF, T_Fem_09) the CI projector removes
!> U^ c; here it does not, and without the CI term the orbitals would rotate
!> by <phi|U|phi> while the coefficients do not follow.
!>
!> Checks at dt = program.timeStep = 0.005 up to program.endTime = 0.25
!> (measured values in parentheses; a scratch build without the CI term in
!> brackets, which fails checks 1 and 2):
!> 1. the induced dipoles d(t) = <z>(t) - <z>(0) of both splits agree at every
!>    step to program.dipoleTolerance = 1e-5 relative to max |d| (1.6e-7
!>    [4.5e-3]);
!> 2. the energies of both splits agree at every step to
!>    program.energyAgreement = 3e-6 Ha (3.4e-7 [1.7e-3]), and the energy
!>    with the reference stays within program.energyTolerance = 1e-5 Ha of
!>    E(0+) (9.4e-7: both splits lose about 1e-6 Ha within the first 0.05 a.u.
!>    at any dt, presumably from the regularized inverse 1-RDM while sigma_u
!>    is nearly empty);
!> 3. N = 2 to 1e-8 (2.8e-10) and <x> = <y> = 0 to 1e-9 (5e-16) in both runs;
!> 4. the CI term acts: at the end its part orthogonal to c (what the CI
!>    projector keeps), |(1 - |c><c|) U^ c|, exceeds MIN_CI_TERM = 1e-2 Ha
!>    (4.5e-2; zero at t = 0, where c is one configuration and
!>    <sigma_u|U|sigma_g> = 0 by symmetry; the smaller natural occupation
!>    reaches 3.1e-3).
program T_Fem_10_H2MctdhfMeanFieldReference
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Grid
  use M_Orbs
  use M_Coeffs
  use M_Method
  use M_Method_Mb_OrbBased
  use testdrive, only: check, error_type

  implicit none

  !> Number of electrons
  real(R64), parameter :: N_ELECTRONS = 2.0_R64
  !> Lower bound of |(1 - |c><c|) U^ c| at the end of the run with the
  !> reference (Ha)
  real(R64), parameter :: MIN_CI_TERM = 1e-2_R64

  type(error_type), allocatable :: error
  complex(R64), allocatable :: state0(:), uOrbs(:, :), uMat(:, :), rdm1(:, :), rdmSpatial(:, :), ciTermVec(:)
  real(R64), allocatable :: dipolePlain(:), dipoleRef(:), energyPlain(:), energyRef(:)
  real(R64) :: kappa, dt, endTime, dipoleTolerance, energyTolerance, energyAgreement
  real(R64) :: maxNormDevPlain, maxNormDevRef, maxTransversePlain, maxTransverseRef, secondsPlain, secondsRef
  real(R64) :: maxResponse, maxDipoleDiff, maxEnergyDiff, maxEnergyDevRef, coupling, ciTerm, occupationSmall
  integer(I32) :: nSteps, iStep

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call SetupFemDynamics("test/simulationTests/Fem/T_Fem_10_H2MctdhfMeanFieldReference.json")

  call Say_Fabricate("program")
  kappa = Json_Get("program.kappa", 0.1_R64)
  dt = Json_Get("program.timeStep", 0.005_R64)
  endTime = Json_Get("program.endTime", 0.25_R64)
  dipoleTolerance = Json_Get("program.dipoleTolerance", 1e-5_R64)
  energyTolerance = Json_Get("program.energyTolerance", 1e-5_R64)
  energyAgreement = Json_Get("program.energyAgreement", 3e-6_R64)
  nSteps = nint(endTime / dt)

  if (.not. Method_Mb_OrbBased_exactFieldFreeMeanFieldQ) &
    error stop "T_Fem_10_H2MctdhfMeanFieldReference: the JSON must set method.mb.orbBased.exactFieldFreeMeanFieldQ"
  if (Orbs_nOrbsInState .ne. 2) &
    error stop "T_Fem_10_H2MctdhfMeanFieldReference: the JSON must declare two stored orbitals (nOrbs = [2, 2], restricted)"

  !==========================================
  call Say_Section("perform action")
  !==========================================

  state0 = Method_state

  call RunSplit(dipolePlain, energyPlain, maxNormDevPlain, maxTransversePlain, secondsPlain, .false.)
  call RunSplit(dipoleRef, energyRef, maxNormDevRef, maxTransverseRef, secondsRef, .true.)

  ! The CI term of the reference (built after the kick) at the end: U^ c
  ! without its component along c (the CI coefficients are the first
  ! Coeffs_nCoeffs entries of the state), |<phi_1|U|phi_2>|, and the smaller
  ! natural occupation (spatial 1-RDM, 2 x 2 Hermitian)
  allocate (uOrbs, mold=Orbs_orbs)
  uOrbs(:, :) = (0.0_R64, 0.0_R64)
  call Method_Mb_OrbBased_ApplyFieldFreeMeanField(uOrbs, Orbs_orbs, h1_=uMat)
  coupling = abs(uMat(1, 2))
  allocate (ciTermVec(Coeffs_nCoeffs))
  ciTermVec(:) = (0.0_R64, 0.0_R64)
  call Coeffs_ApplyH1FillRdm1(Method_state(1:Coeffs_nCoeffs), dCoeffs_=ciTermVec, h1_=uMat)
  call Coeffs_ProjectOnSubspace(ciTermVec, Method_state(1:Coeffs_nCoeffs))
  ciTerm = sqrt(sum(abs(ciTermVec)**2))
  call Method_Mb_OrbBased_FillRdm1(rdm1, Method_state)
  call Method_Mb_OrbBased_FillSpatialRdm1(rdmSpatial, rdm1)
  occupationSmall = 0.5_R64 * (real(rdmSpatial(1, 1) + rdmSpatial(2, 2), R64) &
                               - sqrt(real(rdmSpatial(1, 1) - rdmSpatial(2, 2), R64)**2 + 4.0_R64 * abs(rdmSpatial(1, 2))**2))

  maxResponse = 0.0_R64
  maxDipoleDiff = 0.0_R64
  maxEnergyDiff = 0.0_R64
  maxEnergyDevRef = 0.0_R64
  do iStep = 1, nSteps
    maxResponse = max(maxResponse, abs(dipolePlain(iStep)))
    maxDipoleDiff = max(maxDipoleDiff, abs(dipoleRef(iStep) - dipolePlain(iStep)))
    maxEnergyDiff = max(maxEnergyDiff, abs(energyRef(iStep) - energyPlain(iStep)))
    maxEnergyDevRef = max(maxEnergyDevRef, abs(energyRef(iStep) - energyRef(0)))
  end do

  !==========================================
  call Say_Section("printout results")
  !==========================================

  print *
  print *, "kick kappa, dt, steps:                         ", kappa, dt, nSteps
  print *, "E(0+) plain, reference:                        ", energyPlain(0), energyRef(0)
  print *, "1. max |d_ref(t) - d_plain(t)|, max |d|:       ", maxDipoleDiff, maxResponse
  print *, "     relative:                                 ", maxDipoleDiff / maxResponse
  print *, "2. max |E_ref(t) - E_plain(t)|:                ", maxEnergyDiff
  print *, "   max |E_ref(t) - E_ref(0+)|:                 ", maxEnergyDevRef
  print *, "3. max |N(t) - 2| plain, reference:            ", maxNormDevPlain, maxNormDevRef
  print *, "   max |<x>|, |<y>| plain, reference:          ", maxTransversePlain, maxTransverseRef
  print *, "4. |(1 - |c><c|) U^ c| at the end:              ", ciTerm
  print *, "   |<phi_1|U|phi_2>| at the end:               ", coupling
  print *, "   smaller natural occupation at the end:      ", occupationSmall
  write (*, '(1X, A, F8.4, A, F8.4, A)') "time per step: ", secondsPlain / nSteps, " s (plain), ", &
    secondsRef / nSteps, " s (reference)"

  call check(error, maxDipoleDiff < dipoleTolerance * maxResponse)
  if (allocated(error)) error stop "T_Fem_10_H2MctdhfMeanFieldReference failure: the dipole with the reference &
                                  &leaves the trajectory of the plain split"

  call check(error, maxEnergyDiff < energyAgreement .and. maxEnergyDevRef < energyTolerance)
  if (allocated(error)) error stop "T_Fem_10_H2MctdhfMeanFieldReference failure: energy with the reference not &
                                  &conserved or off the plain split"

  call check(error, maxNormDevPlain < 1e-8_R64 .and. maxNormDevRef < 1e-8_R64 .and. maxTransversePlain < 1e-9_R64 &
             .and. maxTransverseRef < 1e-9_R64)
  if (allocated(error)) error stop "T_Fem_10_H2MctdhfMeanFieldReference failure: particle number or transverse &
                                  &dipole not conserved"

  call check(error, ciTerm > MIN_CI_TERM)
  if (allocated(error)) error stop "T_Fem_10_H2MctdhfMeanFieldReference failure: the CI term of the reference is &
                                  &too small, the run does not test it"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  !> The kick and the propagation of one split from the initial state state0:
  !> plain (meanFieldQ false: exact part h0) or with the frozen mean-field
  !> reference (built from the kicked state at the first exact-part call). A
  !> fresh Cayley setup releases a reference built before. Returns the induced
  !> dipole d(t) = <z>(t) - <z>(0) and E_el(t) at every step (0:nSteps), the
  !> largest deviations of N from 2 and of <x>, <y> from 0, and the seconds
  !> of the propagation.
  subroutine RunSplit(dipoleZ, energies, maxNormDev, maxTransverse, seconds, meanFieldQ)
    use M_Propagator

    real(R64), intent(out), allocatable :: dipoleZ(:)
    real(R64), intent(out), allocatable :: energies(:)
    real(R64), intent(out) :: maxNormDev
    real(R64), intent(out) :: maxTransverse
    real(R64), intent(out) :: seconds
    logical, intent(in) :: meanFieldQ

    real(R64) :: norm, dipole(3), dipoleStart(3), time
    integer(I64) :: count0
    integer(I32) :: iStepRun

    allocate (dipoleZ(0:nSteps), energies(0:nSteps))
    Method_state(:) = state0(:)
    Method_Mb_OrbBased_exactFieldFreeMeanFieldQ = meanFieldQ
    call Method_Mb_OrbBased_SetupFieldFree

    call FemObservables(energies(0), norm, dipoleStart, 0.0_R64, "")
    call KickOrbitalsZ(kappa)
    call FemObservables(energies(0), norm, dipole, 0.0_R64, "")
    dipoleZ(0) = dipole(3) - dipoleStart(3)
    maxNormDev = abs(norm - N_ELECTRONS)
    maxTransverse = maxval(abs(dipole(1:2)))

    time = 0.0_R64
    call system_clock(count0)
    do iStepRun = 1, nSteps
      call Propagator_Propagate(Method_state, time, time + dt)
      time = iStepRun * dt
      call FemObservables(energies(iStepRun), norm, dipole, time, "")
      dipoleZ(iStepRun) = dipole(3) - dipoleStart(3)
      maxNormDev = max(maxNormDev, abs(norm - N_ELECTRONS))
      maxTransverse = max(maxTransverse, maxval(abs(dipole(1:2))))
    end do
    seconds = ElapsedSeconds(count0)

    write (*, '(1X, A, L1, A, ES12.4, A, ES12.4, A, ES10.2)') "split with reference ", meanFieldQ, ": d(T) ", &
      dipoleZ(nSteps), ", E(T) - E(0+) ", energies(nSteps) - energies(0), ", max |N - 2| ", maxNormDev

  end subroutine

  !> The kick phi <- exp(i kick z) phi of every stored orbital at the nodes
  !> (a one-body unitary on the many-body state; Orbs_orbs is the orbital
  !> part of Method_state).
  !> Local copy: including FemKick.inc would leave RunKick unused (-Wunused-function -Werror).
  subroutine KickOrbitalsZ(kick)
    use M_Grid_Fem

    real(R64), intent(in) :: kick

    real(R64), allocatable :: zNodes(:)
    integer(I32) :: iOrb

    call Grid_Fem_FillNodeCoordinate(zNodes, [0.0_R64, 0.0_R64, 1.0_R64])
    do iOrb = 1, Orbs_nOrbsInState
      Orbs_orbs(:, iOrb) = exp(cmplx(0.0_R64, kick * zNodes(:), kind=R64)) * Orbs_orbs(:, iOrb)
    end do

  end subroutine

  include 'FemDynamics.inc'

end program
