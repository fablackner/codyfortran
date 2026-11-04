!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Frozen mean-field reference of the Cayley exact part
!> (method.mb.orbBased.exactFieldFreeMeanFieldQ) in TDHF of H2 (R = 1.4) on the
!> finite-element grid: the exact part propagates h0 + U with U = J0 - K0 of
!> the state at the first propagation (static, Hermitian, K0 by ACE), the
!> nonlinear part the interaction minus U (orbitals -i Q (W - U) phi, CI
!> (H2 - U^) c). Same system, orbital and wiring as T_Fem_05 .. T_Fem_07
!> (restricted MCTDHX with nOrbs = nBodies = [1, 1] from the RHF orbital of
!> T_Fem_03, splitStep order2 with exactPartScheme strang, RK4).
!>
!> Checks (measured values in parentheses):
!> 1. Reference at the ground state (t = 0, no kick): the nonlinear orbital
!>    derivative -i Q (W - U) phi vanishes, relative to -i Q W phi (the
!>    derivative without the reference) below 1e-9 (2.5e-11, the 1-RDM
!>    regularization n/(n^2 + 1e-10) at n = 2); the CI derivative is zero
!>    (TDHF: one configuration, U^ c parallel to c).
!> 2. Ground state at the large step program.groundTimeStep = 0.1 up to
!>    program.groundEndTime = 0.3: stationary, |E(t) - E(0)| < 1e-10 (6e-13),
!>    1 - |<phi(0)|phi(t)>| < 1e-10 (2e-13; the phase rotates), N = 2 to 1e-10
!>    (1e-12); at most MAX_ITERATIONS_GROUND GMRES iterations per Cayley solve
!>    (dt/4 = 0.025: 6).
!> 3. Kick kappa = program.kappa = 1e-3 along z at dt = program.timeStep =
!>    0.02 (the step of T_Fem_06_H2Kick) up to program.endTime = 0.6, with a
!>    fresh reference (built after the kick): the checks of T_Fem_06 (kick
!>    energy, TRK slope, energy conserved to 1e-8, N = 2, <x> = <y> = 0), and
!>    the induced dipole agrees with T_Fem_06's dipoleKick.dat (no reference,
!>    same dt; header checked as in T_Fem_07: kappa, dt, mesh, E_el(0),
!>    complete run) at every common step to program.dipoleTolerance = 1e-5
!>    relative to the largest response (2.7e-7: two splittings of the same
!>    equations, each off by about 6e-6 at dt = 0.02, see 4).
!> 4. Convergence of the splitting with the reference: the kick at dt/2 and
!>    dt/4 up to program.refineEndTime = 0.1; at that time the order
!>    p = log2(|d(dt) - d(dt/2)| / |d(dt/2) - d(dt/4)|) of the induced dipole
!>    d lies in [1.7, 2.3] (2.00), and the errors at dt of both splittings
!>    against the Richardson limit of the refined runs are printed (-5.7e-6
!>    with, -6.4e-6 without the reference, so both converge to the same
!>    trajectory); the error with the reference stays below
!>    program.refineTolerance = 2e-5 relative to the response.
!> 5. A reference built explicitly before the kick, from the ground state (as
!>    app/fem_molecules does), is kept at the first call: at the kicked state
!>    the nonlinear orbital derivative is then the response
!>    -i Q (W(0+) - W0) phi = O(kappa), above MIN_RESPONSE_DERIVATIVE = 1e-4
!>    relative to the ground state's derivative without reference (3.0e-3),
!>    while the default reference, built from the kicked state at that call,
!>    makes it vanish there (below 1e-9; 2.5e-11).
!> Writes dipoleKickMeanField.dat (dt run) and dipoleKickMeanFieldRefine*.dat.
program T_Fem_09_H2KickMeanFieldReference
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Grid
  use M_Orbs
  use M_Coeffs
  use M_Method
  use M_Method_Mb_OrbBased
  use M_Propagator
  use testdrive, only: check, error_type

  implicit none

  !> Number of electrons
  real(R64), parameter :: N_ELECTRONS = 2.0_R64
  !> Bound of the GMRES iterations per Cayley solve of the ground-state run
  !> (dt/4 = 0.025; 6 measured, + 2)
  integer(I32), parameter :: MAX_ITERATIONS_GROUND = 8
  !> Lower bound of the nonlinear orbital derivative at the kicked state with
  !> the reference from the ground state, relative to the ground state's
  !> derivative without reference (check 5)
  real(R64), parameter :: MIN_RESPONSE_DERIVATIVE = 1e-4_R64

  type(error_type), allocatable :: error
  complex(R64), allocatable, target :: groundState(:), dState(:)
  complex(R64), allocatable :: phi0(:)
  real(R64), allocatable :: times(:), dipoleZ(:), timesRef(:), dipoleZRef(:), refineTimeSteps(:)
  real(R64), allocatable :: refineTimes(:), refineDipoleZ(:)
  real(R64) :: kappa, dt, endTime, groundTimeStep, groundEndTime, refineEndTime, dipoleTolerance, refineTolerance
  real(R64) :: derivativeRatio, derivativeCoeffs, energyGround, energy, norm, dipole(3), time, overlap
  real(R64) :: derivativePlain, derivativeExplicit, derivativeDefault
  complex(R64), allocatable :: rdm1(:, :)
  real(R64) :: groundEnergyDev, groundNormDev, groundOverlapDev, energyKick, slopeExact, slopePropagated
  real(R64) :: maxEnergyDev, maxNormDev, maxTransverse, secondsPerStep, kickEnergyDev, trkDev, slopeDev
  real(R64) :: kappaRef, dtRef, endTimeRef, energyGroundRef, maxResponse, maxDiffRef, maxTimeDev
  real(R64) :: dipoleAt(3), order, richardson, errorMeanField, errorPlain, responseRef
  real(R64) :: refineEnergyGround, refineEnergyKick, refineSlopeExact, refineSlopePropagated, refineEnergyDev
  real(R64) :: refineNormDev, refineTransverse, refineSecondsPerStep, maxRefineEnergyDev
  integer(I32) :: nSteps, nStepsRef, nPointsRef, nGround, iStep, iRefine, groundMaxIterations, iC
  character(len=:), allocatable :: headerRef
  character(len=64) :: fileName

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call SetupFemDynamics("test/simulationTests/Fem/T_Fem_09_H2KickMeanFieldReference.json")

  call Say_Fabricate("program")
  kappa = Json_Get("program.kappa", 1e-3_R64)
  dt = Json_Get("program.timeStep", 0.02_R64)
  endTime = Json_Get("program.endTime", 0.6_R64)
  groundTimeStep = Json_Get("program.groundTimeStep", 0.1_R64)
  groundEndTime = Json_Get("program.groundEndTime", 0.3_R64)
  refineEndTime = Json_Get("program.refineEndTime", 0.1_R64)
  refineTimeSteps = Json_Get("program.refineTimeSteps", [0.01_R64, 0.005_R64])
  dipoleTolerance = Json_Get("program.dipoleTolerance", 1e-5_R64)
  refineTolerance = Json_Get("program.refineTolerance", 2e-5_R64)
  nSteps = nint(endTime / dt)
  nGround = nint(groundEndTime / groundTimeStep)

  if (.not. Method_Mb_OrbBased_exactFieldFreeMeanFieldQ) &
    error stop "T_Fem_09_H2KickMeanFieldReference: the JSON must set method.mb.orbBased.exactFieldFreeMeanFieldQ"
  if (size(refineTimeSteps) .ne. 2) error stop "T_Fem_09_H2KickMeanFieldReference: program.refineTimeSteps needs 2 steps"
  if (refineEndTime > endTime) error stop "T_Fem_09_H2KickMeanFieldReference: program.refineEndTime > program.endTime"

  !==========================================
  call Say_Section("perform action")
  !==========================================

  groundState = Method_state
  phi0 = Orbs_orbs(:, 1)
  allocate (dState, mold=Method_state)

  ! T_Fem_06's run (no reference, same kappa and dt): checked before the
  ! propagations
  call ReadKickFile(headerRef, kappaRef, dtRef, endTimeRef, energyGroundRef, nStepsRef, nPointsRef, timesRef, &
                    dipoleZRef, "dipoleKick.dat")
  write (*, '(1X, A, A)') "reference (T_Fem_06_H2Kick): ", headerRef
  if (.not. (abs(kappaRef - kappa) <= 1e-12_R64 * kappa .and. abs(dtRef - dt) <= 1e-12_R64 * dt &
             .and. nPointsRef .eq. Grid_nPoints .and. endTime <= endTimeRef * (1.0_R64 + 1e-12_R64))) &
    error stop "T_Fem_09_H2KickMeanFieldReference: dipoleKick.dat is from another kappa, step, end time or mesh &
               &(stale or foreign file; run T_Fem_06_H2Kick first)"
  if (size(timesRef) .ne. nStepsRef + 1) &
    error stop "T_Fem_09_H2KickMeanFieldReference: dipoleKick.dat is incomplete (T_Fem_06_H2Kick did not finish its run)"

  ! 1. The nonlinear derivative at the ground state without and with the
  !    reference (the second call builds it from the ground state at t = 0)
  Method_Mb_OrbBased_exactFieldFreeMeanFieldQ = .false.
  call Method_Mb_OrbBased_TimeDerivativeCoeffsPlusOrbsNonLin(dState, Method_state, 0.0_R64)
  derivativePlain = OrbitalNorm(dState)
  Method_Mb_OrbBased_exactFieldFreeMeanFieldQ = .true.
  call Method_Mb_OrbBased_TimeDerivativeCoeffsPlusOrbsNonLin(dState, Method_state, 0.0_R64)
  write (*, '(1X, A, ES12.4, A, ES12.4)') "nonlinear orbital derivative at the ground state: without reference ", &
    derivativePlain, ", with ", OrbitalNorm(dState)
  derivativeRatio = OrbitalNorm(dState) / derivativePlain
  derivativeCoeffs = 0.0_R64
  do iC = 1, Coeffs_nCoeffs
    derivativeCoeffs = max(derivativeCoeffs, abs(dState(iC)))
  end do

  ! 2. The ground state at a large step, with that reference
  energyGround = Method_GetEnergy(0.0_R64)
  groundEnergyDev = 0.0_R64
  groundNormDev = 0.0_R64
  groundOverlapDev = 0.0_R64
  time = 0.0_R64
  do iStep = 1, nGround
    call Propagator_Propagate(Method_state, time, time + groundTimeStep)
    time = iStep * groundTimeStep
    call FemObservables(energy, norm, dipole, time, "")
    overlap = abs(Grid_InnerProduct(phi0, Orbs_orbs(:, 1)))
    groundEnergyDev = max(groundEnergyDev, abs(energy - energyGround))
    groundNormDev = max(groundNormDev, abs(norm - N_ELECTRONS))
    groundOverlapDev = max(groundOverlapDev, 1.0_R64 - overlap)
    write (*, '(1X, A, F6.3, A, F18.13, A, ES10.2, A, ES10.2, A, ES10.2)') "ground state, t = ", time, ": E_el ", &
      energy, ", E - E(0) ", energy - energyGround, ", N - 2 ", norm - N_ELECTRONS, ", 1 - |<phi0|phi>| ", 1.0_R64 - overlap
  end do
  groundMaxIterations = Method_Mb_OrbBased_cayleyMaxIterations
  write (*, '(1X, A, F5.2, A, I0, A, I0, A, F6.2, A, I0)') "ground state at dt = ", groundTimeStep, ": ", &
    Method_Mb_OrbBased_cayleyNSolves, " Cayley solves, ", Method_Mb_OrbBased_cayleyNIterations, &
    " GMRES iterations, mean ", real(Method_Mb_OrbBased_cayleyNIterations, R64) / &
    real(max(Method_Mb_OrbBased_cayleyNSolves, 1_I64), R64), ", max ", groundMaxIterations

  ! 3. The kick at dt with a fresh reference (built after the kick)
  call RestartFromGroundState
  call RunKick(energyGround, energyKick, slopeExact, slopePropagated, maxEnergyDev, maxNormDev, maxTransverse, &
               times, dipoleZ, secondsPerStep, kappa, dt, endTime, 10, "dipoleKickMeanField.dat")
  write (*, '(1X, A, I0, A, I0, A, F6.2, A, I0)') "Cayley GMRES of the kick: ", Method_Mb_OrbBased_cayleyNSolves, &
    " solves, ", Method_Mb_OrbBased_cayleyNIterations, " iterations, mean ", &
    real(Method_Mb_OrbBased_cayleyNIterations, R64) / real(max(Method_Mb_OrbBased_cayleyNSolves, 1_I64), R64), &
    ", max ", Method_Mb_OrbBased_cayleyMaxIterations
  if (.not. (abs(energyGroundRef - energyGround) <= 1e-10_R64 * abs(energyGround))) &
    error stop "T_Fem_09_H2KickMeanFieldReference: E_el(0) differs from that of dipoleKick.dat (orbital or mesh &
               &changed since T_Fem_06_H2Kick; run it again)"

  maxResponse = 0.0_R64
  maxDiffRef = 0.0_R64
  maxTimeDev = 0.0_R64
  do iStep = 1, nSteps
    maxTimeDev = max(maxTimeDev, abs(times(iStep) - timesRef(iStep)))
    maxResponse = max(maxResponse, abs(dipoleZRef(iStep) - dipoleZRef(0)))
    maxDiffRef = max(maxDiffRef, abs((dipoleZ(iStep) - dipoleZ(0)) - (dipoleZRef(iStep) - dipoleZRef(0))))
  end do

  ! 4. The kick at dt/2 and dt/4 up to refineEndTime: d(dt), d(dt/2), d(dt/4)
  !    at that time
  dipoleAt(1) = dipoleZ(nint(refineEndTime / dt)) - dipoleZ(0)
  responseRef = dipoleZRef(nint(refineEndTime / dt)) - dipoleZRef(0)
  maxRefineEnergyDev = 0.0_R64
  do iRefine = 1, 2
    call RestartFromGroundState
    write (fileName, '(A, I0, A)') "dipoleKickMeanFieldRefine", iRefine, ".dat"
    call RunKick(refineEnergyGround, refineEnergyKick, refineSlopeExact, refineSlopePropagated, refineEnergyDev, &
                 refineNormDev, refineTransverse, refineTimes, refineDipoleZ, refineSecondsPerStep, kappa, &
                 refineTimeSteps(iRefine), refineEndTime, 10, trim(fileName))
    dipoleAt(iRefine + 1) = refineDipoleZ(ubound(refineDipoleZ, 1)) - refineDipoleZ(0)
    maxRefineEnergyDev = max(maxRefineEnergyDev, refineEnergyDev)
  end do
  order = log(abs(dipoleAt(1) - dipoleAt(2)) / abs(dipoleAt(2) - dipoleAt(3))) / log(2.0_R64)
  richardson = dipoleAt(3) + (dipoleAt(3) - dipoleAt(2)) / 3.0_R64
  errorMeanField = dipoleAt(1) - richardson
  errorPlain = responseRef - richardson

  ! 5. The nonlinear derivative at the kicked state with the reference built
  !    before the kick from the ground state, then with the default reference
  !    (released by the restart, built from the kicked state at the call)
  call RestartFromGroundState
  call Method_Mb_OrbBased_FillRdm1(rdm1, Method_state)
  call Method_Mb_OrbBased_BuildFieldFreeMeanField(Orbs_orbs, rdm1, 0.0_R64)
  call KickAllOrbitals(kappa, [0.0_R64, 0.0_R64, 1.0_R64])
  call Method_Mb_OrbBased_TimeDerivativeCoeffsPlusOrbsNonLin(dState, Method_state, 0.0_R64)
  derivativeExplicit = OrbitalNorm(dState) / derivativePlain
  call RestartFromGroundState
  call KickAllOrbitals(kappa, [0.0_R64, 0.0_R64, 1.0_R64])
  call Method_Mb_OrbBased_TimeDerivativeCoeffsPlusOrbsNonLin(dState, Method_state, 0.0_R64)
  derivativeDefault = OrbitalNorm(dState) / derivativePlain

  !==========================================
  call Say_Section("printout results")
  !==========================================

  kickEnergyDev = (energyKick - energyGround) / kappa**2 - 1.0_R64
  trkDev = slopeExact / (2.0_R64 * kappa) - 1.0_R64
  slopeDev = slopePropagated / slopeExact - 1.0_R64

  print *
  print *, "1. nonlinear orbital derivative with/without reference:", derivativeRatio
  print *, "   max |CI derivative| with reference:          ", derivativeCoeffs
  print *, "2. ground state, dt:                            ", groundTimeStep
  print *, "   max |E(t) - E(0)|:                           ", groundEnergyDev
  print *, "   max 1 - |<phi0|phi(t)>|:                     ", groundOverlapDev
  print *, "   max |N(t) - 2|:                              ", groundNormDev
  print *, "3. kick kappa, dt:                              ", kappa, dt
  print *, "   E(0+) - E(0), relative deviation:            ", energyKick - energyGround, kickEnergyDev
  print *, "   TRK on the grid, relative deviation:         ", trkDev
  print *, "   propagated initial slope, relative deviation:", slopeDev
  print *, "   max |E(t) - E(0+)|:                          ", maxEnergyDev
  print *, "   max |N(t) - 2|:                              ", maxNormDev
  print *, "   max |<x>|, |<y>|:                            ", maxTransverse
  print *, "   max |d(t) - d_T_Fem_06(t)|, max |d|:         ", maxDiffRef, maxResponse
  print *, "     relative:                                  ", maxDiffRef / maxResponse
  print *, "4. d(t) at t =                                  ", refineEndTime
  print *, "   dt, dt/2, dt/4 (reference):                  ", dipoleAt
  print *, "   order p:                                     ", order
  print *, "   Richardson limit:                            ", richardson
  print *, "   error at dt with / without the reference:    ", errorMeanField, errorPlain
  print *, "     relative:                                  ", errorMeanField / richardson, errorPlain / richardson
  print *, "   max |E(t) - E(0+)| of the refined runs:      ", maxRefineEnergyDev
  print *, "5. nonlinear derivative at the kicked state, relative, U from the ground state / the kicked state:"
  print *, "                                                ", derivativeExplicit, derivativeDefault
  write (*, '(1X, A, F8.4, A, F8.4, A)') "time per step: ", secondsPerStep, " s (dt), ", refineSecondsPerStep, " s (dt/4)"

  call check(error, derivativeRatio < 1e-9_R64 .and. derivativeCoeffs < 1e-12_R64)
  if (allocated(error)) error stop "T_Fem_09_H2KickMeanFieldReference failure: the nonlinear part does not vanish &
                                  &at the ground state with the reference"

  call check(error, groundEnergyDev < 1e-10_R64 .and. groundOverlapDev < 1e-10_R64 .and. groundNormDev < 1e-10_R64)
  if (allocated(error)) error stop "T_Fem_09_H2KickMeanFieldReference failure: ground state not stationary at the &
                                  &large time step"

  call check(error, groundMaxIterations <= MAX_ITERATIONS_GROUND)
  if (allocated(error)) error stop "T_Fem_09_H2KickMeanFieldReference failure: too many GMRES iterations with the &
                                  &reference"

  call check(error, abs(kickEnergyDev) < 1e-5_R64 .and. abs(trkDev) < 1e-5_R64 .and. abs(slopeDev) < 1e-4_R64)
  if (allocated(error)) error stop "T_Fem_09_H2KickMeanFieldReference failure: kick energy or initial slope (TRK)"

  call check(error, maxEnergyDev < 1e-8_R64 .and. maxNormDev < 1e-9_R64 .and. maxTransverse < 1e-9_R64)
  if (allocated(error)) error stop "T_Fem_09_H2KickMeanFieldReference failure: energy, particle number or &
                                  &transverse dipole not conserved"

  call check(error, maxTimeDev < 1e-12_R64 .and. maxDiffRef < dipoleTolerance * maxResponse)
  if (allocated(error)) error stop "T_Fem_09_H2KickMeanFieldReference failure: dipole differs from T_Fem_06 &
                                  &(same dt, no reference)"

  call check(error, order > 1.7_R64 .and. order < 2.3_R64 .and. abs(errorMeanField) < refineTolerance * abs(richardson))
  if (allocated(error)) error stop "T_Fem_09_H2KickMeanFieldReference failure: the splitting with the reference does &
                                  &not converge at second order"

  call check(error, derivativeExplicit > MIN_RESPONSE_DERIVATIVE .and. derivativeDefault < 1e-9_R64)
  if (allocated(error)) error stop "T_Fem_09_H2KickMeanFieldReference failure: a reference built before the kick is &
                                  &not kept, or the default one is not built from the kicked state"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  !> Ground state again (orbitals and CI coefficient) with a fresh Cayley
  !> setup: counters reset, reference released (rebuilt from the state at
  !> the next propagation, i.e. after the kick of RunKick).
  subroutine RestartFromGroundState
    Method_state(:) = groundState(:)
    call Method_Mb_OrbBased_SetupFieldFree
  end subroutine

  !> Norm of the orbital part of a state vector in the grid metric.
  function OrbitalNorm(state) result(res)
    complex(R64), intent(in), contiguous :: state(:)
    real(R64) :: res

    integer(I32) :: iOrb, offset

    res = 0.0_R64
    do iOrb = 1, Orbs_nOrbsInState
      offset = Coeffs_nCoeffs + (iOrb - 1) * Grid_nPoints
      res = res + real(Grid_InnerProduct(state(offset + 1:offset + Grid_nPoints), state(offset + 1:offset + Grid_nPoints)), R64)
    end do
    res = sqrt(res)

  end function

  include 'FemKickFile.inc'
  include 'FemKick.inc'
  include 'FemDynamics.inc'

end program
