!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> T_Fem_06_H2Kick with the doubled kick kappa = program.kappa = 2e-3 (same
!> TDHF, propagator, dt = 0.02 and end time 1.0; one run per JSON keeps each
!> test well under a minute) and the linearity of the dipole response: reads
!> the dipole of T_Fem_06 (dipoleKick.dat, written by
!> PrinterObservableFem_DumpDipole; run T_Fem_06 first, as ctest does) and
!> requires <z>_2kappa(t) - <z>(0) = 2 (<z>_kappa(t) - <z>(0)) at every step
!> (relative 1e-3; TDHF is nonlinear, but the even orders of the response
!> vanish by symmetry, so the deviation is O(kappa^2)).
!>
!> The file must be that of T_Fem_06 for this configuration, never stale or
!> foreign data: its header (FemKick.inc, RunKick) must give kappa/2 and the
!> same timeStep, endTime, nSteps and mesh size (checked before the
!> propagation), the same ground-state energy E_el(0) (after it; it
!> fingerprints the loaded orbital) and nSteps + 1 dipole lines (a complete
!> run). The header is printed, with the date of the T_Fem_06 run.
!>
!> Also repeats the single-run checks of T_Fem_06 for this kappa: kick energy
!> N kappa^2/2, TRK initial slope N kappa (semi-discrete and propagated),
!> <x> = <y> = 0, energy and particle number conserved.
program T_Fem_07_H2KickDouble
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Grid
  use M_Method_Mb_OrbBased
  use testdrive, only: check, error_type

  implicit none

  type(error_type), allocatable :: error
  real(R64), allocatable :: times(:), dipoleZ(:), timesRef(:), dipoleZRef(:)
  real(R64) :: kappa, dt, endTime, energyGround, energyKick, slopeExact, slopePropagated
  real(R64) :: maxEnergyDev, maxNormDev, maxTransverse, secondsPerStep, kickEnergyDev, trkDev, slopeDev
  real(R64) :: ratio, maxRatioDev, maxTimeDev, kappaRef, dtRef, endTimeRef, energyGroundRef
  integer(I32) :: outputEvery, iStep, nSteps, nStepsRef, nPointsRef
  character(len=:), allocatable :: headerRef

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call SetupFemDynamics("test/simulationTests/Fem/T_Fem_07_H2KickDouble.json")

  call Say_Fabricate("program")
  kappa = Json_Get("program.kappa", 2e-3_R64)
  dt = Json_Get("program.timeStep", 0.02_R64)
  endTime = Json_Get("program.endTime", 1.0_R64)
  outputEvery = Json_Get("program.outputEvery", 10)
  nSteps = nint(endTime / dt)

  !==========================================
  call Say_Section("perform action")
  !==========================================

  ! Response of T_Fem_06 (kappa/2): its run parameters first, before the
  ! propagation
  call ReadKickFile(headerRef, kappaRef, dtRef, endTimeRef, energyGroundRef, nStepsRef, nPointsRef, timesRef, &
                    dipoleZRef, "dipoleKick.dat")
  write (*, '(1X, A, A)') "reference (T_Fem_06_H2Kick): ", headerRef
  if (.not. (abs(kappaRef - 0.5_R64 * kappa) <= 1e-12_R64 * kappa)) &
    error stop "T_Fem_07_H2KickDouble: dipoleKick.dat is not from a kick with kappa/2 (stale or foreign file; &
               &run T_Fem_06_H2Kick first)"
  if (.not. (abs(dtRef - dt) <= 1e-12_R64 * dt .and. abs(endTimeRef - endTime) <= 1e-12_R64 * endTime &
             .and. nStepsRef .eq. nSteps .and. nPointsRef .eq. Grid_nPoints)) &
    error stop "T_Fem_07_H2KickDouble: dipoleKick.dat is from other steps or another mesh (stale or foreign file; &
               &run T_Fem_06_H2Kick first)"
  if (size(timesRef) .ne. nSteps + 1) &
    error stop "T_Fem_07_H2KickDouble: dipoleKick.dat is incomplete (T_Fem_06_H2Kick did not finish its run)"

  call RunKick(energyGround, energyKick, slopeExact, slopePropagated, maxEnergyDev, maxNormDev, maxTransverse, &
               times, dipoleZ, secondsPerStep, kappa, dt, endTime, outputEvery, "dipoleKickDouble.dat")

  ! Same ground state as in T_Fem_06 (loaded orbital and mesh)
  if (.not. (abs(energyGroundRef - energyGround) <= 1e-10_R64 * abs(energyGround))) &
    error stop "T_Fem_07_H2KickDouble: E_el(0) differs from that of dipoleKick.dat (orbital or mesh changed since &
               &T_Fem_06_H2Kick; run it again)"
  maxTimeDev = maxval(abs(timesRef - times))

  ! Linearity: the responses relative to the <z> at t = 0 of each run
  maxRatioDev = 0.0_R64
  do iStep = 1, ubound(times, 1)
    ratio = (dipoleZ(iStep) - dipoleZ(0)) / (dipoleZRef(iStep) - dipoleZRef(0))
    maxRatioDev = max(maxRatioDev, abs(ratio - kappa / kappaRef))
  end do

  !==========================================
  call Say_Section("printout results")
  !==========================================

  kickEnergyDev = (energyKick - energyGround) / kappa**2 - 1.0_R64
  trkDev = slopeExact / (2.0_R64 * kappa) - 1.0_R64
  slopeDev = slopePropagated / slopeExact - 1.0_R64

  print *
  print *, "kappa:                                ", kappa
  print *, "E(0), E(0+) - E(0), N kappa^2/2:      ", energyGround, energyKick - energyGround, kappa**2
  print *, "  relative deviation:                 ", kickEnergyDev
  print *, "d<z>/dt(0+): semi-discrete, N kappa:  ", slopeExact, 2.0_R64 * kappa
  print *, "  relative deviation (TRK on the grid):", trkDev
  print *, "d<z>/dt(0) of the propagation:        ", slopePropagated
  print *, "  relative deviation from semi-discrete:", slopeDev
  print *, "<z>(t_end) - <z>(0), T_Fem_06 value:  ", dipoleZ(ubound(dipoleZ, 1)) - dipoleZ(0), &
    dipoleZRef(ubound(dipoleZRef, 1)) - dipoleZRef(0)
  print *, "E_el(0), T_Fem_06 value:              ", energyGround, energyGroundRef
  print *, "max |ratio to T_Fem_06 - 2|:          ", maxRatioDev
  print *, "max |E(t) - E(0+)|:                   ", maxEnergyDev
  print *, "max |N(t) - 2|:                       ", maxNormDev
  print *, "max |<x>|, |<y>|:                     ", maxTransverse
  write (*, '(1X, A, I0, A, I0, A, F6.2, A, I0)') "Cayley GMRES during the propagation: ", &
    Method_Mb_OrbBased_cayleyNSolves, " solves, ", Method_Mb_OrbBased_cayleyNIterations, " iterations, mean ", &
    real(Method_Mb_OrbBased_cayleyNIterations, R64) / real(max(Method_Mb_OrbBased_cayleyNSolves, 1_I64), R64), &
    ", max ", Method_Mb_OrbBased_cayleyMaxIterations
  write (*, '(1X, A, I0, A, F8.4, A)') "propagation: ", size(times) - 1, " steps, ", secondsPerStep, " s per step"

  call check(error, maxTimeDev < 1e-12_R64)
  if (allocated(error)) error stop "T_Fem_07_H2KickDouble failure: times differ from T_Fem_06 (same dt and endTime needed)"

  call check(error, maxRatioDev < 1e-3_R64)
  if (allocated(error)) error stop "T_Fem_07_H2KickDouble failure: dipole response not linear in kappa"

  call check(error, abs(kickEnergyDev) < 1e-5_R64)
  if (allocated(error)) error stop "T_Fem_07_H2KickDouble failure: kick energy is not N kappa^2/2"

  call check(error, abs(trkDev) < 1e-5_R64)
  if (allocated(error)) error stop "T_Fem_07_H2KickDouble failure: initial slope d<z>/dt is not N kappa (TRK)"

  call check(error, abs(slopeDev) < 1e-4_R64)
  if (allocated(error)) error stop "T_Fem_07_H2KickDouble failure: propagated initial slope off the semi-discrete value"

  call check(error, maxTransverse < 1e-9_R64)
  if (allocated(error)) error stop "T_Fem_07_H2KickDouble failure: transverse dipole <x>, <y> not zero"

  call check(error, maxEnergyDev < 1e-8_R64)
  if (allocated(error)) error stop "T_Fem_07_H2KickDouble failure: energy not conserved after the kick"

  call check(error, maxNormDev < 1e-9_R64)
  if (allocated(error)) error stop "T_Fem_07_H2KickDouble failure: particle number not conserved"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  include 'FemKickFile.inc'
  include 'FemKick.inc'
  include 'FemDynamics.inc'

end program
