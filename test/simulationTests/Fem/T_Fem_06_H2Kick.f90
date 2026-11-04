!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Dipole kick of H2 (R = 1.4) on the finite-element grid (grid.fem): TDHF
!> (restricted MCTDHX, nOrbs = nBodies = [1, 1]) from the RHF orbital of
!> T_Fem_03, kicked by phi <- exp(i kappa z) phi with kappa = program.kappa =
!> 1e-3 along the molecular axis, then propagated field-free as in T_Fem_05
!> (splitStep order2, Cayley exact part with exactPartScheme strang, RK4) with
!> dt = program.timeStep = 0.02 up to program.endTime = 1.0. The dipole
!> (<x>, <y>, <z>) is written at every step with PrinterObservableFem_DumpDipole
!> to dipoleKick.dat (replaced at the start, first line a header with kappa,
!> the step parameters, E_el(0) and the mesh size), which
!> T_Fem_07_H2KickDouble (kappa doubled) checks and reads for the linearity
!> check.
!>
!> Checks (N = 2 electrons, mass 1):
!> - kick energy: E(0+) - E(0) = N kappa^2/2 (relative 1e-5; the grid value of
!>   <[z,[T,z]]> differs from 1/m by the discretization error, -1.3e-6 here)
!> - TRK / Ehrenfest: d<z>/dt at t = 0+ = <p_z>/m = N kappa. The space-
!>   discretized TDHF gives d<z>/dt = i Tr(rho [T, z]) = 4 Im<z psi|T psi>
!>   exactly (PrinterObservableFem_DumpDipoleVelocity; V and the mean field
!>   commute with z); it matches N kappa to the discretization error
!>   (relative 1e-5; -1.3e-6 here)
!> - the propagated initial slope (Richardson extrapolation of the first two
!>   steps) equals that of the semi-discrete equations to the time
!>   discretization error (relative 1e-4; -1.5e-5 at dt = 0.02, -3.2e-6 at
!>   0.01, -6.7e-7 at 0.005: second order)
!> - <x> = <y> = 0 (kick along the axis, cylindrical symmetry; 1e-9, the
!>   loaded orbital carries ~5e-11 from the SCF convergence), the energy after
!>   the kick conserved (1e-8) and N = Tr(rho S) = 2 (1e-9)
!> Prints the GMRES iteration counts of the Cayley solves and the time per
!> step.
program T_Fem_06_H2Kick
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Method_Mb_OrbBased
  use testdrive, only: check, error_type

  implicit none

  type(error_type), allocatable :: error
  real(R64), allocatable :: times(:), dipoleZ(:)
  real(R64) :: kappa, dt, endTime, energyGround, energyKick, slopeExact, slopePropagated
  real(R64) :: maxEnergyDev, maxNormDev, maxTransverse, secondsPerStep, kickEnergyDev, trkDev, slopeDev
  integer(I32) :: outputEvery

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call SetupFemDynamics("test/simulationTests/Fem/T_Fem_06_H2Kick.json")

  call Say_Fabricate("program")
  kappa = Json_Get("program.kappa", 1e-3_R64)
  dt = Json_Get("program.timeStep", 0.02_R64)
  endTime = Json_Get("program.endTime", 1.0_R64)
  outputEvery = Json_Get("program.outputEvery", 10)

  !==========================================
  call Say_Section("perform action")
  !==========================================

  call RunKick(energyGround, energyKick, slopeExact, slopePropagated, maxEnergyDev, maxNormDev, maxTransverse, &
               times, dipoleZ, secondsPerStep, kappa, dt, endTime, outputEvery, "dipoleKick.dat")

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
  print *, "<z>(t_end) - <z>(0):                  ", dipoleZ(ubound(dipoleZ, 1))
  print *, "max |E(t) - E(0+)|:                   ", maxEnergyDev
  print *, "max |N(t) - 2|:                       ", maxNormDev
  print *, "max |<x>|, |<y>|:                     ", maxTransverse
  write (*, '(1X, A, I0, A, I0, A, F6.2, A, I0)') "Cayley GMRES during the propagation: ", &
    Method_Mb_OrbBased_cayleyNSolves, " solves, ", Method_Mb_OrbBased_cayleyNIterations, " iterations, mean ", &
    real(Method_Mb_OrbBased_cayleyNIterations, R64) / real(max(Method_Mb_OrbBased_cayleyNSolves, 1_I64), R64), &
    ", max ", Method_Mb_OrbBased_cayleyMaxIterations
  write (*, '(1X, A, I0, A, F8.4, A)') "propagation: ", size(times) - 1, " steps, ", secondsPerStep, " s per step"

  call check(error, abs(kickEnergyDev) < 1e-5_R64)
  if (allocated(error)) error stop "T_Fem_06_H2Kick failure: kick energy is not N kappa^2/2"

  call check(error, abs(trkDev) < 1e-5_R64)
  if (allocated(error)) error stop "T_Fem_06_H2Kick failure: initial slope d<z>/dt is not N kappa (TRK)"

  call check(error, abs(slopeDev) < 1e-4_R64)
  if (allocated(error)) error stop "T_Fem_06_H2Kick failure: propagated initial slope off the semi-discrete value"

  call check(error, maxTransverse < 1e-9_R64)
  if (allocated(error)) error stop "T_Fem_06_H2Kick failure: transverse dipole <x>, <y> not zero"

  call check(error, maxEnergyDev < 1e-8_R64)
  if (allocated(error)) error stop "T_Fem_06_H2Kick failure: energy not conserved after the kick"

  call check(error, maxNormDev < 1e-9_R64)
  if (allocated(error)) error stop "T_Fem_06_H2Kick failure: particle number not conserved"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  include 'FemKick.inc'
  include 'FemDynamics.inc'

end program
