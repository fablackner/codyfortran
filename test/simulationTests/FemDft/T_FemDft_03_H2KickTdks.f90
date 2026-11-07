!> ctest-requires: libxc
!> Real-time TDDFT (method.mb.orbBased.tdks, adiabatic LDA) of H2 on the mesh of
!> T_Fem_03 from the Kohn-Sham ground state of T_FemDft_02 (orb01_01.in), with
!> the split-step propagator, the Cayley exact field-free part and the frozen
!> Kohn-Sham reference (exactFieldFreeMeanFieldQ): h0 + U exactly, U =
!> v_KS[rho(0+)] (= v_KS of the ground state, the kick keeps the density), and
!> v_KS[rho(t)] - U by RK4, which vanishes at the ground state:
!>  - kick e^(iκz): E(0+) - E(0) = N κ²/2 (the density, hence E_H and E_xc,
!>    is unchanged by the phase) and the TRK initial slope d<z>/dt = N κ;
!>  - the propagated initial slope, the Kohn-Sham energy (the equations are
!>    its Hamiltonian flow; 2e-11 over 0.6 a.u. at dt = 0.02, 2e-7 with the
!>    plain split h0 | v_KS), the particle number and the vanishing
!>    transverse dipole over 0.6 a.u. (the dipole follows the plain split at
!>    dt = 0.005 to 5e-6 of the response).
program T_FemDft_03_H2KickTdks
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

  call Say_Hello
  call SetupFemDynamics("test/simulationTests/FemDft/T_FemDft_03_H2KickTdks.json")
  kappa = Json_Get("program.kappa", 1e-3_R64)
  dt = Json_Get("program.timeStep", 0.02_R64)
  endTime = Json_Get("program.endTime", 0.6_R64)
  outputEvery = Json_Get("program.outputEvery", 10)

  call RunKick(energyGround, energyKick, slopeExact, slopePropagated, maxEnergyDev, maxNormDev, maxTransverse, &
               times, dipoleZ, secondsPerStep, kappa, dt, endTime, outputEvery, "dipoleKickTdks.dat")

  kickEnergyDev = (energyKick - energyGround) / kappa**2 - 1.0_R64
  trkDev = slopeExact / (2.0_R64 * kappa) - 1.0_R64
  slopeDev = slopePropagated / slopeExact - 1.0_R64

  print *
  print *, "E(0), E(0+) - E(0), N kappa^2/2:      ", energyGround, energyKick - energyGround, kappa**2
  print *, "  relative deviation:                 ", kickEnergyDev
  print *, "d<z>/dt(0+): semi-discrete, N kappa:  ", slopeExact, 2.0_R64 * kappa, ", relative deviation ", trkDev
  print *, "d<z>/dt(0) of the propagation:        ", slopePropagated, ", relative deviation ", slopeDev
  print *, "<z>(t_end) - <z>(0):                  ", dipoleZ(ubound(dipoleZ, 1))
  print *, "max |E(t) - E(0+)|:                   ", maxEnergyDev
  print *, "max |N(t) - 2|:                       ", maxNormDev
  print *, "max |<x>|, |<y>|:                     ", maxTransverse
  write (*, '(1X, A, I0, A, F8.4, A)') "propagation: ", size(times) - 1, " steps, ", secondsPerStep, " s per step"

  call check(error, abs(kickEnergyDev) < 1e-5_R64)
  if (allocated(error)) error stop "T_FemDft_03_H2KickTdks failure: kick energy is not N kappa^2/2"
  call check(error, abs(trkDev) < 1e-5_R64)
  if (allocated(error)) error stop "T_FemDft_03_H2KickTdks failure: initial slope d<z>/dt is not N kappa (TRK)"
  call check(error, abs(slopeDev) < 1e-4_R64)
  if (allocated(error)) error stop "T_FemDft_03_H2KickTdks failure: propagated initial slope"
  call check(error, maxTransverse < 1e-9_R64)
  if (allocated(error)) error stop "T_FemDft_03_H2KickTdks failure: transverse dipole <x>, <y> not zero"
  call check(error, maxEnergyDev < 1e-9_R64)
  if (allocated(error)) error stop "T_FemDft_03_H2KickTdks failure: Kohn-Sham energy not conserved"
  call check(error, maxNormDev < 1e-9_R64)
  if (allocated(error)) error stop "T_FemDft_03_H2KickTdks failure: particle number not conserved"

  call Say_Goodbye

contains

  include '../Fem/FemKick.inc'
  include '../Fem/FemDynamics.inc'

end program
