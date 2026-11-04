!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Kick helpers shared by T_Fem_06_H2Kick, T_Fem_07_H2KickDouble and
! T_Fem_09_H2KickMeanFieldReference (included after `contains`, together with
! FemDynamics.inc; every procedure here is used by each of them).
!
! - RunKick:          kick the loaded ground state and propagate it
! - KickAllOrbitals:  the dipole kick exp(i kappa e.r) of every stored orbital
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Applies the dipole kick phi <- exp(i kappa z) phi to every orbital of the
!> loaded TDHF ground state (KickAllOrbitals: every electron gets the momentum
!> kappa along z) and propagates it with dt up to endTime. The dipole is
!> written at every step to dipoleFile ("time <x> <y> <z>",
!> PrinterObservableFem_DumpDipole),
!> the energy and N = Tr(rho S) every outputEvery steps. The file is replaced
!> at the start; its first line is the header of the run,
!>   "# T_Fem kick (<date> <time>): kappa, timeStep, endTime, E_el(0),
!>    nSteps, nPoints = <4 reals> <2 integers>"
!> (one line), which T_Fem_07_H2KickDouble checks before it compares with the
!> dipole of T_Fem_06_H2Kick.
!>
!> Returns
!>   energyGround     E_el of the loaded ground state
!>   energyKick       E_el right after the kick
!>   slopeExact       d<z>/dt at t = 0+ of the space-discretized TDHF,
!>                    i Tr(rho [T, z]) (PrinterObservableFem_DumpDipoleVelocity;
!>                    V and the mean field are multiplicative and commute with
!>                    z)
!>   slopePropagated  d<z>/dt at t = 0 from the propagated dipole: Richardson
!>                    extrapolation of (<z>(t) - <z>(0))/t from the first two
!>                    steps (the force vanishes at t = 0, so <z>(t) - <z>(0) is
!>                    odd in t up to the splitting error)
!>   maxEnergyDev     max |E(t) - energyKick| at the outputs
!>   maxNormDev       max |N(t) - 2|
!>   maxTransverse    max |<x>|, |<y>| over all steps
!>   times, dipoleZ   t and <z>(t) - <z>(0) at every step (0:nSteps)
!>   secondsPerStep   wall-clock time per step (dipole output included)
subroutine RunKick(energyGround, energyKick, slopeExact, slopePropagated, maxEnergyDev, maxNormDev, maxTransverse, &
                   times, dipoleZ, secondsPerStep, kappa, dt, endTime, outputEvery, dipoleFile)
  use M_Utils_PrinterObservableFem
  use M_Grid
  use M_Orbs
  use M_Method
  use M_Method_Mb_OrbBased
  use M_Propagator

  real(R64), intent(out) :: energyGround
  real(R64), intent(out) :: energyKick
  real(R64), intent(out) :: slopeExact
  real(R64), intent(out) :: slopePropagated
  real(R64), intent(out) :: maxEnergyDev
  real(R64), intent(out) :: maxNormDev
  real(R64), intent(out) :: maxTransverse
  real(R64), intent(out), allocatable :: times(:)
  real(R64), intent(out), allocatable :: dipoleZ(:)
  real(R64), intent(out) :: secondsPerStep
  real(R64), intent(in) :: kappa
  real(R64), intent(in) :: dt
  real(R64), intent(in) :: endTime
  integer(I32), intent(in) :: outputEvery
  character(len=*), intent(in) :: dipoleFile

  complex(R64), allocatable :: rdm1(:, :)
  real(R64) :: energy, norm, dipole(3), dipoleStart(3), velocity(3), time, slope1, slope2
  integer(I64) :: count0
  integer(I32) :: nSteps, iStep, io, dateTime(8)

  nSteps = nint(endTime / dt)
  allocate (times(0:nSteps), dipoleZ(0:nSteps))

  call FemObservables(energyGround, norm, dipoleStart, 0.0_R64, "")

  call KickAllOrbitals(kappa, [0.0_R64, 0.0_R64, 1.0_R64])

  ! d<z>/dt at 0+ of the semi-discrete TDHF, i Tr(rho [T, z]) (for the
  ! doubly occupied psi 4 Im <z psi|T psi>)
  call Method_Mb_OrbBased_FillRdm1(rdm1, Method_state)
  call PrinterObservableFem_DumpDipoleVelocity(rdm1, Orbs_orbs, "", .false., 0.0_R64, velocity_=velocity)
  slopeExact = velocity(3)

  ! Fresh dipole file: the header of the run, then the first line at t = 0
  ! (after the kick)
  call date_and_time(values=dateTime)
  open (newunit=io, file=dipoleFile, status="replace")
  write (io, '(A, I4.4, "-", I2.2, "-", I2.2, " ", I2.2, ":", I2.2, ":", I2.2, A, 4ES24.15E3, 2(1X, I0))') &
    "# T_Fem kick (", dateTime(1:3), dateTime(5:7), "): kappa, timeStep, endTime, E_el(0), nSteps, nPoints = ", &
    kappa, dt, endTime, energyGround, nSteps, Grid_nPoints
  close (io)
  call FemObservables(energyKick, norm, dipole, 0.0_R64, dipoleFile)
  times(0) = 0.0_R64
  dipoleZ(0) = dipole(3) - dipoleStart(3)
  maxEnergyDev = 0.0_R64
  maxNormDev = abs(norm - 2.0_R64)
  maxTransverse = maxval(abs(dipole(1:2)))

  time = 0.0_R64
  call system_clock(count0)
  do iStep = 1, nSteps
    call Propagator_Propagate(Method_state, time, time + dt)
    time = iStep * dt

    if (mod(iStep, outputEvery) .eq. 0 .or. iStep .eq. nSteps) then
      call FemObservables(energy, norm, dipole, time, dipoleFile)
      maxEnergyDev = max(maxEnergyDev, abs(energy - energyKick))
      maxNormDev = max(maxNormDev, abs(norm - 2.0_R64))
      write (*, '(1X, A, F6.3, A, F18.13, A, ES10.2, A, ES10.2, A, 3ES13.5)') "t = ", time, ": E_el ", energy, &
        ", E - E(0+) ", energy - energyKick, ", N - 2 ", norm - 2.0_R64, ", dipole ", dipole
    else
      call Method_Mb_OrbBased_FillRdm1(rdm1, Method_state)
      call PrinterObservableFem_DumpDipole(rdm1, Orbs_orbs, dipoleFile, .false., time, dipole_=dipole)
    end if
    times(iStep) = time
    dipoleZ(iStep) = dipole(3) - dipoleStart(3)
    maxTransverse = max(maxTransverse, maxval(abs(dipole(1:2))))
  end do
  secondsPerStep = ElapsedSeconds(count0) / nSteps

  ! Richardson extrapolation of <z>(t)/t = v0 + c t^2 + ... to t = 0
  slope1 = (dipoleZ(1) - dipoleZ(0)) / times(1)
  slope2 = (dipoleZ(2) - dipoleZ(0)) / times(2)
  slopePropagated = (4.0_R64 * slope1 - slope2) / 3.0_R64

end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> The dipole kick phi <- exp(i kappa e.r) phi of every stored orbital at the
!> nodes (e = direction, a unit vector; e.r from Grid_Fem_FillNodeCoordinate):
!> a one-body unitary on the many-body state that keeps the CI coefficients
!> and gives every electron the momentum kappa along e (Orbs_orbs is the
!> orbital part of Method_state).
subroutine KickAllOrbitals(kappa, direction)
  use M_Grid_Fem
  use M_Orbs

  real(R64), intent(in) :: kappa
  real(R64), intent(in) :: direction(3)

  real(R64), allocatable :: coordinate(:)
  integer(I32) :: iOrb

  call Grid_Fem_FillNodeCoordinate(coordinate, direction)
  do iOrb = 1, Orbs_nOrbsInState
    Orbs_orbs(:, iOrb) = exp(cmplx(0.0_R64, kappa * coordinate(:), kind=R64)) * Orbs_orbs(:, iOrb)
  end do

end subroutine
