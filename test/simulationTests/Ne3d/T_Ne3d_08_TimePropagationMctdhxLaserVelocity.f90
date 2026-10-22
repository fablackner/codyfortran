!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Real-time MCTDHF propagation of Ne (3d, Ylm grid) in a z-polarized
!> velocity-gauge laser pulse, computing the time-dependent dipole moment.
!>
!> Velocity-gauge companion of T_Ne3d_07_TimePropagationMctdhxLaser: the laser
!> enters through the SysGauge coupling A(t)·p + A²/2 instead of the
!> length-gauge E(t)·z potential, while the static -Z/r attraction remains an
!> ordinary SysPotential. The MCSCF ground state is converged first (A = 0 at
!> t = 0, so the gauge term vanishes and the ground state must match
!> T_Ne3d_06/07 exactly). The state is then propagated in real time through
!> the rising edge of a sin² pulse with the split-step (Strang) propagator:
!> the stiff linear one-body part (now including A·p) is integrated by the
!> Short Iterative Lanczos in the FEDVR-weighted metric, the nonlinear
!> mean-field/CI part by classical RK4.
!>
!> Checks:
!> - the MCSCF ground-state energy matches T_Ne3d_06/07 (gauge off at t = 0)
!> - the weighted-metric norm of the state is conserved by the propagation
!>   (requires A·p Hermitian in the weighted metric, see SysGauge_Ylm)
!> - the dipole moment <z> responds to the field (nonzero, regression value)
program T_Ne3d_08_TimePropagationMctdhxLaserVelocity
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Orbs
  use M_Method
  use M_Method_Mb_OrbBased
  use testdrive, only: check, error_type

  implicit none

  real(R64)    :: time, energy, dipole
  real(R64)    :: startNorm, finalNorm
  type(error_type), allocatable :: error

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call SetupMctdhxLaser("test/simulationTests/Ne3d/T_Ne3d_08_TimePropagationMctdhxLaserVelocity.json")

  !==========================================
  call Say_Section("converge mcscf ground state")
  !==========================================

  energy = ConvergedGroundStateEnergy()

  ! Gauge coupling is off at t = 0: must reproduce the T_Ne3d_06/07 mcscf energy.
  call check(error, energy, -127.7449896_R64, thr=1e-7_R64)
  if (allocated(error)) error stop "T_Ne3d_08_TimePropagationMctdhxLaserVelocity failure: ground state off reference"

  !==========================================
  call Say_Section("propagate through laser pulse")
  !==========================================

  startNorm = MetricNorm(Method_state)

  time = Json_Get("program.startTime", 0.0_R64)
  do while (BeforeEndTimeQ(time))
    call PropagateOutputStep(time, dipole, "dipoleVelocity.dat")
  end do

  !==========================================
  call Say_Section("printout results")
  !==========================================

  finalNorm = MetricNorm(Method_state)

  print *
  print *, "start norm:   ", startNorm
  print *, "final norm:   ", finalNorm
  print *, "final dipole: ", dipole

  ! Unitary propagation must conserve the weighted-metric norm of the state
  ! (up to the RK4 integration error of the nonlinear part).
  call check(error, finalNorm, startNorm, thr=1e-7_R64)
  if (allocated(error)) error stop "T_Ne3d_08_TimePropagationMctdhxLaserVelocity failure: norm not conserved"

  ! The dipole moment must respond to the field (regression value). It is not
  ! comparable to T_Ne3d_07: for the 1-cycle pulse the envelope-derivative
  ! part of E = -dA/dt differs substantially from the length-gauge field.
  call check(error, dipole, 1.3142887110e-3_R64, thr=1e-8_R64)
  if (allocated(error)) error stop "T_Ne3d_08_TimePropagationMctdhxLaserVelocity failure: dipole off reference"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  include 'Ne3dMctdhxLaser.inc'

end program
