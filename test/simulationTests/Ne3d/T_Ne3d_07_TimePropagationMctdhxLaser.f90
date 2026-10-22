!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Real-time MCTDHF propagation of Ne (3d, Ylm grid) in a z-polarized
!> length-gauge laser pulse, computing the time-dependent dipole moment.
!>
!> The MCSCF ground state is converged first (laser is zero at t = 0, so the
!> coulombLaser potential reduces to the plain Coulomb potential of
!> T_Ne3d_06_GroundSolverMcscfYlm). The state is then propagated in real time
!> through the rising edge of a sin^2 pulse with the split-step (Strang)
!> propagator: the stiff linear one-body part is integrated by the Short
!> Iterative Lanczos in the FEDVR-weighted metric, the nonlinear
!> mean-field/CI part by classical RK4.
!>
!> Checks:
!> - the MCSCF ground-state energy matches T_Ne3d_06 (laser off at t = 0)
!> - the weighted-metric norm of the state is conserved by the propagation
!> - the dipole moment <z> responds to the field (nonzero, regression value)
program T_Ne3d_07_TimePropagationMctdhxLaser
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

  call SetupMctdhxLaser("test/simulationTests/Ne3d/T_Ne3d_07_TimePropagationMctdhxLaser.json")

  !==========================================
  call Say_Section("converge mcscf ground state")
  !==========================================

  energy = ConvergedGroundStateEnergy()

  ! Laser is off at t = 0: must reproduce the T_Ne3d_06 mcscf energy.
  call check(error, energy, -127.7449896_R64, thr=1e-7_R64)
  if (allocated(error)) error stop "T_Ne3d_07_TimePropagationMctdhxLaser failure: ground state off reference"

  !==========================================
  call Say_Section("propagate through laser pulse")
  !==========================================

  startNorm = MetricNorm(Method_state)

  time = Json_Get("program.startTime", 0.0_R64)
  do while (BeforeEndTimeQ(time))
    call PropagateOutputStep(time, dipole, "dipole.dat")
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
  if (allocated(error)) error stop "T_Ne3d_07_TimePropagationMctdhxLaser failure: norm not conserved"

  ! The dipole moment must respond to the field (regression value).
  call check(error, dipole, 4.4862488373e-4_R64, thr=1e-8_R64)
  if (allocated(error)) error stop "T_Ne3d_07_TimePropagationMctdhxLaser failure: dipole off reference"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  include 'Ne3dMctdhxLaser.inc'

end program
