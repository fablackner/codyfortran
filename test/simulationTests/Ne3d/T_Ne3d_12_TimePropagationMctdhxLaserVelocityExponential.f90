!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> T_Ne3d_09_TimePropagationMctdhxLaserVelocityExactFieldFree with the
!> exponential Runge–Kutta integrator (propagator.exponential.krogstad4): the
!> field-free h0 enters through the phi functions of its eigenbasis
!> (Method_ApplyExactPartPhi), the gauge coupling and the interaction through
!> N, without operator splitting.
!>
!> Checks:
!> - the MCSCF ground-state energy matches T_Ne3d_06/07/08 (gauge off at t = 0)
!> - the weighted-metric norm is conserved up to the integration error
!>   (Krogstad is not unitary; the drift falls ~35x per halved step)
!> - the dipole moment agrees with the converged value 1.3140248e-3 (Krogstad
!>   at a 4x smaller step and the Richardson extrapolation of T_Ne3d_10 at 4x
!>   and 8x smaller steps agree to 3e-10); at this step the error is
!>   5e-8, against 4.3e-7 for the split-step schemes of T_Ne3d_09/10, whose
!>   outer splitting of the linear and nonlinear parts is second order
program T_Ne3d_12_TimePropagationMctdhxLaserVelocityExponential
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

  call SetupMctdhxLaser("test/simulationTests/Ne3d/T_Ne3d_12_TimePropagationMctdhxLaserVelocityExponential.json")

  !==========================================
  call Say_Section("converge mcscf ground state")
  !==========================================

  energy = ConvergedGroundStateEnergy()

  ! Gauge coupling is off at t = 0: must reproduce the T_Ne3d_06/07 mcscf energy.
  call check(error, energy, -127.7449896_R64, thr=1e-7_R64)
  if (allocated(error)) error stop "T_Ne3d_12_TimePropagationMctdhxLaserVelocityExponential failure: ground state off reference"

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

  ! Krogstad conserves the norm only up to its integration error
  call check(error, finalNorm, startNorm, thr=3e-6_R64)
  if (allocated(error)) error stop "T_Ne3d_12_TimePropagationMctdhxLaserVelocityExponential failure: norm not conserved"

  ! Converged dipole (see header): the error at this step is 5e-8
  call check(error, dipole, 1.3140248e-3_R64, thr=1e-7_R64)
  if (allocated(error)) error stop "T_Ne3d_12_TimePropagationMctdhxLaserVelocityExponential failure: dipole off converged value"

  call check(error, dipole, 1.3140769558e-3_R64, thr=1e-11_R64)
  if (allocated(error)) error stop "T_Ne3d_12_TimePropagationMctdhxLaserVelocityExponential failure: dipole off reference"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  include 'Ne3dMctdhxLaser.inc'

end program
