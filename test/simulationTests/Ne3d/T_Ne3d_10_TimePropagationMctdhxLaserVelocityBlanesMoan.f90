!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> T_Ne3d_09_TimePropagationMctdhxLaserVelocityExactFieldFree with the
!> fourth-order composition of the exact field-free propagation and the gauge
!> SIL inside every linear step (propagator.splitStep.exactPartScheme =
!> "blanesMoan6", Blanes & Moan 2002, 7 exact and 6 SIL substeps).
!>
!> Checks:
!> - the MCSCF ground-state energy matches T_Ne3d_06/07/08 (gauge off at t = 0)
!> - the weighted-metric norm is conserved
!> - the dipole moment agrees with T_Ne3d_08 within the time-step error: the
!>   linear step is now accurate to fourth order, the remaining difference is
!>   the midpoint (frozen A(t)) error of the SIL step of T_Ne3d_08 at this
!>   fast field (omega = 0.5); for a slowly varying field both agree ~100x
!>   better than the second-order (Strang) composition of T_Ne3d_09
program T_Ne3d_10_TimePropagationMctdhxLaserVelocityBlanesMoan
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

  call SetupMctdhxLaser("test/simulationTests/Ne3d/T_Ne3d_10_TimePropagationMctdhxLaserVelocityBlanesMoan.json")

  !==========================================
  call Say_Section("converge mcscf ground state")
  !==========================================

  energy = ConvergedGroundStateEnergy()

  ! Gauge coupling is off at t = 0: must reproduce the T_Ne3d_06/07 mcscf energy.
  call check(error, energy, -127.7449896_R64, thr=1e-7_R64)
  if (allocated(error)) error stop "T_Ne3d_10_TimePropagationMctdhxLaserVelocityBlanesMoan failure: ground state off reference"

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
  if (allocated(error)) error stop "T_Ne3d_10_TimePropagationMctdhxLaserVelocityBlanesMoan failure: norm not conserved"

  ! The dipole moment must respond to the field (regression value). It is not
  ! comparable to T_Ne3d_07: for the 1-cycle pulse the envelope-derivative
  ! part of E = -dA/dt differs substantially from the length-gauge field.
  ! Difference to the SIL value of T_Ne3d_08 (1.6e-7 here) is the frozen-A(t)
  ! error of that SIL step at omega = 0.5; it drops fourfold per halved step
  call check(error, dipole, 1.3142887110e-3_R64, thr=3e-7_R64)
  if (allocated(error)) error stop "T_Ne3d_10_TimePropagationMctdhxLaserVelocityBlanesMoan failure: dipole off SIL value"

  call check(error, dipole, 1.3144526364e-3_R64, thr=1e-8_R64)
  if (allocated(error)) error stop "T_Ne3d_10_TimePropagationMctdhxLaserVelocityBlanesMoan failure: dipole off reference"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  include 'Ne3dMctdhxLaser.inc'

end program
