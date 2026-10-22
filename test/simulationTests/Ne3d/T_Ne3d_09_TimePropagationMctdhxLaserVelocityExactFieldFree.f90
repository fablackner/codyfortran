!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> T_Ne3d_08_TimePropagationMctdhxLaserVelocity with exact field-free one-body
!> propagation (method.mb.orbBased.exactFieldFreeQ).
!>
!> Every linear split step is exp(-i h0 dt/2) SIL[gauge](dt) exp(-i h0 dt/2):
!> the field-free h0 = T + V_ext is diagonalized once per l and applied
!> exactly, the Short Iterative Lanczos only propagates the velocity-gauge
!> coupling A(t)·p + A²/2. The grid has lmax = 1 with all m and two body
!> types, which exercises the m-independence check of the h0 blocks.
!>
!> Checks:
!> - the MCSCF ground-state energy matches T_Ne3d_06/07/08 (gauge off at t = 0)
!> - the weighted-metric norm is conserved (exp(-i h0 dt) is unitary)
!> - the dipole moment agrees with T_Ne3d_08 within the O(dt^2) splitting
!>   error of h0 and the gauge coupling
program T_Ne3d_09_TimePropagationMctdhxLaserVelocityExactFieldFree
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

  call SetupMctdhxLaser("test/simulationTests/Ne3d/T_Ne3d_09_TimePropagationMctdhxLaserVelocityExactFieldFree.json")

  !==========================================
  call Say_Section("converge mcscf ground state")
  !==========================================

  energy = ConvergedGroundStateEnergy()

  ! Gauge coupling is off at t = 0: must reproduce the T_Ne3d_06/07 mcscf energy.
  call check(error, energy, -127.7449896_R64, thr=1e-7_R64)
  if (allocated(error)) error stop "T_Ne3d_09_TimePropagationMctdhxLaserVelocityExactFieldFree failure: ground state off reference"

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
  if (allocated(error)) error stop "T_Ne3d_09_TimePropagationMctdhxLaserVelocityExactFieldFree failure: norm not conserved"

  ! The dipole moment must respond to the field (regression value). It is not
  ! comparable to T_Ne3d_07: for the 1-cycle pulse the envelope-derivative
  ! part of E = -dA/dt differs substantially from the length-gauge field.
  ! The splitting of h0 and the gauge coupling adds an O(dt^2) error: the
  ! difference to the SIL value of T_Ne3d_08 (7.6e-8 here) drops fourfold per
  ! halved time step, and both schemes extrapolate to 1.314025e-3
  call check(error, dipole, 1.3142887110e-3_R64, thr=1e-7_R64)
  if (allocated(error)) error stop "T_Ne3d_09_TimePropagationMctdhxLaserVelocityExactFieldFree failure: dipole off SIL value"

  call check(error, dipole, 1.3143642081e-3_R64, thr=1e-8_R64)
  if (allocated(error)) error stop "T_Ne3d_09_TimePropagationMctdhxLaserVelocityExactFieldFree failure: dipole off reference"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  include 'Ne3dMctdhxLaser.inc'

end program
