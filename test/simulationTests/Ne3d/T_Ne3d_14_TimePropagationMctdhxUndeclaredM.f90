!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> T_Ne3d_13_TimePropagationMctdhxDeclaredM (restricted orbitals) without
!> method.mb.orbBased.mOfOrbs: the generic correlation operator and h2 (Gaunt
!> products), the exact field-free part on all (orbital, m) columns, no
!> projection and full eigen-decompositions in orbital space.
!>
!> Checks:
!> - the dipole moment agrees with the declared run (dipoleDeclaredM.dat of
!>   T_Ne3d_13) to round-off after every output step
program T_Ne3d_14_TimePropagationMctdhxUndeclaredM
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Orbs
  use M_Method
  use M_Method_Mb_OrbBased
  use testdrive, only: check, error_type

  implicit none

  real(R64)    :: time, energy, dipole
  real(R64)    :: maxDev, timeRef, dipoleRef
  integer(I32) :: io, ioStat
  type(error_type), allocatable :: error

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call SetupMctdhxLaser("test/simulationTests/Ne3d/T_Ne3d_14_TimePropagationMctdhxUndeclaredM.json")

  !==========================================
  call Say_Section("converge mcscf ground state")
  !==========================================

  energy = ConvergedGroundStateEnergy()

  ! Gauge coupling is off at t = 0: must reproduce the T_Ne3d_06/07 mcscf energy.
  call check(error, energy, -127.7449896_R64, thr=1e-7_R64)
  if (allocated(error)) error stop "T_Ne3d_14_TimePropagationMctdhxUndeclaredM failure: ground state off reference"

  !==========================================
  call Say_Section("propagate through laser pulse")
  !==========================================

  maxDev = 0.0_R64

  open (newunit=io, file="dipoleDeclaredM.dat", status="old", action="read", iostat=ioStat)
  if (ioStat .ne. 0) error stop "T_Ne3d_14_TimePropagationMctdhxUndeclaredM: run T_Ne3d_13 first (dipoleDeclaredM.dat)"

  time = Json_Get("program.startTime", 0.0_R64)
  do while (BeforeEndTimeQ(time))
    call PropagateOutputStep(time, dipole, "dipoleVelocity.dat")
    read (io, *) timeRef, dipoleRef
    if (abs(timeRef - time) > 1e-12_R64) &
      error stop "T_Ne3d_14_TimePropagationMctdhxUndeclaredM: time grid differs from T_Ne3d_13"
    maxDev = max(maxDev, abs(dipole - dipoleRef))
  end do

  close (io)

  !==========================================
  call Say_Section("printout results")
  !==========================================

  print *
  print *, "final dipole: ", dipole
  print *, "max. deviation from the declared run: ", maxDev

  call check(error, .not. Method_Mb_OrbBased_definiteMQ)
  if (allocated(error)) error stop "T_Ne3d_14_TimePropagationMctdhxUndeclaredM failure: mOfOrbs declared"

  ! Same physics, different operators and eigen-decompositions: round-off
  call check(error, maxDev, 0.0_R64, thr=1e-13_R64)
  if (allocated(error)) error stop "T_Ne3d_14_TimePropagationMctdhxUndeclaredM failure: dipole off the declared run"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  include 'Ne3dMctdhxLaser.inc'

end program
