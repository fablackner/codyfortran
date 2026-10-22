!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> T_Ne3d_15_TimePropagationMctdhxSector without coeffs.generic.totalM: the
!> Hamiltonian acts on all 100 configurations, the configurations of total
!> M /= 0 keep round-off amplitudes only.
!>
!> Checks:
!> - the ground-state energy and the dipole moment agree with the sector run
!>   (dipoleSector.dat of T_Ne3d_15) to round-off after every output step
program T_Ne3d_16_TimePropagationMctdhxFullCi
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Method
  use testdrive, only: check, error_type

  implicit none

  real(R64)    :: time, energy, dipole
  real(R64)    :: maxDev, energyRef, timeRef, dipoleRef
  integer(I32) :: io, ioStat
  type(error_type), allocatable :: error

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call SetupMctdhxLaser("test/simulationTests/Ne3d/T_Ne3d_16_TimePropagationMctdhxFullCi.json")

  !==========================================
  call Say_Section("converge mcscf ground state")
  !==========================================

  energy = ConvergedGroundStateEnergy()

  open (newunit=io, file="dipoleSector.dat", status="old", action="read", iostat=ioStat)
  if (ioStat .ne. 0) error stop "T_Ne3d_16_TimePropagationMctdhxFullCi: run T_Ne3d_15 first (dipoleSector.dat)"
  read (io, *) energyRef

  !==========================================
  call Say_Section("propagate through laser pulse")
  !==========================================

  maxDev = 0.0_R64

  time = Json_Get("program.startTime", 0.0_R64)
  do while (BeforeEndTimeQ(time))
    call PropagateOutputStep(time, dipole, "dipoleVelocity.dat")
    read (io, *) timeRef, dipoleRef
    if (abs(timeRef - time) > 1e-12_R64) error stop "T_Ne3d_16_TimePropagationMctdhxFullCi: time grid differs from T_Ne3d_15"
    maxDev = max(maxDev, abs(dipole - dipoleRef))
  end do

  close (io)

  !==========================================
  call Say_Section("printout results")
  !==========================================

  print *
  print *, "ground-state energy deviation from the sector run: ", energy - energyRef
  print *, "final dipole: ", dipole
  print *, "max. dipole deviation from the sector run: ", maxDev

  call check(error, energy, energyRef, thr=1e-10_R64)
  if (allocated(error)) error stop "T_Ne3d_16_TimePropagationMctdhxFullCi failure: ground state off the sector run"

  call check(error, maxDev, 0.0_R64, thr=1e-13_R64)
  if (allocated(error)) error stop "T_Ne3d_16_TimePropagationMctdhxFullCi failure: dipole off the sector run"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  include 'Ne3dMctdhxLaser.inc'

end program
