!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> T_Ne3d_11_TimePropagationMctdhxLaserVelocitySpectralCut with restricted
!> orbitals (orbs.restrictedQ) and the declared m of the orbitals
!> (method.mb.orbBased.mOfOrbs: 1s 2s 2p(-1) 2p(0) 2p(1), equal for both
!> spins), as in the Ne production runs: the node path then works on the
!> stored spatial orbitals (mOfOrbs(1:nOrbsInState), Orbs_blockOfOrb by m
!> only) and duplicates h2 for the second spin.
!>
!> The m of every orbital is conserved by all operators of a z-polarized
!> field; the only projection on it is done once at initialization, and the
!> eigen-decompositions in orbital space (alignment of the ground solver,
!> inverse 1-RDM) act per m block, so the channels of other m must stay
!> exactly zero.
!>
!> Checks:
!> - the channels of other m are exactly zero after every output step
!> - the dipole moment agrees with the unrestricted T_Ne3d_11 (closed shell)
!> Writes the dipole history dipoleDeclaredM.dat, which
!> T_Ne3d_14_TimePropagationMctdhxUndeclaredM compares with the generic path.
program T_Ne3d_13_TimePropagationMctdhxDeclaredM
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Orbs
  use M_Method
  use M_Method_Mb_OrbBased
  use testdrive, only: check, error_type

  implicit none

  real(R64)    :: time, energy, dipole
  integer(I32) :: nOtherM, io
  type(error_type), allocatable :: error

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call SetupMctdhxLaser("test/simulationTests/Ne3d/T_Ne3d_13_TimePropagationMctdhxDeclaredM.json")

  !==========================================
  call Say_Section("converge mcscf ground state")
  !==========================================

  energy = ConvergedGroundStateEnergy()

  ! Gauge coupling is off at t = 0: must reproduce the T_Ne3d_06/07 mcscf energy.
  call check(error, energy, -127.7449896_R64, thr=1e-7_R64)
  if (allocated(error)) error stop "T_Ne3d_13_TimePropagationMctdhxDeclaredM failure: ground state off reference"

  !==========================================
  call Say_Section("propagate through laser pulse")
  !==========================================

  nOtherM = OtherMEntries(Orbs_orbs)

  open (newunit=io, file="dipoleDeclaredM.dat", status="replace")

  time = Json_Get("program.startTime", 0.0_R64)
  do while (BeforeEndTimeQ(time))
    call PropagateOutputStep(time, dipole, "dipoleVelocity.dat")
    nOtherM = nOtherM + OtherMEntries(Orbs_orbs)
    write (io, '(2ES26.17)') time, dipole
  end do

  close (io)

  !==========================================
  call Say_Section("printout results")
  !==========================================

  print *
  print *, "final dipole: ", dipole
  print *, "nonzero entries in channels of other m: ", nOtherM

  call check(error, Method_Mb_OrbBased_definiteMQ)
  if (allocated(error)) error stop "T_Ne3d_13_TimePropagationMctdhxDeclaredM failure: mOfOrbs not declared"

  ! Exact zeros, not round-off
  call check(error, nOtherM, 0)
  if (allocated(error)) error stop "T_Ne3d_13_TimePropagationMctdhxDeclaredM failure: orbitals left their m"

  call check(error, dipole, 1.3143644064e-3_R64, thr=1e-11_R64)
  if (allocated(error)) error stop "T_Ne3d_13_TimePropagationMctdhxDeclaredM failure: dipole off the unrestricted T_Ne3d_11"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  include 'Ne3dMctdhxLaser.inc'

end program
