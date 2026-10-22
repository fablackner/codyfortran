!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> T_Ne3d_09_TimePropagationMctdhxLaserVelocityExactFieldFree with a spectral
!> cut of h0 (method.mb.orbBased.exactFieldFreeCut = 300, keeps 20 of the 22
!> eigenvectors per l, eigenvalues up to 570).
!>
!> The exact part propagates only the eigenvectors below the cut and the
!> linear and nonlinear orbital derivatives are projected onto them, so the
!> orbitals stay in the kept space.
!>
!> The m of the orbitals is declared (method.mb.orbBased.mOfOrbs: 1s 2s
!> 2p(-1) 2p(0) 2p(1) for both spins, unrestricted), which binds the
!> definite-m node quadrature of the correlation operator and h2 and the
!> m-selective exact field-free part. Every operator of the z-polarized field
!> conserves m, the only projection on it is done once at initialization, and
!> the eigen-decompositions in orbital space (alignment of the ground solver,
!> inverse 1-RDM) act per m block, so the channels of other m must stay
!> exactly zero.
!>
!> Checks:
!> - the MCSCF ground-state energy matches T_Ne3d_06/07/08 (gauge off at t = 0)
!> - the weighted-metric norm is conserved (the ground state has negligible
!>   weight above the cut)
!> - the propagated orbitals have no weight above the cut
!> - the dipole moment agrees with T_Ne3d_09 within the (small) effect of the
!>   cut on this grid
!> - the channels of other m are exactly zero after every output step
program T_Ne3d_11_TimePropagationMctdhxLaserVelocitySpectralCut
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
  real(R64)    :: cutResidual
  integer(I32) :: nOtherM
  complex(R64), allocatable :: projectedOrbs(:, :)
  type(error_type), allocatable :: error

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call SetupMctdhxLaser("test/simulationTests/Ne3d/T_Ne3d_11_TimePropagationMctdhxLaserVelocitySpectralCut.json")

  !==========================================
  call Say_Section("converge mcscf ground state")
  !==========================================

  energy = ConvergedGroundStateEnergy()

  ! Gauge coupling is off at t = 0: must reproduce the T_Ne3d_06/07 mcscf energy.
  call check(error, energy, -127.7449896_R64, thr=1e-7_R64)
  if (allocated(error)) error stop "T_Ne3d_11_TimePropagationMctdhxLaserVelocitySpectralCut failure: ground state off reference"

  !==========================================
  call Say_Section("propagate through laser pulse")
  !==========================================

  startNorm = MetricNorm(Method_state)
  nOtherM = OtherMEntries(Orbs_orbs)

  time = Json_Get("program.startTime", 0.0_R64)
  do while (BeforeEndTimeQ(time))
    call PropagateOutputStep(time, dipole, "dipoleVelocity.dat")
    nOtherM = nOtherM + OtherMEntries(Orbs_orbs)
  end do

  !==========================================
  call Say_Section("printout results")
  !==========================================

  finalNorm = MetricNorm(Method_state)

  ! weight of the orbitals above the cut: ||orbs - P orbs|| (unit weights
  ! are fine for a zero check)
  projectedOrbs = Orbs_orbs
  call Method_Mb_OrbBased_ApplySpectralCut(projectedOrbs)
  cutResidual = sqrt(sum(abs(Orbs_orbs - projectedOrbs)**2))

  print *
  print *, "start norm:   ", startNorm
  print *, "final norm:   ", finalNorm
  print *, "final dipole: ", dipole
  print *, "cut residual: ", cutResidual
  print *, "nonzero entries in channels of other m: ", nOtherM

  ! Unitary propagation must conserve the weighted-metric norm of the state
  ! (up to the RK4 integration error of the nonlinear part).
  call check(error, finalNorm, startNorm, thr=1e-7_R64)
  if (allocated(error)) error stop "T_Ne3d_11_TimePropagationMctdhxLaserVelocitySpectralCut failure: norm not conserved"

  ! The cut space is invariant under the propagation
  call check(error, cutResidual, 0.0_R64, thr=1e-12_R64)
  if (allocated(error)) error stop "T_Ne3d_11_TimePropagationMctdhxLaserVelocitySpectralCut failure: weight above the cut"

  ! The cut eigenvectors (E > 300) are barely populated by this field: the
  ! dipole stays within 2e-10 of T_Ne3d_09 (without cut)
  call check(error, dipole, 1.3143642081e-3_R64, thr=1e-9_R64)
  if (allocated(error)) error stop "T_Ne3d_11_TimePropagationMctdhxLaserVelocitySpectralCut failure: dipole off T_Ne3d_09"

  call check(error, dipole, 1.3143644064e-3_R64, thr=1e-11_R64)
  if (allocated(error)) error stop "T_Ne3d_11_TimePropagationMctdhxLaserVelocitySpectralCut failure: dipole off reference"

  call check(error, Method_Mb_OrbBased_definiteMQ)
  if (allocated(error)) error stop "T_Ne3d_11_TimePropagationMctdhxLaserVelocitySpectralCut failure: mOfOrbs not declared"

  ! Exact zeros, not round-off
  call check(error, nOtherM, 0)
  if (allocated(error)) error stop "T_Ne3d_11_TimePropagationMctdhxLaserVelocitySpectralCut failure: orbitals left their m"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  include 'Ne3dMctdhxLaser.inc'

end program
