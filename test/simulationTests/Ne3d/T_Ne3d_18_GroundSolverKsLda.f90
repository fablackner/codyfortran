!> ctest-requires: libxc
!> Kohn-Sham LDA ground state of Ne (Slater exchange + VWN5 correlation) on a
!> Ylm FEDVR grid (lmax 1, mmax 1: sysXc.ylm in the general mode with nPhi = 5)
!> against the NIST LDA reference: E = -128.233481, eps = -30.305855 (1s),
!> -1.322809 (2s), -0.498034 (2p).
program T_Ne3d_18_GroundSolverKsLda
  use M_Utils_Types
  use M_Utils_Say
  use testdrive, only: check, error_type

  implicit none

  type(error_type), allocatable :: error
  real(R64) :: energy
  real(R64), allocatable :: eps(:)
  integer(I32) :: nScf

  call Say_Hello
  call RunKsScfYlm(energy, eps, nScf, "test/simulationTests/Ne3d/T_Ne3d_18_GroundSolverKsLda.json", 60)

  print *
  print '(1X, A, F18.10, A, *(F14.8))', "Ne LDA: E = ", energy, ", eps = ", eps
  print '(1X, A, I0)', "SCF steps ", nScf

  call check(error, energy, -128.233481_R64, thr=1e-6_R64)
  if (allocated(error)) error stop "T_Ne3d_18_GroundSolverKsLda failure: energy (NIST LDA -128.233481)"
  call check(error, eps(1), -30.305855_R64, thr=2e-6_R64)
  if (allocated(error)) error stop "T_Ne3d_18_GroundSolverKsLda failure: eps_1s (NIST LDA -30.305855)"
  call check(error, eps(2), -1.322809_R64, thr=2e-6_R64)
  if (allocated(error)) error stop "T_Ne3d_18_GroundSolverKsLda failure: eps_2s (NIST LDA -1.322809)"
  call check(error, maxval(abs(eps(3:5) + 0.498034_R64)) < 2e-6_R64)
  if (allocated(error)) error stop "T_Ne3d_18_GroundSolverKsLda failure: eps_2p (NIST LDA -0.498034)"

  call Say_Goodbye

contains

  include "../KsScfYlm.inc"

end program
