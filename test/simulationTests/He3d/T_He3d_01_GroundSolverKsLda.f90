!> ctest-requires: libxc
!> Kohn-Sham LDA ground state of He (Slater exchange + VWN5 correlation, the
!> NIST LDA reference: E = -2.834836, eps_1s = -0.570425) on a Ylm FEDVR grid
!> by groundSolver.scf.stdImpl with the Kohn-Sham mean field of sysXc.ylm.
program T_He3d_01_GroundSolverKsLda
  use M_Utils_Types
  use M_Utils_Say
  use testdrive, only: check, error_type

  implicit none

  type(error_type), allocatable :: error
  real(R64) :: energy
  real(R64), allocatable :: eps(:)
  integer(I32) :: nScf

  call Say_Hello
  call RunKsScfYlm(energy, eps, nScf, "test/simulationTests/He3d/T_He3d_01_GroundSolverKsLda.json", 40)

  print *
  print '(1X, A, F16.10, A, F14.8, A, I0)', "He LDA: E = ", energy, ", eps_1s = ", eps(1), ", SCF steps ", nScf

  call check(error, energy, -2.834836_R64, thr=2e-6_R64)
  if (allocated(error)) error stop "T_He3d_01_GroundSolverKsLda failure: energy (NIST LDA -2.834836)"
  call check(error, eps(1), -0.570425_R64, thr=2e-6_R64)
  if (allocated(error)) error stop "T_He3d_01_GroundSolverKsLda failure: eps_1s (NIST LDA -0.570425)"

  call Say_Goodbye

contains

  include "../KsScfYlm.inc"

end program
