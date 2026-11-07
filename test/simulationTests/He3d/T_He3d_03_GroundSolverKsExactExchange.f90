!> ctest-requires: libxc
!> Kohn-Sham machinery with exact exchange only (no Libxc functional,
!> sysXc.ylm.libxc.exactExchangeFraction = 1, ACE in the SCF): the Kohn-Sham
!> energy Σ<g|hφ> + E_H - E_x is the Hartree-Fock energy, so the ground state
!> is the Hartree-Fock limit of He, -2.8616800 (Clementi-Roetti), and the
!> orbital energy -0.9179556.
program T_He3d_03_GroundSolverKsExactExchange
  use M_Utils_Types
  use M_Utils_Say
  use testdrive, only: check, error_type

  implicit none

  type(error_type), allocatable :: error
  real(R64) :: energy
  real(R64), allocatable :: eps(:)
  integer(I32) :: nScf

  call Say_Hello
  call RunKsScfYlm(energy, eps, nScf, "test/simulationTests/He3d/T_He3d_03_GroundSolverKsExactExchange.json", 40)

  print *
  print '(1X, A, F16.10, A, F14.8, A, I0)', "He exact exchange: E = ", energy, ", eps_1s = ", eps(1), &
    ", SCF steps ", nScf

  call check(error, energy, -2.8616800_R64, thr=2e-7_R64)
  if (allocated(error)) error stop "T_He3d_03_GroundSolverKsExactExchange failure: energy (HF limit -2.8616800)"
  call check(error, eps(1), -0.9179556_R64, thr=2e-7_R64)
  if (allocated(error)) error stop "T_He3d_03_GroundSolverKsExactExchange failure: eps_1s (HF -0.9179556)"

  call Say_Goodbye

contains

  include "../KsScfYlm.inc"

end program
