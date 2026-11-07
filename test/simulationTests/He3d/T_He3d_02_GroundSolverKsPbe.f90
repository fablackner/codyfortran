!> ctest-requires: libxc
!> Kohn-Sham PBE ground state of He on a Ylm FEDVR grid (GGA: radial and
!> angular density gradients of sysXc.ylm): E = -2.892935 (the PBE value of
!> numerical atomic codes); the same functional on grid.fem gives the same
!> energy within the mesh error of the FEM nucleus (T_Fem_12_HeKsPbe).
program T_He3d_02_GroundSolverKsPbe
  use M_Utils_Types
  use M_Utils_Say
  use testdrive, only: check, error_type

  implicit none

  type(error_type), allocatable :: error
  real(R64) :: energy
  real(R64), allocatable :: eps(:)
  integer(I32) :: nScf

  call Say_Hello
  call RunKsScfYlm(energy, eps, nScf, "test/simulationTests/He3d/T_He3d_02_GroundSolverKsPbe.json", 40)

  print *
  print '(1X, A, F16.10, A, F14.8, A, I0)', "He PBE: E = ", energy, ", eps_1s = ", eps(1), ", SCF steps ", nScf

  call check(error, energy, -2.892935_R64, thr=2e-6_R64)
  if (allocated(error)) error stop "T_He3d_02_GroundSolverKsPbe failure: energy (PBE -2.892935)"

  call Say_Goodbye

contains

  include "../KsScfYlm.inc"

end program
