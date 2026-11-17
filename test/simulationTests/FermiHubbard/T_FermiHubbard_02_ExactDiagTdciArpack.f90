!> Fermi-Hubbard lattice: ground state of the TDCI Hamiltonian by ARPACK
!> (HubbardExactDiag.inc); the same reference energy as the LAPACK test.
program T_FermiHubbard_02_ExactDiagTdciArpack
  use M_Utils_Types
  use M_Utils_Say

  implicit none

  character(len=*), parameter :: jsonFileName = "test/simulationTests/FermiHubbard/T_FermiHubbard_02_ExactDiagTdciArpack.json"

  call Say_Hello
  call CheckExactDiagGroundState(jsonFileName, -7.29797344_R64)
  call Say_Goodbye

contains

  include 'HubbardExactDiag.inc'

end program
