!> Fermi-Hubbard lattice: ground state of the TDCI Hamiltonian by LAPACK
!> (HubbardExactDiag.inc); the same reference energy as the ARPACK test.
program T_FermiHubbard_01_ExactDiagTdciLapack
  use M_Utils_Types
  use M_Utils_Say

  implicit none

  character(len=*), parameter :: jsonFileName = "test/simulationTests/FermiHubbard/T_FermiHubbard_01_ExactDiagTdciLapack.json"

  call Say_Hello
  call CheckExactDiagGroundState(jsonFileName, -7.29797344_R64)
  call Say_Goodbye

contains

  include 'HubbardExactDiag.inc'

end program
