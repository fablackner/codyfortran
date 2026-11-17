!> Bose-Hubbard lattice: ground state of the TDCI Hamiltonian by ARPACK
!> (HubbardExactDiag.inc); the same reference energy as the LAPACK test.
program T_BoseHubbard_02_ExactDiagTdciArpack
  use M_Utils_Types
  use M_Utils_Say

  implicit none

  character(len=*), parameter :: jsonFileName = "test/simulationTests/BoseHubbard/T_BoseHubbard_02_ExactDiagTdciArpack.json"

  call Say_Hello
  call CheckExactDiagGroundState(jsonFileName, -14.60164824_R64)
  call Say_Goodbye

contains

  include '../FermiHubbard/HubbardExactDiag.inc'

end program
