!> Bose-Hubbard lattice: ground state of the TDCI Hamiltonian by LAPACK
!> (HubbardExactDiag.inc); the same reference energy as the ARPACK test.
program T_BoseHubbard_01_ExactDiagTdciLapack
  use M_Utils_Types
  use M_Utils_Say

  implicit none

  character(len=*), parameter :: jsonFileName = "test/simulationTests/BoseHubbard/T_BoseHubbard_01_ExactDiagTdciLapack.json"

  call Say_Hello
  call CheckExactDiagGroundState(jsonFileName, -14.60164824_R64)
  call Say_Goodbye

contains

  include '../FermiHubbard/HubbardExactDiag.inc'

end program
