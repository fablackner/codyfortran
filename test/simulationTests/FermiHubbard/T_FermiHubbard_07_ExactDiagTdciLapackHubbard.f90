!> Fermi-Hubbard lattice with coeffs.hubbard.noSpinSym: ground state of the
!> TDCI Hamiltonian by LAPACK (HubbardExactDiag.inc), the energy -7.29797344
!> of the coeffs.generic tests T_FermiHubbard_01/02.
program T_FermiHubbard_07_ExactDiagTdciLapackHubbard
  use M_Utils_Types
  use M_Utils_Say

  implicit none

  character(len=*), parameter :: jsonFileName = &
                                 "test/simulationTests/FermiHubbard/T_FermiHubbard_07_ExactDiagTdciLapackHubbard.json"

  call Say_Hello
  call CheckExactDiagGroundState(jsonFileName, -7.29797344_R64)
  call Say_Goodbye

contains

  include 'HubbardExactDiag.inc'

end program
