!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Fermi-Hubbard 2x2 periodic lattice (2 + 2 fermions, U = 1): imaginary-time
!> TDHF to the Hartree-Fock energy -7.0.
program T_FermiHubbard_06_ImagTimePropagationTdhx
  use M_Utils_Types
  use M_Utils_Say

  implicit none

  character(len=*), parameter :: jsonFileName = &
                                 "test/simulationTests/FermiHubbard/T_FermiHubbard_06_ImagTimePropagationTdhx.json"

  call Say_Hello

  call CheckImagTimeGroundState(jsonFileName, -7.0_R64, 1e-7_R64)
  call Say_Goodbye

contains

  include '../ImagTimeGroundState.inc'

end program
