!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Fermi-Hubbard 2x2 periodic lattice (2 + 2 fermions, U = 1): imaginary-time
!> MCTDHF with all 4 orbitals per spin, the exact ground-state energy
!> -7.29797344 of T_FermiHubbard_01/02.
program T_FermiHubbard_05_ImagTimePropagationMctdhx
  use M_Utils_Types
  use M_Utils_Say

  implicit none

  character(len=*), parameter :: jsonFileName = &
                                 "test/simulationTests/FermiHubbard/T_FermiHubbard_05_ImagTimePropagationMctdhx.json"

  call Say_Hello

  call CheckImagTimeGroundState(jsonFileName, -7.29797344_R64, 1e-7_R64, normalizeCoeffsQ_=.true.)
  call Say_Goodbye

contains

  include '../ImagTimeGroundState.inc'

end program
