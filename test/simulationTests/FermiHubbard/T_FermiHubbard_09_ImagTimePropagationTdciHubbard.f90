!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Fermi-Hubbard 2x2 periodic lattice (2 + 2 fermions, U = 1) with
!> coeffs.hubbard.noSpinSym: imaginary-time TDCI to the exact ground-state
!> energy -7.29797344 of T_FermiHubbard_01/02 (ImagTimeGroundState.inc).
program T_FermiHubbard_09_ImagTimePropagationTdciHubbard
  use M_Utils_Types
  use M_Utils_Say

  implicit none

  character(len=*), parameter :: jsonFileName = &
                                 "test/simulationTests/FermiHubbard/T_FermiHubbard_09_ImagTimePropagationTdciHubbard.json"

  call Say_Hello

  call CheckImagTimeGroundState(jsonFileName, -7.29797344_R64, 1e-7_R64, normalizeCoeffsQ_=.true.)

  call Say_Goodbye

contains

  include '../ImagTimeGroundState.inc'

end program
