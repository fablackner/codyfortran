!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Bose-Hubbard 2x2 periodic lattice (4 bosons, U = 1): imaginary-time MCTDHB
!> with all 4 orbitals, the exact ground-state energy -14.60164824 of
!> T_BoseHubbard_01/02.
program T_BoseHubbard_04_ImagTimePropagationMctdhx
  use M_Utils_Types
  use M_Utils_Say

  implicit none

  character(len=*), parameter :: jsonFileName = &
                                 "test/simulationTests/BoseHubbard/T_BoseHubbard_04_ImagTimePropagationMctdhx.json"

  call Say_Hello

  call CheckImagTimeGroundState(jsonFileName, -14.60164824_R64, 1e-7_R64, normalizeCoeffsQ_=.true.)
  call Say_Goodbye

contains

  include '../ImagTimeGroundState.inc'

end program
