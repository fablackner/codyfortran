!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Bose-Hubbard 2x2 periodic lattice (4 bosons, U = 1): imaginary-time TDHX
!> to the mean-field energy -14.5.
program T_BoseHubbard_05_ImagTimePropagationTdhx
  use M_Utils_Types
  use M_Utils_Say

  implicit none

  character(len=*), parameter :: jsonFileName = &
                                 "test/simulationTests/BoseHubbard/T_BoseHubbard_05_ImagTimePropagationTdhx.json"

  call Say_Hello

  call CheckImagTimeGroundState(jsonFileName, -14.5_R64, 1e-7_R64)
  call Say_Goodbye

contains

  include '../ImagTimeGroundState.inc'

end program
