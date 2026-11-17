!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Hydrogen atom on a Ylm grid (lmax 1, rmax 10, FEDVR 10 x 11): imaginary-time
!> propagation of the single-body method from the hydrogen-like 1s orbital,
!> energy -0.5.
program T_H3d_01_ImagTimePropagationSingleBody
  use M_Utils_Types
  use M_Utils_Say

  implicit none

  character(len=*), parameter :: jsonFileName = &
                                 "test/simulationTests/H3d/T_H3d_01_ImagTimePropagationSingleBody.json"

  call Say_Hello

  call CheckImagTimeSingleBody(jsonFileName, -0.5_R64, 1e-5_R64)

  call Say_Goodbye

contains

  include '../ImagTimeSingleBody.inc'

end program
