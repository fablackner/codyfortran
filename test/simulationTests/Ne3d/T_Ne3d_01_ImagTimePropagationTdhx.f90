!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Ne on a coarse Ylm grid (lmax 1, rmax 5, two FEDVR elements): imaginary-time
!> TDHF from hydrogen-like orbitals to the Hartree-Fock energy of this grid,
!> -108.8479585.
program T_Ne3d_01_ImagTimePropagationTdhx
  use M_Utils_Types
  use M_Utils_Say

  implicit none

  character(len=*), parameter :: jsonFileName = &
                                 "test/simulationTests/Ne3d/T_Ne3d_01_ImagTimePropagationTdhx.json"

  call Say_Hello

  call CheckImagTimeGroundState(jsonFileName, -108.8479585_R64, 1e-6_R64)
  call Say_Goodbye

contains

  include '../ImagTimeGroundState.inc'

end program
