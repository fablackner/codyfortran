!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Spin-restricted TDHX imaginary-time propagation for 1D helium.
!>
!> Identical physics to T_He1d_03_ImagTimePropagationTdhx, but with
!> orbs.restrictedQ = true: only the shared spatial orbital set is stored and
!> propagated, and the mean-field (Hartree) term uses the spin-summed density.
!> The converged Hartree-Fock energy must therefore match the unrestricted
!> reference value.
program T_He1d_08_ImagTimePropagationTdhxRestricted
  use M_Utils_Types
  use M_Utils_Say

  implicit none

  character(len=*), parameter :: jsonFileName = &
                                 "test/simulationTests/He1d/T_He1d_08_ImagTimePropagationTdhxRestricted.json"

  call Say_Hello

  call CheckImagTimeGroundState(jsonFileName, -2.4819265615333204_R64, 1e-7_R64)
  call Say_Goodbye

contains

  include '../ImagTimeGroundState.inc'

end program
