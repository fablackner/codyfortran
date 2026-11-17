!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Spin-restricted MCTDHX imaginary-time propagation for 1D helium.
!>
!> Identical physics to T_He1d_04_ImagTimePropagationMctdhx, but with
!> orbs.restrictedQ = true: only the shared spatial orbital set is stored and
!> propagated, and the spin blocks of h1/h2 are duplicated internally. The
!> converged ground-state energy must therefore match the unrestricted
!> reference value.
program T_He1d_07_ImagTimePropagationMctdhxRestricted
  use M_Utils_Types
  use M_Utils_Say

  implicit none

  character(len=*), parameter :: jsonFileName = &
                                 "test/simulationTests/He1d/T_He1d_07_ImagTimePropagationMctdhxRestricted.json"

  call Say_Hello

  call CheckImagTimeGroundState(jsonFileName, -2.4852442075051160_R64, 1e-7_R64, normalizeCoeffsQ_=.true.)
  call Say_Goodbye

contains

  include '../ImagTimeGroundState.inc'

end program
