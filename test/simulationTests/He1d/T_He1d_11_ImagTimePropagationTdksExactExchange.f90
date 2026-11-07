!> ctest-requires: libxc
!> TDKS with exact exchange only (sysXc.linear.libxc: no functional,
!> exactExchangeFraction = 1) is time-dependent Hartree-Fock: imaginary-time
!> propagation of the restricted 1D helium of T_He1d_08 (same grid, operators
!> and orbitals) must reach the TDHX energy -2.4819265615333204. Checks the
!> Kohn-Sham operator (Hartree potential of the xc density on grid.linear,
!> exchange pair potentials) and energy against the independent TDHX path.
program T_He1d_11_ImagTimePropagationTdksExactExchange
  use M_Utils_Types
  use M_Utils_Say
  use testdrive, only: check, error_type

  implicit none

  type(error_type), allocatable :: error
  real(R64) :: energy

  call Say_Hello
  call RunTdksImagTime(energy, "test/simulationTests/He1d/T_He1d_11_ImagTimePropagationTdksExactExchange.json")

  print *
  print *, "TDKS (exact exchange) energy: ", energy, ", TDHX ", -2.4819265615333204_R64, &
    ", difference ", energy + 2.4819265615333204_R64

  call check(error, energy, -2.4819265615333204_R64, thr=1e-9_R64)
  if (allocated(error)) error stop "T_He1d_11_ImagTimePropagationTdksExactExchange failure: energy differs from TDHX"

  call Say_Goodbye

contains

  include "He1dTdksImagTime.inc"

end program
