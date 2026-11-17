!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> H2+ 2p pi_u state via imaginary time propagation on the prolate-spheroidal
!> grid (odd-m channel, axis-adapted xi and eta basis).
!>
!> One electron in the field of two protons at fixed internuclear distance
!> R = 2.0, initialized from the ungerade LCAO 2p (m = +1) guess. Imaginary
!> time keeps m and the inversion parity, so the propagation converges to the
!> lowest pi_u state. The converged electronic energy is compared against the
!> exact two-center value:
!>
!>   E_el(2p pi_u, R = 2) = -0.4287718199
!>
!> On the plain Lagrange basis (without the sqrt(xi^2-1) axis factor) the error
!> did not decrease with the xi grid (0.06 to 0.6 mHa, partly below the exact value).
program T_H2_03_ImagTimePropagationSingleBodyPi
  use M_Utils_Types
  use M_Utils_Say

  implicit none

  character(len=*), parameter :: jsonFileName = &
                                 "test/simulationTests/H2/T_H2_03_ImagTimePropagationSingleBodyPi.json"

  call Say_Hello

  ! exact electronic 2p pi_u energy of H2+ at R = 2.0 (the converged
  ! grid reproduces it to ~4e-12)
  call CheckImagTimeSingleBody(jsonFileName, -0.42877182_R64, 1e-8_R64)

  call Say_Goodbye

contains

  include '../ImagTimeSingleBody.inc'

end program
