!> ctest-requires: libxc
!> 1D helium with soft-Coulomb interactions 1/sqrt(x^2 + 1) and the matching
!> 1D LDA of Libxc (lda_x_1d_soft with beta = 1, lda_c_1d_csc with the
!> soft-Coulomb interaction and beta = 1): imaginary-time TDKS reaches the
!> Kohn-Sham energy of the SCF of T_He1d_12 (two independent solvers of the
!> same equations). The imaginary-time fixed point of RK4 with the nonlinear
!> Kohn-Sham potential carries an O(dtau^2) bias: about 6e-5 at dtau = 0.05
!> (this test), 2.6e-6 at dtau = 0.01 and 6.4e-7 at dtau = 0.005, whose
!> Richardson extrapolation, -2.2013755803, equals the SCF energy to 1e-9.
program T_He1d_13_ImagTimePropagationTdksLda1d
  use M_Utils_Types
  use M_Utils_Say
  use testdrive, only: check, error_type

  implicit none

  real(R64), parameter :: E_SCF = -2.2014358365_R64

  type(error_type), allocatable :: error
  real(R64) :: energy

  call Say_Hello
  call RunTdksImagTime(energy, "test/simulationTests/He1d/T_He1d_13_ImagTimePropagationTdksLda1d.json")

  print *
  print *, "1D LDA imaginary-time TDKS energy: ", energy, ", SCF ", E_SCF, ", difference ", energy - E_SCF

  ! the RK4 stages of the unshifted imaginary-time derivative see rescaled
  ! densities: the fixed point is biased by about 2.5e-2 dt^2 (6e-5 at dt 0.05)
  call check(error, energy, E_SCF, thr=1e-4_R64)
  if (allocated(error)) error stop "T_He1d_13_ImagTimePropagationTdksLda1d failure: energy differs from the SCF"

  call Say_Goodbye

contains

  include "He1dTdksImagTime.inc"

end program
