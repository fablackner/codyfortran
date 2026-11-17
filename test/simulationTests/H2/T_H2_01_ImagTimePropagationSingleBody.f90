!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> H2+ molecular ion ground state via imaginary time propagation on the
!> prolate-spheroidal grid.
!>
!> One electron in the field of two protons at fixed internuclear distance
!> R = 2.0 (foci of the prolate grid), initialized from the gerade LCAO
!> 1s sigma_g guess. The converged electronic energy is compared against the
!> exact two-center value (Bates et al. / standard tables):
!>
!>   E_el(1s sigma_g, R = 2) = -1.10263421449...
!>
!> (the constant nuclear repulsion 1/R = 0.5 is not part of the electronic
!> Hamiltonian and is excluded from the comparison).
program T_H2_01_ImagTimePropagationSingleBody
  use M_Utils_Types
  use M_Utils_Say

  implicit none

  character(len=*), parameter :: jsonFileName = &
                                 "test/simulationTests/H2/T_H2_01_ImagTimePropagationSingleBody.json"

  call Say_Hello

  ! exact electronic 1s sigma_g energy of H2+ at R = 2.0 (the converged
  ! grid reproduces it to ~4e-12)
  call CheckImagTimeSingleBody(jsonFileName, -1.10263421449_R64, 1e-8_R64)

  call Say_Goodbye

contains

  include '../ImagTimeSingleBody.inc'

end program
