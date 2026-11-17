!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> H2 Hartree-Fock ground state via imaginary time propagation (TDHX) on the
!> prolate-spheroidal grid.
!>
!> Two electrons (spin up/down as two fermionic body types with one orbital
!> each) in the field of two protons at the fixed equilibrium distance
!> R = 1.4, starting from the gerade LCAO 1s sigma_g guess. The converged
!> electronic energy is compared against the Hartree-Fock limit
!>
!>   E_HF(total, R = 1.4) = -1.13362957 (Kolos/Roothaan-quality HF limit)
!>   E_el = E_HF - 1/R    = -1.13362957 - 0.71428571 = -1.84791528
!>
!> (the constant nuclear repulsion 1/R is not part of the electronic
!> Hamiltonian and is excluded from the comparison).
program T_H2_02_ImagTimePropagationTdhx
  use M_Utils_Types
  use M_Utils_Say

  implicit none

  character(len=*), parameter :: jsonFileName = &
                                 "test/simulationTests/H2/T_H2_02_ImagTimePropagationTdhx.json"

  call Say_Hello

  ! Hartree-Fock limit of the electronic energy at R = 1.4 (see header);
  ! the converged grid value -1.84791526025 agrees with it to ~2e-8
  call CheckImagTimeGroundState(jsonFileName, -1.8479152602_R64, 1e-7_R64)
  call Say_Goodbye

contains

  include '../ImagTimeGroundState.inc'

end program
