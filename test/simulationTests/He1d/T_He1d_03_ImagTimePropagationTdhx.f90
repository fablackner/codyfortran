!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> 1D helium: imaginary-time TDHF (method.mb.orbBased.tdhx) from the
!> hydrogen-like start to the Hartree-Fock energy -2.4819265615333204.
!> Writes the Hartree-Fock orbitals (orb*.in) that T_He1d_04 starts from.
program T_He1d_03_ImagTimePropagationTdhx
  use M_Utils_Types
  use M_Utils_Say
  use M_Orbs

  implicit none

  character(len=*), parameter :: jsonFileName = &
                                 "test/simulationTests/He1d/T_He1d_03_ImagTimePropagationTdhx.json"

  call Say_Hello

  call CheckImagTimeGroundState(jsonFileName, -2.4819265615333204_R64, 1e-7_R64)

  !==========================================
  call Say_Section("save data")
  !==========================================

  call Orbs_SaveOrbs(Orbs_orbs)

  call Say_Goodbye

contains

  include '../ImagTimeGroundState.inc'

end program
