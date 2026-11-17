!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> 1D helium: imaginary-time MCTDHF (4 orbitals per spin) from the
!> Hartree-Fock orbitals written by T_He1d_03, energy -2.4852442075051160.
!> Writes the MCTDHF ground state (orb*.in, coeffs.in and rdm2Bt*.in) that
!> T_He1d_05, the T_He1d_06 tests, T_He1d_09, T_He1d_10 and T_He1d_15 read.
program T_He1d_04_ImagTimePropagationMctdhx
  use M_Utils_Types
  use M_Utils_Say
  use M_Orbs
  use M_Coeffs

  implicit none

  character(len=*), parameter :: jsonFileName = &
                                 "test/simulationTests/He1d/T_He1d_04_ImagTimePropagationMctdhx.json"

  call Say_Hello

  call CheckImagTimeGroundState(jsonFileName, -2.4852442075051160_R64, 1e-7_R64, normalizeCoeffsQ_=.true.)

  !==========================================
  call Say_Section("save data")
  !==========================================

  call Coeffs_SaveTwoRdm(Coeffs_coeffs)
  call Coeffs_SaveCoeffs(Coeffs_coeffs)
  call Orbs_SaveOrbs(Orbs_orbs)

  call Say_Goodbye

contains

  include '../ImagTimeGroundState.inc'

end program
