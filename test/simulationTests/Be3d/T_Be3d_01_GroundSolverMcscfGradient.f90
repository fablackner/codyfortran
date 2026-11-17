!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> MCSCF ground state of Be with groundSolver.mcscf.gradient: 5 restricted
!> orbitals 1s 2s 2p of declared m, full CI in the total-M = 0 sector, Ylm
!> grid lmax = mmax = 1 (rmax 40, FEDVR 20 x 15), start from hydrogen-like
!> orbitals (Z = 3) and a unary CI vector.
!>
!> The valence minimum 1s²(2s² + 2p²) has E = -14.6168557502 on this grid
!> (2p occupation 0.0325 per orbital; CASSCF in an s,p-only uncontracted
!> cc-pV5Z basis: -14.616842, imaginary-time MCTDHF on this grid
!> -14.6168557170 at dtau 0.005, the time-step bias). The core-correlated
!> local minimum, where the 2p shell correlates 1s², has E = -14.5956: a
!> regression into it fails the energy check. Also checks the orbital
!> energy gradient at convergence and the iteration count (12 at the time of
!> writing).
!>
!> Writes the state (orbXX_YY.in, coeffs.in, rdm2BtXX_YY.in and the up-up-down 3-RDM
!> rdm3Bt01_01_02.in) for T_Be3d_03 and T_Be3d_04.
program T_Be3d_01_GroundSolverMcscfGradient
  use M_Utils_Types
  use M_Utils_Say
  use M_Orbs
  use M_Coeffs
  use M_Reconstruction_Singlet_Anchor

  implicit none

  character(len=*), parameter :: jsonFileName = &
                                 "test/simulationTests/Be3d/T_Be3d_01_GroundSolverMcscfGradient.json"

  call Say_Hello
  call CheckMcscfGradientGroundState(jsonFileName, -14.6168557502_R64, 20)

  !==========================================
  call Say_Section("save data")
  !==========================================

  call Orbs_SaveOrbs(Orbs_orbs)
  call Coeffs_SaveCoeffs(Coeffs_coeffs)
  call Coeffs_SaveTwoRdm(Coeffs_coeffs)
  call Reconstruction_Singlet_Anchor_SaveReference(reshape(Coeffs_coeffs, [Coeffs_nCoeffs, 1]))

  call Say_Goodbye

contains

  include 'Be3dMcscfGradient.inc'

end program
