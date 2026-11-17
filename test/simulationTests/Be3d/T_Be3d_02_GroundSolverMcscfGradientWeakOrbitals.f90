!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> MCSCF ground state of Be with groundSolver.mcscf.gradient and weakly
!> occupied orbitals: 9 restricted orbitals 1s 2s 2p 3s 3p of declared m, full
!> CI in the total-M = 0 sector (ARPACK), Ylm grid lmax = mmax = 1 (rmax 40,
!> FEDVR 20 x 15), start from hydrogen-like orbitals (Z = 3) and a unary CI
!> vector, plain linear mixing with alpha = 1 (no DIIS), so that the bare
!> iteration has to be stable.
!>
!> The natural occupations per spin are 0.99844 (1s), 0.90366 (2s), 0.032112
!> (2p), 6.25e-4 (3s) and 3.15e-4 (3p): 2s/3s and 2p/3p share a symmetry
!> block, so the stored orbitals are mixtures of natural orbitals of very
!> different occupation, and the 3s 3p feel a strong CI response. Without the
!> natural-orbital frame of the step the iteration stalls (E = -14.374 after
!> 150 iterations); with the CI vector of the previous orbitals in the
!> residual (one pass per iteration) it diverges in a period-2 oscillation
!> (with the natural-orbital frame) or creeps (without: gradient 7e-4 after
!> 150 iterations).
!>
!> E = -14.6539148604 on this grid (CASSCF(4,9) in the uncontracted s,p
!> primitives of cc-pV5Z + cc-pCVQZ: -14.6538931 with the same occupations
!> to 4 digits; that basis misses 1.5e-5 already in Hartree-Fock). Also
!> checks the orbital energy gradient at convergence and the iteration count
!> (26 at the time of writing).
program T_Be3d_02_GroundSolverMcscfGradientWeakOrbitals
  use M_Utils_Types
  use M_Utils_Say

  implicit none

  character(len=*), parameter :: jsonFileName = &
                                 "test/simulationTests/Be3d/T_Be3d_02_GroundSolverMcscfGradientWeakOrbitals.json"

  call Say_Hello
  call CheckMcscfGradientGroundState(jsonFileName, -14.6539148604_R64, 40)
  call Say_Goodbye

contains

  include 'Be3dMcscfGradient.inc'

end program
