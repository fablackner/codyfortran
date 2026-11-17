!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Ne on a coarse Ylm grid (lmax 1, rmax 5, FEDVR 2 x 12): Hartree-Fock ground
!> state with groundSolver.scf.stdImpl, linear mixing of the density
!> (alpha 0.95), the energy -127.7449896 of T_Ne3d_02 (Ne3dGroundSolverScf.inc).
program T_Ne3d_05_GroundSolverScfDensity
  use M_Utils_Types
  use M_Utils_Say

  implicit none

  character(len=*), parameter :: jsonFileName = &
                                 "test/simulationTests/Ne3d/T_Ne3d_05_GroundSolverScfDensity.json"

  call Say_Hello
  call CheckScfGroundState(jsonFileName, -127.7449896_R64)
  call Say_Goodbye

contains

  include 'Ne3dGroundSolverScf.inc'

end program
