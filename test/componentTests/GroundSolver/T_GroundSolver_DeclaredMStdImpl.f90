!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> A standard ground solver together with declared orbital m
!> (method.mb.orbBased.mOfOrbs): the standard solvers diagonalize on the full
!> grid and would mix degenerate orbitals of different m, so
!> GroundSolver_Fabricate must stop with an error. The test passes only if
!> the program stops with that error:
!> ctest-pass-regex: mOfOrbs requires groundSolver.scf.ylmOpt or groundSolver.mcscf.ylmOpt
program T_GroundSolver_DeclaredMStdImpl
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Grid
  use M_Method
  use M_GroundSolver

  implicit none

  character(len=:), allocatable :: jsonFileName

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  jsonFileName = "test/componentTests/GroundSolver/T_GroundSolver_DeclaredMStdImpl.json"
  call Json_LoadJsonFile(manualFileName_=jsonFileName, relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call Method_Fabricate

  !==========================================
  call Say_Section("fabricate the ground solver (must stop with an error)")
  !==========================================

  call GroundSolver_Fabricate

  error stop "T_GroundSolver_DeclaredMStdImpl failure: groundSolver.scf.stdImpl was accepted with mOfOrbs"

end program T_GroundSolver_DeclaredMStdImpl
