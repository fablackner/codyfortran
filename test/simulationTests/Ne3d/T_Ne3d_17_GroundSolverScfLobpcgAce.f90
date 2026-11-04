!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Neon SCF ground state as in T_Ne3d_02_GroundSolverScf (same grid, same
!> scf.stdImpl solver with potential mixing, same reference energy), but with
!> - the preconditioned block eigensolver diagonalizerList.lobpcg, working in
!>   the grid metric (metricWeights = Grid_Ylm_weights) and warm-started from
!>   the current orbitals and then from the previous eigenvectors; the
!>   orbital preconditioner GroundSolver_PreconditionAction (bound by
!>   GroundSolver_Fabricate for every backend) is the identity here (the Ylm
!>   kinetic backend provides no SysKinetic_SolveShifted), so the eigensolves
!>   of this stiff FEDVR operator take O(100) iterations (maxIter is raised
!>   accordingly)
!> - the adaptively compressed exchange operator
!>   (groundSolver.scf.stdImpl.aceQ), exact on the occupied space, so the SCF
!>   fixed point and hence the energy are those of T_Ne3d_02.
!> Besides the energy it asserts that every LOBPCG solve converged (no
!> maxIter, no stagnation exit: the unpreconditioned solves have plateaus of
!> 10-30 iterations over which the residual drops by less than 10 % but
!> fluctuates by 40 % or more, which the stagnation test must tolerate) and
!> that the SCF converges within MAX_SCF_ITERATIONS steps (measured: 14).
program T_Ne3d_17_GroundSolverScfLobpcgAce
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Grid
  use M_Grid_Ylm, only: Grid_Ylm_weights
  use M_SysKinetic
  use M_SysPotential
  use M_SysInteraction
  use M_Orbs
  use M_OrbsInit
  use M_ConfigList
  use M_Coeffs
  use M_CoeffsInit
  use M_Method
  use M_DiagonalizerList
  use M_DiagonalizerList_Lobpcg, only: T_DiagonalizerList_E_Lobpcg
  use M_Mixing
  use M_GroundSolver
  use M_GroundSolver_Scf
  use testdrive, only: check, error_type

  implicit none

  !> Largest accepted number of SCF steps (measured: 14). The margin is needed:
  !> the stop criterion |dE| <= convThresh = 1e-12 sits at the round-off floor
  !> of eigensolves converged to tol 1e-8, so round-off (compiler flags, thread
  !> count, FMA) sets the step count (-O0 DiagonalizerList objects: 18 steps).
  integer(I32), parameter :: MAX_SCF_ITERATIONS = 30

  real(R64)    :: energyNew = 1e10_R64
  real(R64)    :: conv = 1e10_R64
  real(R64)    :: energyOld
  real(R64)    :: convThresh
  integer(I32) :: iStep, nTimeSteps, innerStep, nMatVecTotal
  type(error_type), allocatable :: error
  character(len=:), allocatable :: jsonFileName
  type(T_DiagonalizerList_FabricateInput) :: DiagonalizerListInput(1)
  logical :: allConvergedQ

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  jsonFileName = "test/simulationTests/Ne3d/T_Ne3d_17_GroundSolverScfLobpcgAce.json"
  call Json_LoadJsonFile(manualFileName_=jsonFileName, relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call SysKinetic_Fabricate
  call SysPotential_Fabricate
  call SysInteraction_Fabricate
  call Method_Fabricate
  call Orbs_Fabricate
  call OrbsInit_Fabricate
  call ConfigList_Fabricate
  call Coeffs_Fabricate
  call CoeffsInit_Fabricate
  call Mixing_Fabricate
  call GroundSolver_Fabricate
  DiagonalizerListInput(1) % ApplyMatOnVec => GroundSolver_Scf_HartreeFockAction
  DiagonalizerListInput(1) % ApplyPrecOnVec => GroundSolver_PreconditionAction
  DiagonalizerListInput(1) % dim = Grid_nPoints
  call DiagonalizerList_Fabricate(DiagonalizerListInput)

  call Say_Fabricate("program")
  convThresh = Json_Get("program.convThresh", 1e-13_R64)
  nTimeSteps = Json_Get("program.nTimeSteps", 10)

  call Grid_Setup
  call SysKinetic_Setup
  call SysPotential_Setup
  call SysInteraction_Setup
  call Orbs_Setup
  call OrbsInit_Setup
  call ConfigList_Setup
  call Coeffs_Setup
  call CoeffsInit_Setup
  call Method_Setup
  call GroundSolver_Setup
  call DiagonalizerList_Setup

  ! LOBPCG works in the grid metric (diagonal FEDVR quadrature weights)
  DiagonalizerList(1) % e % metricWeights = Grid_Ylm_weights

  !==========================================
  call Say_Section("perform action")
  !==========================================

  iStep = 0
  nMatVecTotal = 0
  allConvergedQ = .true.

  do while (abs(conv) > convThresh)
    iStep = iStep + 1
    if (iStep > MAX_SCF_ITERATIONS) error stop "T_Ne3d_17_GroundSolverScfLobpcgAce failure: SCF not converged within "// &
      "MAX_SCF_ITERATIONS steps"

    ! Explicit loop to execute approach nTimeSteps times
    do innerStep = 1, nTimeSteps
      call GroundSolver_Approach(Method_state, time=0.0_R64)
    end do

    energyOld = energyNew
    energyNew = Method_GetEnergy(0.0_R64)

    conv = energyNew - energyOld
    select type (e => DiagonalizerList(1) % e)
    type is (T_DiagonalizerList_E_Lobpcg)
      nMatVecTotal = nMatVecTotal + e % nMatVec
      allConvergedQ = allConvergedQ .and. e % convergedQ
      write (*, *) "convergence: ", energyNew, conv, " step:", iStep, " lobpcg iterations:", e % nIter, &
        " matvecs:", e % nMatVec, " converged: ", e % convergedQ
    class default
      error stop "T_Ne3d_17_GroundSolverScfLobpcgAce: diagonalizerList must be lobpcg"
    end select
  end do

  !==========================================
  call Say_Section("printout results")
  !==========================================

  print *
  print *, "final energy: ", energyNew
  print *, "SCF iterations: ", iStep, ", Fock matrix-vector products: ", nMatVecTotal

  call check(error, energyNew, -127.7449896_R64, thr=1e-8_R64)
  if (allocated(error)) error stop "T_Ne3d_17_GroundSolverScfLobpcgAce failure"

  call check(error, allConvergedQ, "every LOBPCG solve must converge")
  if (allocated(error)) error stop "T_Ne3d_17_GroundSolverScfLobpcgAce failure: a LOBPCG solve did not converge"

  call check(error, iStep <= MAX_SCF_ITERATIONS, "too many SCF steps")
  if (allocated(error)) error stop "T_Ne3d_17_GroundSolverScfLobpcgAce failure: SCF iteration count"

  !==========================================
  call Say_Goodbye
  !==========================================

end program
