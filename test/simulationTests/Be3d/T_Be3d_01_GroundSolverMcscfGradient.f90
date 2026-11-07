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
program T_Be3d_01_GroundSolverMcscfGradient
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Grid
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
  use M_Mixing
  use M_GroundSolver
  use M_GroundSolver_Mcscf
  use M_GroundSolver_Mcscf_Gradient
  use testdrive, only: check, error_type

  implicit none

  real(R64)    :: energy, energyOld, convThresh, gradThresh
  integer(I32) :: maxIter, iter
  type(error_type), allocatable :: error
  character(len=:), allocatable :: jsonFileName
  type(T_DiagonalizerList_FabricateInput) :: DiagonalizerListInput(1)

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  jsonFileName = "test/simulationTests/Be3d/T_Be3d_01_GroundSolverMcscfGradient.json"
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
  DiagonalizerListInput(1) % ApplyMatOnVec => GroundSolver_Mcscf_HamiltonianAction
  DiagonalizerListInput(1) % dim = Coeffs_nCoeffs
  call DiagonalizerList_Fabricate(DiagonalizerListInput)

  call Say_Fabricate("program")
  convThresh = Json_Get("program.convThresh", 1e-12_R64)
  gradThresh = Json_Get("program.gradThresh", 1e-6_R64)
  maxIter = Json_Get("program.maxIter", 100)

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

  call Orbs_Orthonormalize(Orbs_orbs)

  !==========================================
  call Say_Section("perform action")
  !==========================================

  energy = Method_GetEnergy(0.0_R64)

  do iter = 1, maxIter
    call GroundSolver_Approach(Method_state, 0.0_R64)
    energyOld = energy
    energy = GroundSolver_Mcscf_Gradient_energy
    write (*, '(A,I4,A,F20.12,A,ES10.2)') " iteration", iter, "  energy =", energy, &
      "  gradient =", GroundSolver_Mcscf_Gradient_residualNorm
    if (abs(energy - energyOld) < convThresh .and. GroundSolver_Mcscf_Gradient_residualNorm < gradThresh) exit
  end do

  !==========================================
  call Say_Section("printout results")
  !==========================================

  energy = Method_GetEnergy(0.0_R64)

  print *
  print *, "final energy: ", energy

  call check(error, energy, -14.6168557502_R64, thr=1e-9_R64)
  if (allocated(error)) error stop "T_Be3d_01_GroundSolverMcscfGradient failure: energy"

  call check(error, GroundSolver_Mcscf_Gradient_residualNorm < gradThresh)
  if (allocated(error)) error stop "T_Be3d_01_GroundSolverMcscfGradient failure: gradient"

  call check(error, iter <= 20)
  if (allocated(error)) error stop "T_Be3d_01_GroundSolverMcscfGradient failure: more than 20 iterations"

  !==========================================
  call Say_Goodbye
  !==========================================

end program
