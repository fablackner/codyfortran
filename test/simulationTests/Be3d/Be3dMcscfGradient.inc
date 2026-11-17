!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MCSCF ground state of Be with groundSolver.mcscf.gradient, shared by
! T_Be3d_01 and T_Be3d_02; included after `contains`.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Loads jsonFileName and repeats GroundSolver_Approach (CI eigenvector by
!> the first diagonalizer, then the orbital step) until the energy changes
!> by less than program.convThresh and the orbital residual norm is below
!> program.gradThresh, at most program.maxIter times. Stops unless the final
!> energy agrees with referenceEnergy to 1e-9, the residual norm is below
!> gradThresh and at most maxIterations iterations were needed. The modules
!> stay set up, so the caller can save the ground state.
subroutine CheckMcscfGradientGroundState(jsonFileName, referenceEnergy, maxIterations)
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

  character(len=*), intent(in) :: jsonFileName
  real(R64), intent(in) :: referenceEnergy
  integer(I32), intent(in) :: maxIterations

  real(R64) :: energy, energyOld, convThresh, gradThresh
  integer(I32) :: maxIter, iter
  type(error_type), allocatable :: error
  type(T_DiagonalizerList_FabricateInput) :: DiagonalizerListInput(1)

  !==========================================
  call Say_Section("start")
  !==========================================

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

  call check(error, energy, referenceEnergy, thr=1e-9_R64)
  if (allocated(error)) error stop "CheckMcscfGradientGroundState: energy off the reference"

  call check(error, GroundSolver_Mcscf_Gradient_residualNorm < gradThresh)
  if (allocated(error)) error stop "CheckMcscfGradientGroundState: orbital gradient above program.gradThresh"

  call check(error, iter <= maxIterations)
  if (allocated(error)) error stop "CheckMcscfGradientGroundState: more than maxIterations iterations"

end subroutine
