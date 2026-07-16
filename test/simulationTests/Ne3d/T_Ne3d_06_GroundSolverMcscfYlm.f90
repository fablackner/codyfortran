program T_Ne3d_06_GroundSolverMcscfYlm
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_Actions
  use M_Grid
  use M_Grid_Ylm
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
  use M_GroundSolver_Mcscf_YlmOpt
  use testdrive, only: check, error_type

  implicit none

  real(R64)    :: energyNew = 1e10_R64
  real(R64)    :: conv = 1e10_R64
  real(R64)    :: energyOld
  real(R64)    :: convThresh
  real(R64)    :: time
  integer(I32) :: iStep, nTimeSteps, innerStep
  type(error_type), allocatable :: error
  character(len=:), allocatable :: jsonFileName
  type(T_DiagonalizerList_FabricateInput) :: DiagonalizerListInput(3)

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  jsonFileName = "test/simulationTests/Ne3d/T_Ne3d_06_GroundSolverMcscfYlm.json"
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
  DiagonalizerListInput(2) % ApplyMatOnVec => ApplyMatOnVecL0
  DiagonalizerListInput(2) % dim = Grid_Ylm_nRadial
  DiagonalizerListInput(3) % ApplyMatOnVec => ApplyMatOnVecL1
  DiagonalizerListInput(3) % dim = Grid_Ylm_nRadial
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

  !==========================================
  call Say_Section("perform action")
  !==========================================

  iStep = 0
  time = 0.0_R64  ! Initialize simulation time

  do while (abs(conv) > convThresh)
    iStep = iStep + 1

    ! Explicit loop to execute approach nTimeSteps times
    do innerStep = 1, nTimeSteps
      call GroundSolver_Approach(Method_state, time=0.0_R64)
    end do

    energyOld = energyNew
    energyNew = Method_GetEnergy(time)  ! Use current time for energy calculation

    conv = energyNew - energyOld
    write (*, *) "convergence: ", energyNew, conv, " step:", iStep
  end do

  !==========================================
  call Say_Section("printout results")
  !==========================================

  print *
  print *, "final energy: ", energyNew

  call check(error, energyNew, -127.7449896_R64, thr=1e-8_R64)
  if (allocated(error)) error stop "T_Ne3d_06_GroundSolverMcscfYlm failure"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  subroutine ApplyMatOnVecL0(dOrbLm, orbLm, time)
    use M_GroundSolver_Mcscf_YlmOpt

    complex(R64), intent(out), contiguous, target :: dOrbLm(:)
    complex(R64), intent(in), contiguous, target :: orbLm(:)
    real(R64), intent(in) :: time

    call GroundSolver_Mcscf_YlmOpt_FockAction(dOrbLm, orbLm, l=0, time=time)

  end subroutine

  subroutine ApplyMatOnVecL1(dOrbLm, orbLm, time)
    use M_GroundSolver_Mcscf_YlmOpt

    complex(R64), intent(out), contiguous, target :: dOrbLm(:)
    complex(R64), intent(in), contiguous, target :: orbLm(:)
    real(R64), intent(in) :: time

    call GroundSolver_Mcscf_YlmOpt_FockAction(dOrbLm, orbLm, l=1, time=time)

  end subroutine

end program
