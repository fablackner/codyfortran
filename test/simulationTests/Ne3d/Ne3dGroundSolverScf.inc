!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Hartree-Fock ground state of Ne with groundSolver.scf, shared by T_Ne3d_02
! (potential mixing) and T_Ne3d_05 (density mixing); included after
! `contains`.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Loads jsonFileName and repeats GroundSolver_Approach (program.nTimeSteps
!> calls per step, the Fock operator GroundSolver_Scf_HartreeFockAction
!> diagonalized by the first diagonalizer) until the energy changes by less
!> than program.convThresh in one step. Stops unless the final energy agrees
!> with referenceEnergy to 1e-8.
subroutine CheckScfGroundState(jsonFileName, referenceEnergy)
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
  use M_GroundSolver_Scf
  use testdrive, only: check, error_type

  character(len=*), intent(in) :: jsonFileName
  real(R64), intent(in) :: referenceEnergy

  real(R64) :: energyNew, energyOld, conv, convThresh
  integer(I32) :: iStep, nTimeSteps, innerStep
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
  DiagonalizerListInput(1) % ApplyMatOnVec => GroundSolver_Scf_HartreeFockAction
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

  !==========================================
  call Say_Section("perform action")
  !==========================================

  iStep = 0
  energyNew = 1e10_R64
  conv = 1e10_R64

  do while (abs(conv) > convThresh)
    iStep = iStep + 1

    do innerStep = 1, nTimeSteps
      call GroundSolver_Approach(Method_state, time=0.0_R64)
    end do

    energyOld = energyNew
    energyNew = Method_GetEnergy(0.0_R64)

    conv = energyNew - energyOld
    write (*, *) "convergence: ", energyNew, conv, " step:", iStep
  end do

  !==========================================
  call Say_Section("printout results")
  !==========================================

  print *
  print *, "final energy: ", energyNew

  call check(error, energyNew, referenceEnergy, thr=1e-8_R64)
  if (allocated(error)) error stop "CheckScfGroundState: energy off the reference"

end subroutine
