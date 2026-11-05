!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Hartree-Fock ground state of a closed-shell atom on a Ylm grid with the
!> l-channel SCF solver (groundSolver.scf.ylmOpt): per iteration the radial
!> Fock operator of every l = 0..lmax (lmax <= 1, s and p orbitals) is
!> diagonalized with its own diagonalizer (diagonalizerList entries in the
!> order l = 0, 1), and the eigenvectors are placed in the (l, m) channels of
!> orbsInit.ylm.hydrogenLike (n, l, m per orbital). Iterates until the energy
!> changes by less than program.convThresh per program.nTimeSteps steps.
!>
!> Output: orbXX_YY.in, coeffs.in
program P_GroundScfYlm
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
  use M_GroundSolver_Scf
  use M_GroundSolver_Scf_YlmOpt

  implicit none

  real(R64)    :: energyNew = 1e10_R64
  real(R64)    :: conv = 1e10_R64
  real(R64)    :: energyOld
  real(R64)    :: convThresh
  real(R64)    :: time
  integer(I32) :: iStep, nTimeSteps, innerStep
  type(T_DiagonalizerList_FabricateInput) :: DiagonalizerListInput(2)
  integer(I32) :: nChannelsL

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call Json_LoadJsonFile()

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
  DiagonalizerListInput(1) % ApplyMatOnVec => ApplyMatOnVecL0
  DiagonalizerListInput(1) % dim = Grid_Ylm_nRadial
  DiagonalizerListInput(2) % ApplyMatOnVec => ApplyMatOnVecL1
  DiagonalizerListInput(2) % dim = Grid_Ylm_nRadial
  ! one radial diagonalizer per l channel
  nChannelsL = Grid_Ylm_lmax + 1
  if (2 < nChannelsL) error stop "P_GroundScfYlm: grid.ylm.lmax must be 0 or 1"
  call DiagonalizerList_Fabricate(DiagonalizerListInput(1:nChannelsL))

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

  !==========================================
  call Say_Section("save data")
  !==========================================

  call Coeffs_SaveCoeffs(Coeffs_coeffs)
  call Orbs_SaveOrbs(Orbs_orbs)

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  subroutine ApplyMatOnVecL0(dOrbLm, orbLm, time)
    use M_Grid_Ylm
    use M_GroundSolver_Scf

    complex(R64), intent(out), contiguous, target :: dOrbLm(:)
    complex(R64), intent(in), contiguous, target :: orbLm(:)
    real(R64), intent(in) :: time

    call GroundSolver_Scf_YlmOpt_HartreeFockAction(dOrbLm, orbLm, l=0, time=time)

  end subroutine

  subroutine ApplyMatOnVecL1(dOrbLm, orbLm, time)
    use M_Grid_Ylm
    use M_GroundSolver_Scf

    complex(R64), intent(out), contiguous, target :: dOrbLm(:)
    complex(R64), intent(in), contiguous, target :: orbLm(:)
    real(R64), intent(in) :: time

    call GroundSolver_Scf_YlmOpt_HartreeFockAction(dOrbLm, orbLm, l=1, time=time)

  end subroutine

end program
