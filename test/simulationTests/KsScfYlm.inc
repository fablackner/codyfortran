!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Kohn-Sham SCF of a closed-shell atom on a Ylm grid, shared by the DFT
! ground-state tests of the atom folders (included after `contains`):
!
! - RunKsScfYlm: fabricate and set up the run of a JSON file
!   (method.mb.orbBased.tdks, sysXc.ylm, groundSolver.scf.stdImpl with the
!   Kohn-Sham mean field, LOBPCG in the Ylm metric with the exact kinetic
!   preconditioner of the FEDVR blocks), iterate to |dE| < program.convThresh
!   and return the Kohn-Sham energy, its parts and the orbital energies
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine RunKsScfYlm(energy, eps, nScf, jsonFileName, maxScfIterations)
  use M_Utils_Json
  use M_Grid
  use M_SysKinetic
  use M_SysPotential
  use M_SysInteraction
  use M_SysXc
  use M_Orbs
  use M_OrbsInit
  use M_ConfigList
  use M_Coeffs
  use M_CoeffsInit
  use M_FieldFree
  use M_Method
  use M_Method_Mb_OrbBased_Tdks
  use M_DiagonalizerList
  use M_DiagonalizerList_Lobpcg, only: T_DiagonalizerList_E_Lobpcg
  use M_Mixing
  use M_GroundSolver
  use M_GroundSolver_Scf

  real(R64), intent(out) :: energy
  real(R64), intent(out), allocatable :: eps(:)
  integer(I32), intent(out) :: nScf
  character(len=*), intent(in) :: jsonFileName
  integer(I32), intent(in) :: maxScfIterations

  type(T_DiagonalizerList_FabricateInput) :: DiagonalizerListInput(1)
  real(R64) :: energyOld, convThresh
  integer(I64) :: count0, count1, countRate

  call Json_LoadJsonFile(manualFileName_=jsonFileName, relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call SysKinetic_Fabricate
  call SysPotential_Fabricate
  call SysInteraction_Fabricate
  call SysXc_Fabricate
  call FieldFree_Fabricate
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
  convThresh = Json_Get("program.convThresh", 1e-11_R64)

  call Grid_Setup
  call SysKinetic_Setup
  call SysPotential_Setup
  call SysInteraction_Setup
  call SysXc_Setup
  call Orbs_Setup
  call OrbsInit_Setup
  call ConfigList_Setup
  call Coeffs_Setup
  call CoeffsInit_Setup
  call FieldFree_Setup
  call Method_Setup
  call GroundSolver_Setup
  call DiagonalizerList_Setup

  ! LOBPCG in the Ylm metric; OrbsInit only normalizes
  DiagonalizerList(1) % e % metricWeights = Grid_weights
  call Orbs_Orthonormalize(Orbs_orbs)

  energy = huge(1.0_R64)
  nScf = 0
  do
    nScf = nScf + 1
    if (nScf > maxScfIterations) error stop "RunKsScfYlm (KsScfYlm.inc): SCF not converged"
    energyOld = energy
    call system_clock(count0, countRate)
    call GroundSolver_Approach(Method_state, time=0.0_R64)
    energy = Method_GetEnergy(0.0_R64)
    call system_clock(count1)
    select type (e => DiagonalizerList(1) % e)
    type is (T_DiagonalizerList_E_Lobpcg)
      write (*, '(1X, A, I3, A, F20.12, A, ES10.2, A, F16.10, A, F16.10, A, I4, A, F8.3, A)') "ks scf ", nScf, ": E ", &
        energy, ", dE ", energy - energyOld, ", E_H ", Method_Mb_OrbBased_Tdks_hartreeEnergy, ", E_xc ", &
        Method_Mb_OrbBased_Tdks_xcEnergy, ", lobpcg iterations ", e % nIter, ", ", &
        real(count1 - count0, kind=R64) / real(countRate, kind=R64), " s"
    class default
      error stop "RunKsScfYlm (KsScfYlm.inc): diagonalizerList must be lobpcg"
    end select
    if (abs(energy - energyOld) < convThresh) exit
  end do

  eps = DiagonalizerList(1) % e % evals(1:Orbs_nOrbsInState)
  write (*, '(1X, A, *(F14.8))') "orbital energies: ", eps

end subroutine
