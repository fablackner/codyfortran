!> ctest-requires: libxc
!> Kohn-Sham SCF (groundSolver.scf.stdImpl with sysXc.linear, dense LAPACK
!> diagonalization) of 1D helium with soft-Coulomb interactions
!> 1/sqrt(x^2 + 1) and the matching 1D LDA of Libxc (lda_x_1d_soft,
!> lda_c_1d_csc with beta = 1); the energy is the reference of the
!> imaginary-time TDKS of T_He1d_13.
program T_He1d_12_GroundSolverKsLda1d
  use M_Utils_Types
  use M_Utils_Say
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
  use M_Method
  use M_Method_Mb_OrbBased_Tdks
  use M_DiagonalizerList
  use M_Mixing
  use M_GroundSolver
  use M_GroundSolver_Scf
  use testdrive, only: check, error_type

  implicit none

  type(error_type), allocatable :: error
  type(T_DiagonalizerList_FabricateInput) :: DiagonalizerListInput(1)
  real(R64) :: energy, energyOld, convThresh
  integer(I32) :: iScf

  call Say_Hello
  call Json_LoadJsonFile(manualFileName_="test/simulationTests/He1d/T_He1d_12_GroundSolverKsLda1d.json", &
                         relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call SysKinetic_Fabricate
  call SysPotential_Fabricate
  call SysInteraction_Fabricate
  call SysXc_Fabricate
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
  convThresh = Json_Get("program.convThresh", 1e-12_R64)

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
  call Method_Setup
  call GroundSolver_Setup
  call DiagonalizerList_Setup
  call Orbs_Orthonormalize(Orbs_orbs)

  energy = huge(1.0_R64)
  do iScf = 1, 100
    energyOld = energy
    call GroundSolver_Approach(Method_state, time=0.0_R64)
    energy = Method_GetEnergy(0.0_R64)
    write (*, '(1X, A, I3, A, F20.14, A, ES10.2, A, 2F14.10)') "ks scf ", iScf, ": E ", energy, ", dE ", &
      energy - energyOld, ", E_H, E_xc ", Method_Mb_OrbBased_Tdks_hartreeEnergy, Method_Mb_OrbBased_Tdks_xcEnergy
    if (abs(energy - energyOld) < convThresh) exit
  end do

  print *
  print *, "1D LDA Kohn-Sham energy: ", energy, ", eps ", DiagonalizerList(1) % e % evals(1)

  call check(error, abs(energy - energyOld) < convThresh)
  if (allocated(error)) error stop "T_He1d_12_GroundSolverKsLda1d failure: SCF not converged"
  call check(error, energy, -2.2014358365_R64, thr=1e-8_R64)
  if (allocated(error)) error stop "T_He1d_12_GroundSolverKsLda1d failure: energy (regression)"

  call Say_Goodbye

end program
