!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> MCSCF ground state of Be (5 restricted orbitals 1s 2s 2p of definite m,
!> full CI in the total-M = 0 sector, Ylm grid) with the ground solver
!> groundSolver.mcscf.gradient: the same state as the imaginary-time MCTDHF
!> program of ../mctdhf, without its time-step bias (any atom with the
!> corresponding input).
!>
!> Start: the Hartree-Fock orbitals 1s 2s of ../hf (orbsInit.load with
!> nOrbsIn 2 and lmaxIn 0) and hydrogen-like 2p (orbsInit.ylm.hydrogenLike,
!> listing only the 2p). Every iteration (GroundSolver_Approach) starts with
!> the CI eigenvector of the current orbitals (diagonalizerList(1)), so the
!> coeffsInit vector only seeds the first CI solve.
!>
!> Converged when the change of the MCSCF energy per iteration
!> (GroundSolver_Mcscf_Gradient_energy) is below program.convThresh and the
!> orbital energy gradient below program.gradThresh (at most program.maxIter
!> iterations). The final energy is that of the saved state.
!>
!> Output: orbXX_YY.in, coeffs.in, rdm2BtXX_YY.in (as ../mctdhf, so that the
!>         dynamics can load either), relax.dat (iteration, energy, gradient
!>         norm)
program P_GroundMcscf
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_RdmDiagonalize
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
  use M_Method_Mb_OrbBased
  use M_DiagonalizerList
  use M_Mixing
  use M_GroundSolver
  use M_GroundSolver_Mcscf
  use M_GroundSolver_Mcscf_Gradient

  implicit none

  real(R64)    :: energy, energyOld, convThresh, gradThresh
  integer(I32) :: maxIter, iter, io
  real(R64), allocatable :: natOcc(:)
  complex(R64), allocatable :: rdm1(:, :), natOrbs(:, :)
  type(T_DiagonalizerList_FabricateInput) :: DiagonalizerListInput(1)

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
  DiagonalizerListInput(1) % ApplyMatOnVec => GroundSolver_Mcscf_HamiltonianAction
  DiagonalizerListInput(1) % dim = Coeffs_nCoeffs
  call DiagonalizerList_Fabricate(DiagonalizerListInput)

  call Say_Fabricate("program")
  convThresh = Json_Get("program.convThresh", 1e-12_R64)
  gradThresh = Json_Get("program.gradThresh", 1e-6_R64)
  maxIter = Json_Get("program.maxIter", 1000)

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
  call Say_Section("mcscf iteration")
  !==========================================

  open (newunit=io, file="relax.dat", status="replace")

  energy = Method_GetEnergy(0.0_R64)
  write (*, '(A,F22.14)') " start energy: ", energy

  do iter = 1, maxIter

    call GroundSolver_Approach(Method_state, 0.0_R64)

    energyOld = energy
    energy = GroundSolver_Mcscf_Gradient_energy

    write (io, '(I8, F24.16, ES14.6)') iter, energy, GroundSolver_Mcscf_Gradient_residualNorm
    flush (io)
    write (*, '(A,I6,A,F22.14,A,ES10.2,A,ES10.2)') " iteration", iter, "  energy =", energy, &
      "  change =", energy - energyOld, "  gradient =", GroundSolver_Mcscf_Gradient_residualNorm

    if (abs(energy - energyOld) < convThresh .and. GroundSolver_Mcscf_Gradient_residualNorm < gradThresh) exit

  end do

  close (io)

  if (maxIter < iter) write (*, '(A,I0,A)') " warning: not converged within ", maxIter, " iterations"

  !==========================================
  call Say_Section("printout results")
  !==========================================

  call Method_Mb_OrbBased_FillRdm1(rdm1, Method_state)
  call RdmDiagonalize_Rdm1(natOcc, natOrbs, rdm1(1:Orbs_nOrbsInState, 1:Orbs_nOrbsInState))

  energy = Method_GetEnergy(0.0_R64)

  print *
  write (*, '(A,F22.14)') " final energy: ", energy
  write (*, '(A,*(ES12.4))') " natural occupations (spin up): ", natOcc

  !==========================================
  call Say_Section("save data")
  !==========================================

  call Orbs_SaveOrbs(Orbs_orbs)
  call Coeffs_SaveCoeffs(Coeffs_coeffs)
  call Coeffs_SaveTwoRdm(Coeffs_coeffs)

  !==========================================
  call Say_Goodbye
  !==========================================

end program
