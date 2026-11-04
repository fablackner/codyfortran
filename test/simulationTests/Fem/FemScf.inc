!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Helpers shared by the SCF tests T_Fem_03_H2ScfRhf and T_Fem_04_H2oScfRhf
! (included after `contains`; every procedure here is used by each of them).
!
! - SetupFemScf:  fabricate and set up the restricted Hartree-Fock run of a
!                 JSON file (method tdhx, groundSolver.scf.stdImpl and the
!                 preconditioned LOBPCG in the FEM metric)
! - RunScf:       SCF loop until the energy and the orbital energies converge
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Fabricates and sets up the SCF run of jsonFileName: the eigensolver of
!> groundSolver.scf.stdImpl gets the Fock action GroundSolver_Scf_HartreeFockAction
!> and the preconditioner GroundSolver_PreconditionAction ((T + max(-theta,
!> sigma_min))^-1 by the exact FDM solve), and works in the FEM metric
!> (metricWeights = Grid_Fem_weights: diagonal mass, BLAS-3 Gram matrices).
!> The start orbitals are orthonormalized at the end (OrbsInit only
!> normalizes).
subroutine SetupFemScf(jsonFileName)
  use M_Utils_Json
  use M_Grid
  use M_Grid_Fem
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

  character(len=*), intent(in) :: jsonFileName

  type(T_DiagonalizerList_FabricateInput) :: DiagonalizerListInput(1)

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

  ! LOBPCG works in the FEM metric (diagonal mass, BLAS-3 Gram matrices)
  DiagonalizerList(1) % e % metricWeights = Grid_Fem_weights

  ! OrbsInit only normalizes: orthonormalize the start orbitals
  call Orbs_Orthonormalize(Orbs_orbs)

end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> SCF loop from the current orbitals, one GroundSolver_Approach per step,
!> until |dE| < convThresh and max |d eps| < epsThresh (the orbital energies
!> converge only like sqrt(dE)); stops after maxScfIterations steps. energy
!> and eps (the size(eps) lowest orbital energies) hold the values of the
!> start orbitals on input (huge if unknown) and E_el and the orbital
!> energies on output. Prints one line per step ("SCF (label)") and a
!> summary. Also returns the number of steps, the largest LOBPCG iteration
!> count of the warm-started steps (all steps if firstStepWarmQ, else all but
!> the first) and whether every LOBPCG solve converged.
subroutine RunScf(nScf, maxWarmIter, allConvergedQ, energy, eps, label, maxScfIterations, convThresh, epsThresh, &
                  firstStepWarmQ)
  use M_Method
  use M_DiagonalizerList
  use M_DiagonalizerList_Lobpcg, only: T_DiagonalizerList_E_Lobpcg
  use M_GroundSolver

  integer(I32), intent(out) :: nScf
  integer(I32), intent(out) :: maxWarmIter
  logical, intent(out) :: allConvergedQ
  real(R64), intent(inout) :: energy
  real(R64), intent(inout) :: eps(:)
  character(len=*), intent(in) :: label
  integer(I32), intent(in) :: maxScfIterations
  real(R64), intent(in) :: convThresh
  real(R64), intent(in) :: epsThresh
  logical, intent(in) :: firstStepWarmQ

  real(R64), allocatable :: epsOld(:)
  real(R64) :: energyOld
  integer(I64) :: count0, count1, countRate
  integer(I32) :: nIter, nMatVec, nMatVecTotal
  logical :: convergedQ

  write (*, *)
  write (*, *) "SCF (", label, ")"

  call system_clock(count0, countRate)
  nScf = 0
  maxWarmIter = 0
  nMatVecTotal = 0
  allConvergedQ = .true.

  do
    nScf = nScf + 1
    if (nScf > maxScfIterations) error stop "RunScf (FemScf.inc): SCF not converged within program.maxScfIterations"

    energyOld = energy
    epsOld = eps
    call GroundSolver_Approach(Method_state, time=0.0_R64)
    energy = Method_GetEnergy(0.0_R64)
    eps(:) = DiagonalizerList(1) % e % evals(1:size(eps))

    select type (e => DiagonalizerList(1) % e)
    type is (T_DiagonalizerList_E_Lobpcg)
      nIter = e % nIter
      nMatVec = e % nMatVec
      convergedQ = e % convergedQ
    class default
      error stop "RunScf (FemScf.inc): diagonalizerList must be lobpcg"
    end select
    allConvergedQ = allConvergedQ .and. convergedQ
    if (nScf > 1 .or. firstStepWarmQ) maxWarmIter = max(maxWarmIter, nIter)
    nMatVecTotal = nMatVecTotal + nMatVec

    write (*, '(1X, A, I3, A, F18.12, A, ES10.2, A, ES10.2, A, I3, A, I4, A, L1)') "scf ", nScf, ": E_el ", energy, &
      ", dE ", energy - energyOld, ", max|deps| ", maxval(abs(eps - epsOld)), ", lobpcg iterations ", nIter, &
      ", matvecs ", nMatVec, ", converged ", convergedQ

    if (abs(energy - energyOld) < convThresh .and. maxval(abs(eps - epsOld)) < epsThresh) exit
  end do

  call system_clock(count1)
  write (*, '(1X, 3A, I0, A, I0, A, F8.2, A)') "SCF (", label, "): ", nScf, " iterations, ", nMatVecTotal, &
    " Fock matrix-vector products, ", real(count1 - count0, R64) / real(countRate, R64), " s"

end subroutine
