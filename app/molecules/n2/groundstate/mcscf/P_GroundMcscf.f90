!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> MCSCF ground state of N2 (R = 2.068) on the prolate spheroidal grid: 10
!> restricted orbitals (the 7 occupied of ../scf plus the valence 1πg and
!> 3σu), full CI of 14 electrons (14 400 configurations), ground solver
!> groundSolver.mcscf.gradient (the program app/atoms/P_GroundMcscf.f90
!> plus the spectral-anchor output).
!>
!> Start: the Hartree-Fock orbitals of ../scf (orbsInit.load with nOrbsIn 7)
!> and LCAO 2p combinations for 1πg (m = ±1, gerade) and 3σu (m = 0,
!> ungerade) (orbsInit.prolate.lcao, listing only these three). Every
!> iteration (GroundSolver_Approach) starts with the CI eigenvector of the
!> current orbitals (diagonalizerList(1), ARPACK), so the coeffsInit vector
!> only seeds the first CI solve.
!>
!> Converged when the change of the MCSCF energy per iteration
!> (GroundSolver_Mcscf_Gradient_energy) is below program.convThresh and the
!> orbital energy gradient below program.gradThresh (at most program.maxIter
!> iterations). The final energy is that of the saved state.
!>
!> Output: orbXX_YY.in, coeffs.in, rdm2BtXX_YY.in (read by
!>         ../../dynamics/staticOrbs), relax.dat (iteration, energy, gradient
!>         norm); with a program.anchor block the reference of an anchored
!>         TD-2RDM reconstruction (reconstruction.singlet.anchor):
!>         rdm3Bt01_01_02.in, the up-up-down 3-RDM, and with
!>         program.anchor.nTransitions = K > 0 the transition RDMs
!>         <Psi|..|chi_k> (trdm2Bt01_02.in, trdm3Bt01_01_02.in) of the K
!>         lowest Ritz vectors chi_k of H in the Krylov space
!>         span{Z Psi, H Z Psi, ...} (program.anchor.nLanczos vectors, default
!>         60, orthogonal to Psi; Z = sum z_ij a+_i a_j the dipole along the
!>         bond, H the CI Hamiltonian of the final orbitals): the spectral
!>         anchor, which makes the TD-2RDM linear response to a field along
!>         the bond exact for these (nearly exact) excited states.
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
  logical      :: saveAnchorQ
  integer(I32) :: nTransitions
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
  saveAnchorQ = Json_GetExistence("program.anchor")
  nTransitions = Json_Get("program.anchor.nTransitions", 0)

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
  if (saveAnchorQ) call SaveAnchor(nTransitions)

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  !> Reference 3-RDM and, for nTransitions > 0, the transition RDMs of the
  !> Krylov vectors of Z Psi (header).
  subroutine SaveAnchor(nTransitions)
    use M_Utils_LapackLib
    use M_Reconstruction_Singlet_Anchor
    use M_Grid_Prolate

    integer(I32), intent(in) :: nTransitions

    complex(R64), allocatable :: kets(:, :), zMat(:, :), zPsi(:), x(:), q(:, :), hq(:, :), hProj(:, :)
    complex(R64), allocatable :: ritzVectors(:, :)
    real(R64), allocatable :: ritzValues(:), z(:)
    integer(I32) :: i, k, pass, nLanczos

    nLanczos = Json_Get("program.anchor.nLanczos", 60)
    allocate (kets(Coeffs_nCoeffs, nTransitions + 1))
    kets(:, 1) = Coeffs_coeffs

    if (nTransitions > 0) then
      call Grid_Prolate_FillNodeCoordinate(z, [0.0_R64, 0.0_R64, 1.0_R64])
      call Method_Mb_OrbBased_FillOneBodyMatrix(zMat, Orbs_orbs, spread(z, 2, Orbs_nOrbsInState) * Orbs_orbs)

      ! Lanczos vectors q of Z Psi (orthonormal, orthogonal to Psi; full
      ! reorthogonalization) and the projected Hamiltonian q^+ H q
      nLanczos = max(nLanczos, nTransitions)
      allocate (q(Coeffs_nCoeffs, nLanczos), hq(Coeffs_nCoeffs, nLanczos), zPsi(Coeffs_nCoeffs))
      zPsi = 0.0_R64
      call Coeffs_ApplyH1FillRdm1(Coeffs_coeffs, dCoeffs_=zPsi, h1_=zMat)
      x = zPsi
      do k = 1, nLanczos
        do pass = 1, 2
          x = x - Coeffs_coeffs * dot_product(Coeffs_coeffs, x)
          do i = 1, k - 1
            x = x - q(:, i) * dot_product(q(:, i), x)
          end do
        end do
        q(:, k) = x / norm2(abs(x))
        call Method_Mb_OrbBased_ApplyCiHamiltonian(hq(:, k), q(:, k), 0.0_R64)
        x = hq(:, k)
      end do
      hProj = matmul(conjg(transpose(q)), hq)
      hProj = 0.5_R64 * (hProj + conjg(transpose(hProj)))
      call LapackLib_DiagonalizeGeneric(ritzValues, ritzVectors, hProj, .true.)

      ! the nTransitions lowest Ritz vectors (converge first to the low-lying
      ! dipole-allowed eigenstates)
      write (*, '(A)') " spectral anchor: Ritz states (excitation energy, |<chi|Z|Psi>|^2)"
      do k = 1, nTransitions
        kets(:, k + 1) = matmul(q, ritzVectors(:, k))
        write (*, '(I6, F14.6, ES14.4)') k, ritzValues(k) - energy, abs(dot_product(kets(:, k + 1), zPsi))**2
      end do
    end if

    call Reconstruction_Singlet_Anchor_SaveReference(kets)

  end subroutine

end program
