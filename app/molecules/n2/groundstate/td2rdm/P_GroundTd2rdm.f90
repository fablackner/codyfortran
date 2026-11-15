!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> TD-2RDM ground state of N2 in a fixed active space by adiabatic switching
!> (td2rdmStaticOrbs, no CI anywhere): a test of the closure, since the exact
!> answer in these orbitals is the CI energy of ../mcscf.
!>
!> The active orbitals are those of orbsInit (the MCSCF orbitals of ../mcscf).
!> The program rotates them to the Hartree-Fock orbitals of the active space
!> (Fock operator of the first nBodies orbitals, iterated in the matrices h1
!> and h2 of the method), so that the determinant of twoRdmInit.excited is an
!> eigenstate of h + U0, and switches the correlation on over -T < t < 0
!> (method.mb.orbBased.adiabaticSwitching.switchTime T). program.relaxTime of
!> propagation with the physical Hamiltonian follows, which shows how
!> stationary the state is under the TD-2RDM equations.
!>
!> Output <program.outputName>.dat: t, lambda, energy <H(t)>, min eigenvalues
!> of D and Q (up-down block), natural occupations (spin up); with
!> program.saveQ the rotated orbitals and the 2-RDM at t = 0 (orbXX_YY.in,
!> rdm2BtXX_YY.in).
program P_GroundTd2rdm
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_LapackLib
  use M_Utils_RdmBasisChange
  use M_Utils_RdmDiagonalize
  use M_Utils_RdmHole
  use M_Utils_RdmTrace
  use M_Grid
  use M_SysKinetic
  use M_SysPotential
  use M_SysInteraction
  use M_Orbs
  use M_OrbsInit
  use M_TwoRdm
  use M_TwoRdmInit
  use M_Reconstruction
  use M_Method
  use M_Method_Mb
  use M_Method_Mb_OrbBased
  use M_IntegratorList
  use M_Propagator

  implicit none

  real(R64)    :: time, relaxTime, outputStep, timeStep, energy0
  integer(I32) :: nTimeSteps, innerStep, io, nSteps
  integer(I64) :: c0, c1, cRate
  logical      :: saveQ
  character(len=:), allocatable :: outputName
  complex(R64), allocatable :: rdm2(:, :, :, :)
  type(T_IntegratorList_FabricateInput) :: IntegratorListInput(1)

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
  call TwoRdm_Fabricate
  call TwoRdmInit_Fabricate
  call Reconstruction_Fabricate
  IntegratorListInput(1) % TimeDerivative => Method_TimeDerivative
  call IntegratorList_Fabricate(IntegratorListInput)
  call Propagator_Fabricate()

  call Say_Fabricate("program")
  relaxTime = Json_Get("program.relaxTime", 10.0_R64)
  outputStep = Json_Get("program.outputStep", 0.2_R64)
  nTimeSteps = Json_Get("program.nTimeSteps", 5)
  outputName = Json_Get("program.outputName", "switching")
  saveQ = Json_Get("program.saveQ", .true.)
  if (.not. (Method_Mb_OrbBased_staticOrbsQ .and. Method_Mb_OrbBased_switchTime > 0.0_R64)) &
    error stop "P_GroundTd2rdm needs td2rdmStaticOrbs with method.mb.orbBased.adiabaticSwitching"

  call Grid_Setup
  call SysKinetic_Setup
  call SysPotential_Setup
  call SysInteraction_Setup
  call Orbs_Setup
  call OrbsInit_Setup
  call TwoRdm_Setup
  call TwoRdmInit_Setup
  call Method_Setup
  call Reconstruction_Setup
  call IntegratorList_Setup
  call Propagator_Setup

  !==========================================
  call Say_Section("Hartree-Fock orbitals of the active space")
  !==========================================

  call RotateToActiveSpaceHartreeFock
  call Method_Mb_OrbBased_SetSwitchingReference

  !==========================================
  call Say_Section("adiabatic switching and relaxation")
  !==========================================

  open (newunit=io, file=outputName//".dat", status="replace")
  timeStep = outputStep / nTimeSteps
  time = -Method_Mb_OrbBased_switchTime
  write (*, '(A,F22.12)') " Hartree-Fock determinant: energy ", Method_GetEnergy(0.0_R64)
  call WriteProgress(time)
  nSteps = 0
  call system_clock(c0, cRate)

  do while (time < relaxTime - 0.5_R64 * timeStep)
    do innerStep = 1, nTimeSteps
      call Propagator_Propagate(Method_state, time, time + timeStep)
      time = time + timeStep
      nSteps = nSteps + 1
    end do
    call WriteProgress(time)

    if (abs(time) < 0.5_R64 * timeStep) then
      energy0 = Method_GetEnergy(0.0_R64)
      write (*, '(A,F22.12)') " ground state (t = 0): energy ", energy0
      if (saveQ) then
        call TwoRdm_FillRdm2FromTwoRdmState(rdm2, TwoRdm_twoRdmState)
        call Orbs_SaveOrbs(Orbs_orbs)
        call Method_Mb_OrbBased_SaveRdm2Blocks(rdm2)
      end if
    end if
  end do
  close (io)

  call system_clock(c1)
  write (*, '(A,F22.12,A,ES10.2)') " energy after the relaxation: ", Method_GetEnergy(time), &
    ",  drift ", Method_GetEnergy(time) - energy0
  write (*, '(A,F10.2,A,F8.4,A)') " wall time ", real(c1 - c0, R64) / cRate, " s, ", &
    real(c1 - c0, R64) / cRate / max(nSteps, 1), " s per time step"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  !> Restricted Hartree-Fock in the span of the active orbitals: iterates the
  !> Fock matrix f = h + sum_kl (2 <ik|jl> - <ik|lj>) rdm1(l, k) of the
  !> nBodies(1) lowest orbitals (start: the first ones), then rotates the
  !> orbitals and the stored h1 and h2 of the method to its eigenvectors.
  subroutine RotateToActiveSpaceHartreeFock
    complex(R64), allocatable :: h(:, :), w(:, :, :, :), rot(:, :), fock(:, :), rdm1(:, :), rdm1Old(:, :), toNew(:, :)
    real(R64), allocatable :: evals(:)
    real(R64) :: change
    integer(I32), parameter :: maxIterations = 200
    integer(I32) :: n, nOcc, i, j, k, l, iteration

    n = Method_Mb_OrbBased_nOrbs(1)
    nOcc = Method_Mb_nBodies(1)
    if (size(Orbs_orbs, 2) .ne. n) error stop "P_GroundTd2rdm needs restricted orbitals"

    h = Method_Mb_OrbBased_h1(1:n, 1:n)
    w = 2.0_R64 * Method_Mb_OrbBased_h2(1:n, n + 1:, 1:n, n + 1:)

    allocate (rdm1(n, n), fock(n, n))
    rdm1 = 0.0_R64
    do i = 1, nOcc
      rdm1(i, i) = 1.0_R64
    end do
    do iteration = 1, maxIterations
      fock = h
      do j = 1, n
        do i = 1, n
          do l = 1, n
            do k = 1, n
              fock(i, j) = fock(i, j) + (2.0_R64 * w(i, k, j, l) - w(i, k, l, j)) * rdm1(l, k)
            end do
          end do
        end do
      end do
      call LapackLib_DiagonalizeGeneric(evals, rot, fock, .true.)
      rdm1Old = rdm1
      rdm1 = matmul(rot(:, :nOcc), conjg(transpose(rot(:, :nOcc))))
      change = maxval(abs(rdm1 - rdm1Old))
      if (change < 1e-12_R64) exit
    end do
    if (.not. (change < 1e-12_R64)) error stop "P_GroundTd2rdm: the Hartree-Fock iteration in the active space did not converge"
    write (*, '(A,I0,A)') " Hartree-Fock in the active space: ", iteration, " iterations, orbital energies"
    write (*, '(5F14.6)') evals

    ! orbitals phi rot; matrices and integrals transform with rot^+ on the
    ! bra orbitals (RdmBasisChange_Transform: first half of the indices)
    toNew = conjg(transpose(rot))
    Orbs_orbs(:, :) = matmul(Orbs_orbs, rot)
    h = matmul(toNew, matmul(h, rot))
    call RdmBasisChange_Transform(w, 4, toNew)
    w = 0.5_R64 * w

    Method_Mb_OrbBased_h1(:, :) = 0.0_R64
    Method_Mb_OrbBased_h1(1:n, 1:n) = h
    Method_Mb_OrbBased_h1(n + 1:, n + 1:) = h
    Method_Mb_OrbBased_h2(:, :, :, :) = 0.0_R64
    Method_Mb_OrbBased_h2(1:n, 1:n, 1:n, 1:n) = w
    Method_Mb_OrbBased_h2(n + 1:, 1:n, n + 1:, 1:n) = w
    Method_Mb_OrbBased_h2(1:n, n + 1:, 1:n, n + 1:) = w
    Method_Mb_OrbBased_h2(n + 1:, n + 1:, n + 1:, n + 1:) = w

  end subroutine

  subroutine WriteProgress(t)
    real(R64), intent(in) :: t

    real(R64), allocatable :: natOcc(:)
    complex(R64), allocatable :: rdm1(:, :), natOrbs(:, :), rdm2(:, :, :, :), rdm2ab(:, :, :, :)
    real(R64) :: energy, minD, minQ
    integer(I32) :: n

    n = Method_Mb_OrbBased_nOrbs(1)
    call TwoRdm_FillRdm1FromTwoRdmState(rdm1, TwoRdm_twoRdmState)
    call TwoRdm_FillRdm2FromTwoRdmState(rdm2, TwoRdm_twoRdmState)
    call RdmDiagonalize_Rdm1(natOcc, natOrbs, rdm1(1:n, 1:n))
    call RdmTrace_FillRdm2abFromRdm2(rdm2ab, rdm2, n, n)
    call RdmHole_MinEigenvaluesDQab(minD, minQ, rdm2ab)
    energy = Method_GetEnergy(t)

    write (io, '(F10.4, F8.4, F24.14, 2ES13.4, *(ES15.7))') t, Method_Mb_OrbBased_SwitchFactor(t), energy, &
      minD, minQ, natOcc
    flush (io)
    if (abs(t - nint(t)) < 1e-9_R64) &
      write (*, '(A,F9.3,A,F7.4,A,F22.12,A,2ES11.3)') " t =", t, "  lambda =", Method_Mb_OrbBased_SwitchFactor(t), &
      "  energy =", energy, "  min D, Q:", minD, minQ

  end subroutine

end program
