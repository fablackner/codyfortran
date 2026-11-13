!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Ethylene in a fixed orbital basis on grid.fem with TD-2RDM
!> (td2rdmStaticOrbs): correlated ground state, dipole kick, free evolution;
!> the absorption spectrum from the dipole.
!>
!> Basis: an active space of ethylene in a finite orbital basis (HF occupied
!> orbitals of ../../groundstate/scf plus correlating virtual orbitals,
!> make_inputs.py); the C 1s core is frozen and folded into the one-body
!> operator. The Hamiltonian is that of program.hamiltonianFile (stream,
!> written by make_inputs.py from the integrals of ../../analysis/integrals):
!>   int32 n, real64 core energy, complex128 h1eff(n, n), h2(n, n, n, n),
!>   dip(n, n, 3)
!> with h2 = <ij|kl>/2 (the convention of Method_Mb_OrbBased_FillH2). They
!> replace h1 and h2 of the setup in every spin block; the orbitals that
!> orbsInit loads only give the method its shape (n orbitals on grid.fem).
!>
!> Initial state (twoRdmInit):
!>   excited  the Hartree-Fock determinant (the closures are exact for it),
!>            then the correlation switched on over -T < t < 0
!>            (method.mb.orbBased.adiabaticSwitching.switchTime T):
!>            H = h1eff + (1 - lambda) U0 + lambda h2, U0 the mean field of the
!>            determinant in h2, so that at lambda = 0 it is an eigenstate of
!>            its Fock operator. No CI anywhere.
!>   load     the 2-RDM of twoRdmInit.load (an exact ground state).
!> Then the kick exp(-i program.kick D_a) (a = program.polarization, 1..3 =
!> x, y, z) and program.endTime of field-free propagation.
!>
!> Output <program.outputName>.dat: t (from the kick; negative during the
!> switching), switching factor lambda, active electronic energy (with the
!> core: + core energy), dipole <D_a> of the active electrons, min eigenvalue
!> of D and Q (up-down block), natural occupations (spin up). Wall time per
!> derivative on stdout.
program P_KickTd2rdmFem
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_Constants
  use M_Utils_LapackLib
  use M_Utils_RdmDiagonalize
  use M_Utils_RdmHole
  use M_Utils_RdmObservables
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
  use M_Purification
  use M_Method
  use M_Method_Mb
  use M_Method_Mb_OrbBased
  use M_IntegratorList
  use M_Propagator

  implicit none

  real(R64)    :: time, endTime, outputStep, timeStep, kick, coreEnergy
  integer(I32) :: nTimeSteps, innerStep, io, n, nOcc, polarization, nDerivatives
  integer(I64) :: c0, c1, cRate
  logical      :: purifyQ
  character(len=:), allocatable :: outputName, hamiltonianFile
  complex(R64), allocatable :: h1eff(:, :), h2eff(:, :, :, :), dip(:, :, :), dMat(:, :), h1Full(:, :), h2Full(:, :, :, :)
  type(T_IntegratorList_FabricateInput) :: IntegratorListInput(1)

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call Json_LoadJsonFile()
  purifyQ = Json_GetExistence("purification")

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
  if (purifyQ) call Purification_Fabricate
  IntegratorListInput(1) % TimeDerivative => Method_TimeDerivative
  call IntegratorList_Fabricate(IntegratorListInput)
  call Propagator_Fabricate()

  call Say_Fabricate("program")
  hamiltonianFile = Json_Get("program.hamiltonianFile", "hamiltonian.bin")
  kick = Json_Get("program.kick", 0.01_R64)
  polarization = Json_Get("program.polarization", 3)
  endTime = Json_Get("program.endTime", 80.0_R64)
  outputStep = Json_Get("program.outputStep", 0.2_R64)
  nTimeSteps = Json_Get("program.nTimeSteps", 1)
  outputName = Json_Get("program.outputName", "td2rdm")

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
  if (purifyQ) call Purification_Setup
  call IntegratorList_Setup
  call Propagator_Setup
  if (Json_GetExistence("sysGauge") .or. .not. (SysKinetic_timeIndependentQ .and. SysPotential_timeIndependentQ)) &
    error stop "P_KickTd2rdmFem: program.hamiltonianFile needs a static h1 (no sysGauge, static sysPotential)"

  !==========================================
  call Say_Section("Hamiltonian of the active space")
  !==========================================

  n = Method_Mb_OrbBased_nOrbs(1)
  nOcc = Method_Mb_nBodies(1)
  call ReadHamiltonian(hamiltonianFile)

  ! h1 = h1eff in both spin blocks, h2 in the four spin blocks (s1 s2 s1 s2)
  Method_Mb_OrbBased_h1(:, :) = 0.0_R64
  Method_Mb_OrbBased_h1(1:n, 1:n) = h1eff
  Method_Mb_OrbBased_h1(n + 1:, n + 1:) = h1eff
  Method_Mb_OrbBased_h2(:, :, :, :) = 0.0_R64
  Method_Mb_OrbBased_h2(1:n, 1:n, 1:n, 1:n) = h2eff
  Method_Mb_OrbBased_h2(n + 1:, 1:n, n + 1:, 1:n) = h2eff
  Method_Mb_OrbBased_h2(1:n, n + 1:, 1:n, n + 1:) = h2eff
  Method_Mb_OrbBased_h2(n + 1:, n + 1:, n + 1:, n + 1:) = h2eff
  h1Full = Method_Mb_OrbBased_h1
  h2Full = Method_Mb_OrbBased_h2
  ! the switching reference of this h2 (the setup built it from the grid h2)
  call Method_Mb_OrbBased_SetSwitchingReference
  allocate (dMat(2 * n, 2 * n))
  dMat = 0.0_R64
  dMat(1:n, 1:n) = dip(:, :, polarization)
  dMat(n + 1:, n + 1:) = dip(:, :, polarization)
  write (*, '(A,I0,A,I0,A,F18.10)') " active orbitals ", n, ", electrons per spin ", nOcc, ", core energy ", coreEnergy

  open (newunit=io, file=outputName//".dat", status="replace")
  nDerivatives = 0
  timeStep = outputStep / nTimeSteps
  call system_clock(c0, cRate)

  if (Json_GetExistence("twoRdmInit.excited")) then
    if (.not. (Method_Mb_OrbBased_switchTime > 0.0_R64)) &
      error stop "P_KickTd2rdmFem: the determinant start needs method.mb.orbBased.adiabaticSwitching"
    call AdiabaticSwitching
  end if

  !==========================================
  call Say_Section("kick and free evolution")
  !==========================================

  write (*, '(A,F22.12)') " ground state (before the kick): active energy ", ActiveEnergy()
  call KickTwoRdm(kick)

  time = 0.0_R64
  call WriteObservables(time, 1.0_R64)
  do while (time < endTime - 0.5_R64 * timeStep)
    do innerStep = 1, nTimeSteps
      call Advance(time)
    end do
    call WriteObservables(time, 1.0_R64)
  end do
  close (io)

  call system_clock(c1)
  write (*, '(A,I0,A,F10.4,A)') " active orbitals: ", n, ",  wall time per time derivative: ", &
    real(c1 - c0, R64) / cRate / max(nDerivatives, 1), " s"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  subroutine ReadHamiltonian(filename)
    character(len=*), intent(in) :: filename

    integer(I32) :: nFile, unit

    open (newunit=unit, file=filename, access="stream", form="unformatted", status="old")
    read (unit) nFile, coreEnergy
    if (nFile .ne. n) error stop "P_KickTd2rdmFem: program.hamiltonianFile has another number of orbitals"
    allocate (h1eff(n, n), h2eff(n, n, n, n), dip(n, n, 3))
    read (unit) h1eff
    read (unit) h2eff
    read (unit) dip
    close (unit)

  end subroutine

  !> From the HF determinant of twoRdmInit.excited, H(lambda) switched on
  !> over -T < t < 0.
  subroutine AdiabaticSwitching
    real(R64) :: tSwitch
    integer(I32) :: step

    tSwitch = -Method_Mb_OrbBased_switchTime
    write (*, '(A,F22.12)') " Hartree-Fock determinant: active energy ", ActiveEnergy()
    call WriteObservables(tSwitch, 0.0_R64)
    do while (tSwitch < -0.5_R64 * timeStep)
      do step = 1, nTimeSteps
        call Advance(tSwitch)
      end do
      ! the row at t = 0 follows the kick
      if (tSwitch < -0.5_R64 * timeStep) call WriteObservables(tSwitch, Method_Mb_OrbBased_SwitchFactor(tSwitch))
    end do

  end subroutine

  !> One time step from t (four derivatives: rk.o4Expl, make_inputs.py), then
  !> the purification.
  subroutine Advance(t)
    real(R64), intent(inout) :: t

    call Propagator_Propagate(Method_state, t, t + timeStep)
    t = t + timeStep
    nDerivatives = nDerivatives + 4
    if (purifyQ) call Purification_PurifyTwoRdmState(TwoRdm_twoRdmState, Method_Mb_OrbBased_h2)

  end subroutine

  !> Energy of the full active Hamiltonian.
  real(R64) function ActiveEnergy()
    complex(R64), allocatable :: rdm1(:, :), rdm2(:, :, :, :)

    call Method_Mb_OrbBased_FillRdm1(rdm1, Method_state)
    call Method_Mb_OrbBased_FillRdm2(rdm2, Method_state)
    ActiveEnergy = RdmObservables_Energy(rdm1, h1Full, rdm2, h2Full)

  end function

  !> rdm2 of exp(-i kick D)|Psi> (ket-first: m = exp(-i kick d) on the ket,
  !> conjg(m) on the bra indices).
  subroutine KickTwoRdm(kick)
    use M_Utils_RdmBasisChange

    real(R64), intent(in) :: kick

    complex(R64), allocatable :: rdm2(:, :, :, :), m(:, :), evecs(:, :)
    real(R64), allocatable :: evals(:)
    integer(I32) :: k

    call LapackLib_DiagonalizeGeneric(evals, evecs, 0.5_R64 * (dMat + conjg(transpose(dMat))), .true.)
    allocate (m, mold=evecs)
    do k = 1, size(evals)
      m(:, k) = evecs(:, k) * exp(cmplx(0.0_R64, -kick * evals(k), R64))
    end do
    m = matmul(m, conjg(transpose(evecs)))

    call TwoRdm_FillRdm2FromTwoRdmState(rdm2, TwoRdm_twoRdmState)
    call RdmBasisChange_Transform(rdm2, 4, m)
    call TwoRdm_FillTwoRdmStateFromRdm2(TwoRdm_twoRdmState, rdm2)

  end subroutine

  subroutine WriteObservables(t, sw)
    real(R64), intent(in) :: t
    real(R64), intent(in) :: sw

    complex(R64), allocatable :: rdm1(:, :), rdm2(:, :, :, :), rdm2ab(:, :, :, :), natOrbs(:, :)
    real(R64), allocatable :: natOcc(:)
    real(R64) :: minD, minQ

    call Method_Mb_OrbBased_FillRdm1(rdm1, Method_state)
    call Method_Mb_OrbBased_FillRdm2(rdm2, Method_state)
    call RdmDiagonalize_Rdm1(natOcc, natOrbs, rdm1(1:n, 1:n))
    call RdmTrace_FillRdm2abFromRdm2(rdm2ab, rdm2, n, n)
    call RdmHole_MinEigenvaluesDQab(minD, minQ, rdm2ab)

    write (io, '(F10.4, F8.4, F22.12, F18.10, 2ES13.4, *(ES15.7))') t, sw, &
      RdmObservables_Energy(rdm1, h1Full, rdm2, h2Full), real(sum(dMat * transpose(rdm1)), R64), minD, minQ, natOcc
    flush (io)

  end subroutine

end program
