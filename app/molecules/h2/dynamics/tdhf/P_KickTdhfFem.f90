!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Delta-kick linear response of a molecule with TDHF on the finite-element
!> grid (grid.fem), started from the orbitals of groundstate/scf/P_GroundScfFem. Run
!> commands in app/README.md. The orbitals of the SCF run and its femMesh.dat
!> are read from orbsInit.load.dir; every other file is read from and written
!> to the working directory.
!>
!> TDHF is restricted MCTDHX with nOrbs = nBodies (a single configuration,
!> coeffsInit.unary), started from the orbitals of P_GroundScfFem (orbsInit.load of
!> orbXX_YY.in). They must come from the IDENTICAL grid.fem mesh: the program
!> stops if the mesh differs from femMesh.dat of the SCF run (when present).
!>
!> Kick: at t = 0 every orbital is multiplied by exp(i κ r_a), κ =
!> program.kappa, a = program.kickDirection (1, 2, 3 for x, y, z). Every
!> electron receives the momentum κ along a, the response to the field
!> E(t) = −κ e_a δ(t); the density, hence ⟨r⟩, is unchanged at t = 0+. The
!> program prints the kick energy E(0+) − E(0) (= N κ²/2) and the initial
!> slope d⟨r_a⟩/dt(0+) = i Tr(ρ [T, r_a]) of the space-discretized equations
!> (PrinterObservableFem_DumpDipoleVelocity; = N κ, Thomas–Reiche–Kuhn on the
!> grid).
!>
!> Propagation (wiring of test/simulationTests/Fem/T_Fem_05 .. T_Fem_09):
!> propagator.splitStep.order2 with exactPartScheme strang; the field-free
!> one-body part h0 = T + V_nuc exactly by the unitary Cayley step
!> (method.mb.orbBased.exactFieldFreeQ, GMRES preconditioned by the FDM
!> inverse of 1 + i dt T/2), the mean field by RK4 (integrator 2; integrator 1
!> integrates only a sysGauge term, none for a kick), program.timeStep up to
!> program.endTime in blocks of program.outputStep. With
!> method.mb.orbBased.exactFieldFreeMeanFieldQ (set in the inputs) the Cayley
!> step propagates h0 + U and RK4 only the change of the mean field, W - U:
!> without it a heavy core (the O of H2O) needs tiny time steps
!> (src/Method/A_Method.md, "Frozen Mean-Field Reference"). U = J0 - K0 is the Fock mean field of the
!> ground state, built here before the kick
!> (Method_Mb_OrbBased_BuildFieldFreeMeanField; Method would otherwise build
!> it from the kicked state at the first step): the remainder W - U then has
!> no static part, and the energy drifts far less. The
!> ground state it is built from is saved
!> (meanFieldReference.dat) for continued runs. An absorber block
!> (absorber.fem.cosinus) is applied after every step; with
!> absorber.keepOrthonormalQ the orbitals are then re-orthonormalized and the
!> absorbed norm moves into the CI coefficient.
!>
!> Output, one line at t = startTime and after every program.outputStep:
!>   dipole.dat       time, <x>, <y>, <z> (PrinterObservableFem_DumpDipole:
!>                    position expectation summed over the electrons); the
!>                    header records κ, the direction and the electron number
!>                    for absorption_spectrum.py
!>   observables.dat  time, E_el, E_tot = E_el + E_nn, N = Tr(ρS)
!> Checkpoints: with program.saveStep > 0 (a multiple of program.outputStep)
!> the orbitals and the CI vector are saved (orbXX_YY.in, coeffs.in;
!> checkpoint.dat holds their time) after every saveStep and at the end. A run
!> continues from the last checkpoint in the same directory with
!> program.startTime set to that time (no kick, no orthonormalization; with an
!> absorber and keepOrthonormalQ also coeffsInit.load, which restores the
!> absorbed norm); the program stops unless checkpoint.dat holds startTime.
!> dipole.dat and observables.dat keep their lines up to startTime and are
!> continued; for that, orbsInit.load must read the checkpoint orbitals of the
!> run directory (no dir). With exactFieldFreeMeanFieldQ a continued run rebuilds U from
!> meanFieldReference.dat of the first run (the program stops if it is missing
!> or from another mesh or orbital space), so it reproduces the uninterrupted
!> run; rebuilding U from the checkpoint state instead would leave a static
!> O(kappa) part in W - U, and the energy would grow quadratically. program.endTime − program.startTime must be a
!> multiple of program.outputStep.
program P_KickTdhfFem
  use, intrinsic :: iso_fortran_env, only: output_unit
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Utils_Constants
  use M_Utils_PrinterObservableFem
  use M_Grid
  use M_Grid_Fem
  use M_SysKinetic
  use M_SysPotential
  use M_SysPotential_Fem_Coulomb, only: SysPotential_Fem_Coulomb_nuclearRepulsion
  use M_SysGauge
  use M_SysInteraction
  use M_Orbs
  use M_OrbsInit
  use M_ConfigList
  use M_Coeffs
  use M_CoeffsInit
  use M_Method
  use M_Method_Mb
  use M_Method_Mb_OrbBased
  use M_IntegratorList
  use M_Propagator
  use M_Absorber

  implicit none

  character(len=*), parameter :: DIPOLE_FILE = "dipole.dat"
  character(len=*), parameter :: OBSERVABLES_FILE = "observables.dat"
  character(len=*), parameter :: CHECKPOINT_FILE = "checkpoint.dat"
  character(len=:), allocatable :: meshFile
  character(len=*), parameter :: REFERENCE_FILE = "meanFieldReference.dat"
  character(len=1), parameter :: AXIS_NAMES(3) = ["x", "y", "z"]

  type(T_IntegratorList_FabricateInput) :: IntegratorListInput(2)
  complex(R64), allocatable :: rdm1(:, :)
  real(R64), allocatable :: coordinate(:)
  real(R64) :: kappa, startTime, endTime, timeStep, outputStep, saveStep, time, nElectrons, kickVector(3)
  real(R64) :: energyGround, energyKick, energyStart, slope, dipoleGround(3), velocity(3), seconds, secondsBlock
  integer(I64) :: nStepsDone, count0, countBlock
  integer(I32) :: kickDirection, nStepsPerOutput, nOutputsPerSave, nOutputsTotal, nOutputs, iStep, io
  logical :: gaugeQ, restartQ

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call Json_LoadJsonFile()
  gaugeQ = Json_GetExistence("sysGauge")

  call Grid_Fabricate
  call SysKinetic_Fabricate
  call SysPotential_Fabricate
  if (gaugeQ) call SysGauge_Fabricate
  call SysInteraction_Fabricate
  call Method_Fabricate
  call Orbs_Fabricate
  call OrbsInit_Fabricate
  call ConfigList_Fabricate
  call Coeffs_Fabricate
  call CoeffsInit_Fabricate
  IntegratorListInput(1) % TimeDerivative => Method_Mb_OrbBased_TimeDerivativeOrbsLin
  IntegratorListInput(2) % TimeDerivative => Method_Mb_OrbBased_TimeDerivativeCoeffsPlusOrbsNonLin
  call IntegratorList_Fabricate(IntegratorListInput)
  call Propagator_Fabricate()
  call Absorber_Fabricate

  call Say_Fabricate("program")
  kappa = Json_Get("program.kappa", 1e-3_R64)
  kickDirection = Json_Get("program.kickDirection", 3)
  startTime = Json_Get("program.startTime", 0.0_R64)
  endTime = Json_Get("program.endTime", 10.0_R64)
  timeStep = Json_Get("program.timeStep", 0.05_R64)
  outputStep = Json_Get("program.outputStep", 0.1_R64)
  saveStep = Json_Get("program.saveStep", 0.0_R64)

  if (kickDirection < 1 .or. kickDirection > 3) error stop "P_KickTdhfFem: program.kickDirection must be 1 (x), 2 (y) or 3 (z)"
  if (.not. (timeStep > 0.0_R64)) error stop "P_KickTdhfFem: program.timeStep must be positive"
  nStepsPerOutput = nint(outputStep / timeStep)
  if (nStepsPerOutput < 1 .or. abs(nStepsPerOutput * timeStep - outputStep) > 1e-9_R64 * outputStep) &
    error stop "P_KickTdhfFem: program.outputStep must be a positive multiple of program.timeStep"
  nOutputsPerSave = 0
  if (saveStep > 0.0_R64) then
    nOutputsPerSave = nint(saveStep / outputStep)
    if (nOutputsPerSave < 1 .or. abs(nOutputsPerSave * outputStep - saveStep) > 1e-9_R64 * saveStep) &
      error stop "P_KickTdhfFem: program.saveStep must be a multiple of program.outputStep"
  end if
  if (startTime < 0.0_R64) error stop "P_KickTdhfFem: program.startTime must not be negative"
  nOutputsTotal = nint((endTime - startTime) / outputStep)
  if (nOutputsTotal < 0 .or. &
      abs(nOutputsTotal * outputStep - (endTime - startTime)) > 1e-9_R64 * max(outputStep, abs(endTime))) &
    error stop "P_KickTdhfFem: program.endTime - program.startTime must be a non-negative multiple of program.outputStep"
  restartQ = startTime > 0.0_R64
  if (restartQ) call CheckCheckpointTime(CHECKPOINT_FILE, startTime)

  call Grid_Setup
  ! Before the orbitals are loaded (Method_Setup): LoadData only notices a
  ! file shorter than the mesh
  meshFile = "femMesh.dat"
  if (len_trim(Json_Get("orbsInit.load.dir", "")) > 0) meshFile = Json_Get("orbsInit.load.dir", "")//"/femMesh.dat"
  call CheckFemMeshFile(meshFile)
  call SysKinetic_Setup
  call SysPotential_Setup
  if (gaugeQ) call SysGauge_Setup
  call SysInteraction_Setup
  call Orbs_Setup
  call OrbsInit_Setup
  call ConfigList_Setup
  call Coeffs_Setup
  call CoeffsInit_Setup
  call Method_Setup
  call IntegratorList_Setup
  call Propagator_Setup
  call Absorber_Setup

  nElectrons = real(Method_Mb_nBodiesSum, kind=R64)
  ! r_a at the nodes, a = kickDirection
  kickVector(:) = 0.0_R64
  kickVector(kickDirection) = 1.0_R64
  call Grid_Fem_FillNodeCoordinate(coordinate, kickVector)

  !==========================================
  call Say_Section("kick")
  !==========================================

  write (*, '(1X, A, I0, A, 3(I0, 1X), A, I0, A, I0, A)') "mesh: N = ", Grid_nPoints, " = ", Grid_Fem_nx, Grid_Fem_ny, &
    Grid_Fem_nz, "(x y z); ", Method_Mb_nBodiesSum, " electrons in ", Orbs_nOrbsInState, " stored orbitals"

  if (.not. restartQ) then
    ! Loaded SCF orbitals: orthonormal up to the eigensolver tolerance
    call Orbs_Orthonormalize(Orbs_orbs)

    energyGround = Method_GetEnergy(0.0_R64)
    call Method_Mb_OrbBased_FillRdm1(rdm1, Method_state)
    call PrinterObservableFem_DumpDipole(rdm1, Orbs_orbs, "", .false., 0.0_R64, dipole_=dipoleGround)

    ! Frozen mean-field reference from the ground state (before the kick),
    ! saved for continued runs
    if (Method_Mb_OrbBased_exactFieldFreeMeanFieldQ) call BuildAndSaveReference(REFERENCE_FILE, Orbs_orbs, rdm1)

    call ApplyKick(Orbs_orbs, kappa, coordinate)

    energyKick = Method_GetEnergy(0.0_R64)
    call Method_Mb_OrbBased_FillRdm1(rdm1, Method_state)
    call PrinterObservableFem_DumpDipoleVelocity(rdm1, Orbs_orbs, "", .false., 0.0_R64, velocity_=velocity)
    slope = velocity(kickDirection)

    write (*, '(1X, A, ES12.4, 3A)') "kick phi <- exp(i kappa r_a) phi: kappa = ", kappa, ", direction ", &
      AXIS_NAMES(kickDirection), " (program.kickDirection)"
    write (*, '(1X, A, 3F18.12)') "ground state <r> (x, y, z):          ", dipoleGround
    write (*, '(1X, A, F22.14, A, F22.14)') "E_el before the kick ", energyGround, ", after ", energyKick
    if (abs(kappa) > 0.0_R64) then
      write (*, '(1X, A, ES16.8, A, ES16.8, A, ES10.2)') "kick energy E(0+) - E(0) = ", energyKick - energyGround, &
        ", N kappa^2/2 = ", 0.5_R64 * nElectrons * kappa**2, ", relative deviation ", &
        (energyKick - energyGround) / (0.5_R64 * nElectrons * kappa**2) - 1.0_R64
      write (*, '(1X, 3A, ES16.8, A, ES16.8, A, ES10.2)') "TRK on the grid: d<", AXIS_NAMES(kickDirection), &
        ">/dt(0+) = ", slope, ", N kappa = ", nElectrons * kappa, ", relative deviation ", &
        slope / (nElectrons * kappa) - 1.0_R64
    else
      write (*, '(1X, A)') "no kick (program.kappa = 0): field-free propagation of the ground state"
    end if

    ! Fresh output files; the dipole header carries the kick for the spectrum
    open (newunit=io, file=DIPOLE_FILE, status="replace", action="write")
    write (io, '(A, ES16.8, A, I0, 3A, I0)') "# P_KickTdhfFem: kick exp(i kappa r_a) at t = 0, kappa = ", kappa, &
      ", direction = ", kickDirection, " (", AXIS_NAMES(kickDirection), "), electrons = ", Method_Mb_nBodiesSum
    write (io, '(A)') "# time  <x>  <y>  <z>   (position expectation summed over the electrons, a.u.)"
    close (io)
    open (newunit=io, file=OBSERVABLES_FILE, status="replace", action="write")
    write (io, '(A)') "# time  E_el  E_tot = E_el + E_nn  N = Tr(rho S)"
    close (io)
  else
    write (*, '(1X, A, F12.4, A)') "continuing from the checkpoint at t = ", startTime, &
      " (orbXX_YY.in, coeffs.in): no kick"
    if (Method_Mb_OrbBased_exactFieldFreeMeanFieldQ) call LoadAndBuildReference(REFERENCE_FILE)
    call TruncateTimeSeries(DIPOLE_FILE, startTime - 0.5_R64 * timeStep)
    call TruncateTimeSeries(OBSERVABLES_FILE, startTime - 0.5_R64 * timeStep)
  end if

  !==========================================
  call Say_Section("propagate")
  !==========================================

  write (*, '(1X, A, F8.4, A, F8.4, A, F10.3, A, F10.3)') "dt = ", timeStep, ", output every ", outputStep, &
    ", from t = ", startTime, " to ", endTime
  call WriteObservables(energyStart, startTime, 0.0_R64, .true.)

  Method_Mb_OrbBased_cayleyNSolves = 0
  Method_Mb_OrbBased_cayleyNIterations = 0
  Method_Mb_OrbBased_cayleyMaxIterations = 0

  time = startTime
  nStepsDone = 0
  nOutputs = 0
  call system_clock(count0)
  do while (time < endTime - 0.5_R64 * timeStep)

    call system_clock(countBlock)
    do iStep = 1, nStepsPerOutput
      call Propagator_Propagate(Method_state, time, time + timeStep)
      nStepsDone = nStepsDone + 1
      time = startTime + real(nStepsDone, kind=R64) * timeStep
      call Absorber_ApplyAbsorber(Orbs_orbs)
      if (Absorber_keepOrthonormalQ) call Method_Mb_OrbBased_OrthonormalizeKeepingState(Method_state)
    end do
    secondsBlock = SecondsSince(countBlock)

    call WriteObservables(energyStart, time, secondsBlock / nStepsPerOutput, .false.)

    nOutputs = nOutputs + 1
    if (nOutputsPerSave > 0) then
      if (mod(nOutputs, nOutputsPerSave) .eq. 0) call SaveCheckpoint(time)
    end if

  end do
  seconds = SecondsSince(count0)

  !==========================================
  call Say_Section("printout results")
  !==========================================

  write (*, '(1X, A, I0, A, F12.2, A, F9.4, A, F10.2, A)') "propagation: ", nStepsDone, " steps, ", seconds, " s (", &
    seconds / max(nStepsDone, 1_I64), " s per step, ", seconds / max(time - startTime, timeStep), &
    " s per a.u., observables included)"
  write (*, '(1X, A, I0, A, I0, A, F6.2, A, I0)') "Cayley GMRES: ", Method_Mb_OrbBased_cayleyNSolves, " solves, ", &
    Method_Mb_OrbBased_cayleyNIterations, " iterations, mean ", &
    real(Method_Mb_OrbBased_cayleyNIterations, R64) / real(max(Method_Mb_OrbBased_cayleyNSolves, 1_I64), R64), &
    ", max ", Method_Mb_OrbBased_cayleyMaxIterations

  !==========================================
  call Say_Section("save data")
  !==========================================

  call SaveCheckpoint(time)

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> The kick phi <- exp(i kappa r_a) phi of every stored orbital, a one-body
  !> unitary that keeps the CI coefficients (orbs is the orbital part of
  !> Method_state).
  subroutine ApplyKick(orbs, kappa, coordinate)
    complex(R64), intent(inout), contiguous :: orbs(:, :)
    real(R64), intent(in) :: kappa
    real(R64), intent(in), contiguous :: coordinate(:)

    complex(R64), allocatable :: phase(:)
    integer(I32) :: iOrb

    phase = exp(IU * kappa * coordinate(:))
    do iOrb = 1, size(orbs, 2)
      orbs(:, iOrb) = phase(:) * orbs(:, iOrb)
    end do

  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Appends the dipole (dipole.dat) and E_el, E_tot, N (observables.dat) at
  !> time and prints them; energyStart is set on the first call (firstQ) and
  !> the energy change since then is printed.
  subroutine WriteObservables(energyStart, time, secondsPerStep, firstQ)
    real(R64), intent(inout) :: energyStart
    real(R64), intent(in) :: time
    real(R64), intent(in) :: secondsPerStep
    logical, intent(in) :: firstQ

    complex(R64), allocatable :: rdm1(:, :)
    real(R64) :: energy, norm, dipole(3)
    integer(I32) :: io

    energy = Method_GetEnergy(time)
    if (firstQ) energyStart = energy
    call Method_Mb_OrbBased_FillRdm1(rdm1, Method_state)
    call PrinterObservableFem_DumpDipole(rdm1, Orbs_orbs, DIPOLE_FILE, .false., time, dipole_=dipole)
    call PrinterObservableFem_DumpNorm(rdm1, Orbs_orbs, "", .false., time, norm_=norm)

    open (newunit=io, file=OBSERVABLES_FILE, status="unknown", position="append", action="write")
    write (io, '(*(ES24.15E3))') time, energy, energy + SysPotential_Fem_Coulomb_nuclearRepulsion, norm
    close (io)

    write (*, '(1X, A, F9.3, A, F20.12, A, ES10.2, A, F14.10, A, 3ES16.8, A, F8.3, A)') "t =", time, "  E_el ", energy, &
      "  dE ", energy - energyStart, "  N ", norm, "  <r> ", dipole, "  (", secondsPerStep, " s/step)"
    flush (output_unit)

  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Builds the frozen mean-field reference of the Cayley step
  !> (exactFieldFreeMeanFieldQ) from the stored orbitals orbs and the 1-RDM
  !> rdm1 of the ground state, and saves both to filename (unformatted
  !> stream: Grid_nPoints, Orbs_nOrbsInState, size(rdm1, 1), orbs, rdm1) for
  !> LoadAndBuildReference. Method then does not rebuild it at the first step.
  subroutine BuildAndSaveReference(filename, orbs, rdm1)
    character(len=*), intent(in) :: filename
    complex(R64), intent(in), contiguous :: orbs(:, :)
    complex(R64), intent(in), contiguous :: rdm1(:, :)

    integer(I32) :: io

    call Method_Mb_OrbBased_BuildFieldFreeMeanField(orbs, rdm1, 0.0_R64)
    open (newunit=io, file=filename, access="stream", form="unformatted", status="replace", action="write")
    write (io) Grid_nPoints, Orbs_nOrbsInState, int(size(rdm1, 1), I32)
    write (io) orbs, rdm1
    close (io)
    write (*, '(1X, 3A)') "frozen mean-field reference: built from the ground state before the kick, saved to ", &
      filename, " (for continued runs)"

  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Rebuilds the frozen mean-field reference of a continued run from the
  !> state saved by BuildAndSaveReference in the first run (the same U as
  !> there, so the continued run reproduces the uninterrupted one); stops if
  !> filename is missing, unreadable or from another mesh or orbital space.
  subroutine LoadAndBuildReference(filename)
    character(len=*), intent(in) :: filename

    complex(R64), allocatable :: orbsRef(:, :), rdm1Ref(:, :)
    integer(I32) :: io, ios, nPointsRef, nOrbsRef, nRdmRef
    logical :: existsQ

    inquire (file=filename, exist=existsQ)
    if (.not. existsQ) then
      write (*, '(1X, 3A)') "P_KickTdhfFem: continuing with exactFieldFreeMeanFieldQ needs ", filename, &
        " of the first run (the state the reference U is built from)"
      error stop "P_KickTdhfFem: the frozen mean-field reference of the first run is missing"
    end if

    open (newunit=io, file=filename, access="stream", form="unformatted", status="old", action="read")
    read (io, iostat=ios) nPointsRef, nOrbsRef, nRdmRef
    if (ios .ne. 0 .or. nPointsRef .ne. Grid_nPoints .or. nOrbsRef .ne. Orbs_nOrbsInState &
        .or. nRdmRef .ne. Method_Mb_OrbBased_nOrbsSum) then
      close (io)
      error stop "P_KickTdhfFem: the frozen mean-field reference file is unreadable or from another mesh or orbital space"
    end if
    allocate (orbsRef(nPointsRef, nOrbsRef), rdm1Ref(nRdmRef, nRdmRef))
    read (io, iostat=ios) orbsRef, rdm1Ref
    close (io)
    if (ios .ne. 0) error stop "P_KickTdhfFem: the frozen mean-field reference file is truncated"

    call Method_Mb_OrbBased_BuildFieldFreeMeanField(orbsRef, rdm1Ref, 0.0_R64)
    write (*, '(1X, 3A)') "frozen mean-field reference: rebuilt from ", filename, " (the state of the first run)"

  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Saves orbitals and CI vector (orbXX_YY.in, coeffs.in) and their time
  !> (checkpoint.dat).
  subroutine SaveCheckpoint(time)
    real(R64), intent(in) :: time

    integer(I32) :: io

    call Orbs_SaveOrbs(Orbs_orbs)
    call Coeffs_SaveCoeffs(Coeffs_coeffs)
    open (newunit=io, file=CHECKPOINT_FILE, status="replace", action="write")
    write (io, '(A)') "# time of the state in orbXX_YY.in and coeffs.in (restart: program.startTime)"
    write (io, '(F0.12)') time
    close (io)
    write (*, '(1X, A, F10.4)') "checkpoint saved at t =", time

  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Stops unless filename (written by SaveCheckpoint) exists and holds
  !> startTime: a continued run must start at the time of the state it loads.
  subroutine CheckCheckpointTime(filename, startTime)
    character(len=*), intent(in) :: filename
    real(R64), intent(in) :: startTime

    real(R64) :: checkpointTime
    integer(I32) :: io, ios
    logical :: existsQ

    inquire (file=filename, exist=existsQ)
    if (.not. existsQ) then
      write (*, '(1X, 3A)') "P_KickTdhfFem: program.startTime > 0 continues a run, but ", filename, " is missing"
      error stop "P_KickTdhfFem: no checkpoint to continue from (program.startTime > 0)"
    end if

    open (newunit=io, file=filename, status="old", action="read")
    read (io, *, iostat=ios)
    if (ios .eq. 0) read (io, *, iostat=ios) checkpointTime
    close (io)
    if (ios .ne. 0) error stop "P_KickTdhfFem: cannot read the time of the checkpoint"

    if (abs(checkpointTime - startTime) > 1e-6_R64 * max(1.0_R64, abs(startTime))) then
      write (*, '(1X, A, F0.6, 3A, F0.6)') "P_KickTdhfFem: program.startTime = ", startTime, ", but ", filename, &
        " holds t = ", checkpointTime
      error stop "P_KickTdhfFem: program.startTime must be the time of the checkpoint (orbXX_YY.in, coeffs.in)"
    end if
    write (*, '(1X, 3A, F0.6)') "checkpoint check: ", filename, " holds program.startTime = ", startTime

  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Stops unless grid.fem has the nLocals and element bounds written to
  !> filename by P_GroundScfFem (the mesh of the loaded orbitals); only notes the
  !> missing check if the file is absent.
  subroutine CheckFemMeshFile(filename)
    character(len=*), intent(in) :: filename

    real(R64), allocatable :: bounds(:)
    integer(I32) :: io, nLocals, nElements, iAxis
    logical :: existsQ

    inquire (file=filename, exist=existsQ)
    if (.not. existsQ) then
      write (*, '(1X, 3A)') "note: ", filename, " (P_GroundScfFem) not found, the mesh of the loaded orbitals is not checked"
      return
    end if

    open (newunit=io, file=filename, status="old", action="read")
    read (io, *)
    read (io, *)
    read (io, *) nLocals
    if (nLocals .ne. Grid_Fem_nLocals) call StopMeshMismatch("nLocals")
    do iAxis = 1, 3
      read (io, *) nElements
      if (nElements .ne. Grid_Fem_axes(iAxis) % nElements) call StopMeshMismatch("number of elements along "//AXIS_NAMES(iAxis))
      allocate (bounds(0:nElements))
      read (io, *) bounds
      if (maxval(abs(bounds - Grid_Fem_axes(iAxis) % bounds)) > 1e-10_R64 * max(1.0_R64, maxval(abs(bounds)))) &
        call StopMeshMismatch("element bounds along "//AXIS_NAMES(iAxis))
      deallocate (bounds)
    end do
    close (io)
    write (*, '(1X, 3A)') "mesh check: grid.fem is the mesh of the loaded orbitals (", filename, ")"

  end subroutine

  !> Error stop for a mesh that differs from the mesh of the loaded orbitals.
  subroutine StopMeshMismatch(what)
    character(len=*), intent(in) :: what

    write (*, '(1X, 4A)') "P_KickTdhfFem: the ", what, " of grid.fem differs from ", meshFile
    error stop "P_KickTdhfFem: the orbitals were computed on another mesh; use the grid block of the SCF input"

  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Keeps the lines of a time series file up to time tCut (first column);
  !> comment lines stay. Nothing happens if the file does not exist.
  subroutine TruncateTimeSeries(filename, tCut)
    character(len=*), intent(in) :: filename
    real(R64), intent(in) :: tCut

    character(len=4096), allocatable :: lines(:)
    character(len=4096) :: line
    real(R64) :: t
    integer(I32) :: io, ios, i
    logical :: existsQ

    inquire (file=filename, exist=existsQ)
    if (.not. existsQ) return

    allocate (lines(0))
    open (newunit=io, file=filename, status="old", action="read")
    do
      read (io, '(A)', iostat=ios) line
      if (ios .ne. 0) exit
      if (line(1:1) .ne. "#") then
        read (line, *) t
        if (t > tCut) exit
      end if
      lines = [lines, line]
    end do
    close (io)

    open (newunit=io, file=filename, status="replace", action="write")
    do i = 1, size(lines)
      write (io, '(A)') trim(lines(i))
    end do
    close (io)

  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Wall-clock seconds since count (system_clock counts).
  function SecondsSince(count) result(res)
    integer(I64), intent(in) :: count
    real(R64) :: res

    integer(I64) :: countNow, rate

    call system_clock(countNow, rate)
    res = real(countNow - count, kind=R64) / real(rate, kind=R64)

  end function

end program
