!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Delta-kick linear response of a molecule with real-time TDDFT on the
!> finite-element grid (grid.fem), started from the Kohn–Sham orbitals of
!> groundstate/dft/P_GroundKsFem. Run commands in app/README.md. The orbitals
!> of the ground-state run and its femMesh.dat are read from orbsInit.load.dir;
!> every other file is written to the working directory.
!>
!> TDDFT is method.mb.orbBased.tdks (adiabatic Kohn–Sham potential of the
!> functional of sysXc, the same block as the ground-state run) with restricted
!> orbitals, started from the Kohn–Sham orbitals of P_GroundKsFem (orbsInit.load
!> of orbXX_YY.in). They must come from the IDENTICAL grid.fem mesh: the
!> program stops if the mesh differs from femMesh.dat of the ground-state run
!> (when present).
!>
!> Kick: at t = 0 every orbital is multiplied by exp(i κ r_a), κ =
!> program.kappa, a = program.kickDirection (1, 2, 3 for x, y, z): the
!> response to the field E(t) = −κ e_a δ(t); the density, hence ⟨r⟩, E_H and
!> E_xc, are unchanged at t = 0+. The program prints the kick energy
!> E(0+) − E(0) (= N κ²/2) and the initial slope d⟨r_a⟩/dt(0+) = N κ
!> (Thomas–Reiche–Kuhn on the grid).
!>
!> Propagation (wiring of test/simulationTests/FemDft/T_FemDft_03):
!> propagator.splitStep.order2 with exactPartScheme strang; the field-free
!> one-body part by the unitary Cayley step (fieldFree.cayley)
!> and, with fieldFree.cayley.meanFieldQ (set in the inputs),
!> h0 + U with the frozen Kohn–Sham reference U = v_KS[ρ0] − a K0 of the
!> ground state, built here before the kick
!> (FieldFree_BuildMeanField); RK4 (integrator 2) then only
!> propagates the change of the Kohn–Sham potential, v_KS(t) − a K(t) − U,
!> which is O(κ). Without the reference the splitting error of a heavy core
!> limits the time step (src/FieldFree/A_FieldFree.md, "Frozen Mean-Field
!> Reference"). The inputs use LDA: a GGA (e.g. PBE) makes this split step
!> stable only with a ~10x smaller time step on grid.fem (H2 PBE: 0.005;
!> src/Method/A_Method.md, "Time step with a GGA"). program.timeStep up to
!> program.endTime in blocks of program.outputStep.
!>
!> Output, one line at t = 0 and after every program.outputStep:
!>   dipole.dat       time, <x>, <y>, <z> (position expectation summed over
!>                    the electrons); the header records κ, the direction and
!>                    the electron number for ../../analysis/absorption_spectrum.py
!>   observables.dat  time, E_el, E_tot = E_el + E_nn, N = Tr(ρS)
program P_KickTddftFem
  use M_FieldFree
  use M_FieldFree_Cayley
  use, intrinsic :: iso_fortran_env, only: output_unit
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Printer_Observable_Fem
  use M_Grid
  use M_Grid_Fem
  use M_SysKinetic
  use M_SysPotential
  use M_SysPotential_Fem_Coulomb, only: SysPotential_Fem_Coulomb_nuclearRepulsion
  use M_SysInteraction
  use M_SysXc
  use M_Orbs
  use M_OrbsInit
  use M_OrbsInit_Load, only: OrbsInit_Load_dir
  use M_ConfigList
  use M_Coeffs
  use M_CoeffsInit
  use M_Method
  use M_Method_Mb
  use M_Method_Mb_OrbBased
  use M_IntegratorList
  use M_Propagator

  implicit none

  character(len=*), parameter :: DIPOLE_FILE = "dipole.dat"
  character(len=*), parameter :: OBSERVABLES_FILE = "observables.dat"
  character(len=1), parameter :: AXIS_NAMES(3) = ["x", "y", "z"]

  type(T_IntegratorList_FabricateInput) :: IntegratorListInput(2)
  complex(R64), allocatable :: rdm1(:, :)
  real(R64) :: kappa, timeStep, outputStep, endTime, time, energyGround, energyKick, energyStart, nElectrons
  real(R64) :: velocity(3), direction(3), secondsPerStep
  integer(I64) :: count0
  integer(I32) :: kickDirection, nStepsPerOutput, nOutputs, iOut, iStep, io
  character(len=:), allocatable :: meshFile

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call Json_LoadJsonFile()

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
  IntegratorListInput(1) % TimeDerivative => Method_Mb_OrbBased_TimeDerivativeOrbsLin
  IntegratorListInput(2) % TimeDerivative => Method_Mb_OrbBased_TimeDerivativeCoeffsPlusOrbsNonLin
  call IntegratorList_Fabricate(IntegratorListInput)
  call Propagator_Fabricate()

  call Say_Fabricate("program")
  kappa = Json_Get("program.kappa", 1e-3_R64)
  kickDirection = Json_Get("program.kickDirection", 3)
  timeStep = Json_Get("program.timeStep", 0.05_R64)
  outputStep = Json_Get("program.outputStep", 0.1_R64)
  endTime = Json_Get("program.endTime", 10.0_R64)
  if (.not. Json_GetExistence("method.mb.orbBased.tdks")) error stop "P_KickTddftFem: method.mb.orbBased.tdks is required"
  if (kickDirection < 1 .or. kickDirection > 3) error stop "P_KickTddftFem: program.kickDirection must be 1 (x), 2 (y) or 3 (z)"
  if (.not. (timeStep > 0.0_R64)) error stop "P_KickTddftFem: program.timeStep must be positive"
  nStepsPerOutput = nint(outputStep / timeStep)
  if (nStepsPerOutput < 1 .or. abs(nStepsPerOutput * timeStep - outputStep) > 1e-9_R64 * outputStep) &
    error stop "P_KickTddftFem: program.outputStep must be a positive multiple of program.timeStep"
  nOutputs = nint(endTime / outputStep)
  if (nOutputs < 0 .or. abs(nOutputs * outputStep - endTime) > 1e-9_R64 * max(endTime, outputStep)) &
    error stop "P_KickTddftFem: program.endTime must be a non-negative multiple of program.outputStep"

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
  call IntegratorList_Setup
  call Propagator_Setup

  ! The loaded orbitals must live on this mesh
  meshFile = "femMesh.dat"
  if (len_trim(OrbsInit_Load_dir) > 0) meshFile = trim(OrbsInit_Load_dir)//"/femMesh.dat"
  call CheckFemMesh(meshFile)
  call Orbs_Orthonormalize(Orbs_orbs)

  !==========================================
  call Say_Section("kick")
  !==========================================

  call Method_Mb_OrbBased_FillRdm1(rdm1, Method_state)
  call Printer_Observable_Fem_DumpNorm(rdm1, Orbs_orbs, "", .false., 0.0_R64, norm_=nElectrons)
  energyGround = Method_GetEnergy(0.0_R64)

  ! The frozen Kohn-Sham reference of the ground state (for semilocal
  ! functionals equal to that of the kicked state, whose density is the same)
  if (FieldFree_meanFieldQ) &
    call FieldFree_BuildMeanField(Orbs_orbs, rdm1, 0.0_R64)

  direction(:) = 0.0_R64
  direction(kickDirection) = 1.0_R64
  call KickOrbitals(Orbs_orbs, kappa, direction)
  energyKick = Method_GetEnergy(0.0_R64)
  call Printer_Observable_Fem_DumpDipoleVelocity(rdm1, Orbs_orbs, "", .false., 0.0_R64, velocity_=velocity)

  write (*, '(1X, A, F20.12, A, F20.12)') "E_el(0) ", energyGround, ", E_el(0+) ", energyKick
  if (abs(kappa) > 0.0_R64) then
    write (*, '(1X, A, ES16.8, A, ES16.8, A, ES10.2)') "kick energy E(0+) - E(0) = ", energyKick - energyGround, &
      ", N kappa^2/2 = ", 0.5_R64 * nElectrons * kappa**2, ", relative deviation ", &
      (energyKick - energyGround) / (0.5_R64 * nElectrons * kappa**2) - 1.0_R64
    write (*, '(1X, 3A, ES16.8, A, ES16.8, A, ES10.2)') "TRK on the grid: d<", AXIS_NAMES(kickDirection), &
      ">/dt(0+) = ", velocity(kickDirection), ", N kappa = ", nElectrons * kappa, ", relative deviation ", &
      velocity(kickDirection) / (nElectrons * kappa) - 1.0_R64
  else
    write (*, '(1X, A)') "no kick (program.kappa = 0): field-free propagation of the ground state"
  end if

  ! Fresh output files; the dipole header carries the kick for the spectrum
  open (newunit=io, file=DIPOLE_FILE, status="replace", action="write")
  write (io, '(A, ES16.8, A, I0, 3A, I0)') "# P_KickTddftFem: kick exp(i kappa r_a) at t = 0, kappa = ", kappa, &
    ", direction = ", kickDirection, " (", AXIS_NAMES(kickDirection), "), electrons = ", Method_Mb_nBodiesSum
  write (io, '(A)') "# time  <x>  <y>  <z>   (position expectation summed over the electrons, a.u.)"
  close (io)
  open (newunit=io, file=OBSERVABLES_FILE, status="replace", action="write")
  write (io, '(A)') "# time  E_el  E_tot = E_el + E_nn  N = Tr(rho S)"
  close (io)

  !==========================================
  call Say_Section("propagation")
  !==========================================

  time = 0.0_R64
  call WriteObservables(energyStart, time, 0.0_R64, .true.)
  do iOut = 1, nOutputs
    call system_clock(count0)
    do iStep = 1, nStepsPerOutput
      call Propagator_Propagate(Method_state, time, time + timeStep)
      time = time + timeStep
    end do
    secondsPerStep = SecondsSince(count0) / nStepsPerOutput
    call WriteObservables(energyStart, time, secondsPerStep, .false.)
  end do

  write (*, *)
  write (*, '(1X, A, I0, A, I0, A, F8.2)') "Cayley GMRES: ", FieldFree_Cayley_nSolves, " solves, mean ", &
    nint(real(FieldFree_Cayley_nIterations, R64) / real(max(FieldFree_Cayley_nSolves, 1_I64), R64)), &
    " iterations; last block ", secondsPerStep
  write (*, '(1X, A)') "spectrum: python3 ../../analysis/absorption_spectrum.py dipole.dat"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> orbs(:, i) = exp(i kappa direction · r) orbs(:, i).
  subroutine KickOrbitals(orbs, kappa, direction)
    complex(R64), intent(inout), contiguous :: orbs(:, :)
    real(R64), intent(in) :: kappa
    real(R64), intent(in) :: direction(3)

    real(R64), allocatable :: coordinate(:)
    integer(I32) :: iOrb

    call Grid_Fem_FillNodeCoordinate(coordinate, direction)
    do iOrb = 1, size(orbs, 2)
      orbs(:, iOrb) = exp(cmplx(0.0_R64, kappa * coordinate(:), kind=R64)) * orbs(:, iOrb)
    end do

  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Appends the dipole (dipole.dat) and E_el, E_tot, N (observables.dat) at
  !> time and prints them; energyStart is set on the first call (firstQ).
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
    call Printer_Observable_Fem_DumpDipole(rdm1, Orbs_orbs, DIPOLE_FILE, .false., time, dipole_=dipole)
    call Printer_Observable_Fem_DumpNorm(rdm1, Orbs_orbs, "", .false., time, norm_=norm)

    open (newunit=io, file=OBSERVABLES_FILE, status="unknown", position="append", action="write")
    write (io, '(*(ES24.15E3))') time, energy, energy + SysPotential_Fem_Coulomb_nuclearRepulsion, norm
    close (io)

    write (*, '(1X, A, F9.3, A, F20.12, A, ES10.2, A, F14.10, A, 3ES16.8, A, F8.3, A)') "t =", time, "  E_el ", energy, &
      "  dE ", energy - energyStart, "  N ", norm, "  <r> ", dipole, "  (", secondsPerStep, " s/step)"
    flush (output_unit)

  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Stops unless grid.fem is the mesh of femMesh.dat (P_GroundKsFem); a note
  !> if the file does not exist.
  subroutine CheckFemMesh(filename)
    character(len=*), intent(in) :: filename

    real(R64), allocatable :: bounds(:)
    integer(I32) :: io, nLocals, nElements, iAxis
    logical :: existsQ

    inquire (file=filename, exist=existsQ)
    if (.not. existsQ) then
      write (*, '(1X, 3A)') "note: ", filename, " (P_GroundKsFem) not found, the mesh of the loaded orbitals is not checked"
      return
    end if

    open (newunit=io, file=filename, status="old", action="read")
    read (io, *)
    read (io, *)
    read (io, *) nLocals
    if (nLocals .ne. Grid_Fem_nLocals) error stop "P_KickTddftFem: nLocals differs from the mesh of the loaded orbitals"
    do iAxis = 1, 3
      read (io, *) nElements
      if (nElements .ne. Grid_Fem_axes(iAxis) % nElements) &
        error stop "P_KickTddftFem: the number of elements differs from the mesh of the loaded orbitals"
      allocate (bounds(0:nElements))
      read (io, *) bounds
      if (maxval(abs(bounds - Grid_Fem_axes(iAxis) % bounds)) > 1e-10_R64 * max(1.0_R64, maxval(abs(bounds)))) &
        error stop "P_KickTddftFem: the element bounds differ from the mesh of the loaded orbitals"
      deallocate (bounds)
    end do
    close (io)
    write (*, '(1X, 3A)') "mesh check: grid.fem is the mesh of the loaded orbitals (", filename, ")"

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
