!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Field-free TDHF of H2 (R = 1.4) on the finite-element grid (grid.fem):
!> stationarity of the Hartree–Fock ground state under real-time propagation
!> with the Cayley exact part.
!>
!> TDHF is restricted MCTDHX with nOrbs = nBodies = [1, 1], started from the
!> RHF orbital of T_Fem_03 (orbsInit.load of orb01_01.in, same mesh, N =
!> 130331) and the single configuration (coeffsInit.unary). Propagator:
!> propagator.splitStep.order2 with method.mb.orbBased.exactFieldFreeQ, so
!> every linear half step is Cayley(h0, dt/4) · RK4[gauge](dt/2) ·
!> Cayley(h0, dt/4) (exactPartScheme strang; h0 = T + V_nuc, no sysGauge here,
!> so the RK4 substep is the identity) and the interaction part is RK4;
!> dt = program.timeStep = 0.01 up to program.endTime = 0.5.
!>
!> Checks:
!> 1. the loaded state reproduces the SCF energy of T_Fem_03 (E_el against the
!>    prototype value, 1e-8)
!> 2. the Cayley step itself, on the loaded orbital with dt =
!>    program.checkTimeStep (0.05): (1 + i dt h0/2) phi_new = (1 - i dt h0/2)
!>    phi to the GMRES tolerance (h0 applied independently here), unitarity
!>    (norm and <h0> unchanged), the group property Cayley(-dt) Cayley(dt)
!>    = 1 and the solver efficiency: at most MAX_ITERATIONS_CHECK GMRES
!>    iterations (7 measured; a wrong Krylov operator or initial residual is
!>    corrected by the true-residual restarts, so it shows only here)
!> 3. propagation: |E(t) - E(0)| < 1e-6, N(t) = Tr(rho S) = 2 to 1e-10,
!>    <x>, <y>, <z> = 0 (the mesh is symmetric under x, y, z -> -x, -y, -z; the
!>    loaded orbital carries ~5e-11 from the SCF convergence, so: |<r_a>| <
!>    1e-9 and unchanged to 1e-11) and the orbital stays the ground state,
!>    |<phi(0)|phi(t)>| = 1 to 1e-6 (its phase rotates); every Cayley solve
!>    of the propagation (dt/4 = 0.0025) takes at most
!>    MAX_ITERATIONS_PROPAGATION GMRES iterations (5 measured)
!> The dipole is written to dipoleFieldFree.dat (replaced at the start).
!> Prints the GMRES iteration counts of the Cayley solves and the time per
!> step.
program T_Fem_05_H2FieldFree
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Grid
  use M_Orbs
  use M_Method
  use M_Method_Mb_OrbBased
  use M_Propagator
  use testdrive, only: check, error_type

  implicit none

  !> E_el of the T_Fem_03 orbital (prototype on the same mesh)
  real(R64), parameter :: E_EL_PROTOTYPE = -1.8479538609_R64
  !> Number of electrons
  real(R64), parameter :: N_ELECTRONS = 2.0_R64
  !> Bounds of the GMRES iterations per Cayley solve (measured + 2): the
  !> contract check at dt = 0.05 (7) and the propagation at dt/4 = 0.0025 (5)
  integer(I32), parameter :: MAX_ITERATIONS_CHECK = 9
  integer(I32), parameter :: MAX_ITERATIONS_PROPAGATION = 7
  !> Dipole output of the run
  character(len=*), parameter :: DIPOLE_FILE = "dipoleFieldFree.dat"

  type(error_type), allocatable :: error
  complex(R64), allocatable :: phi0(:)
  real(R64) :: dt, endTime, checkTimeStep, time, energy0, energy, norm0, norm, dipole(3), dipole0(3)
  real(R64) :: maxEnergyDev, maxNormDev, maxDipole, maxDipoleDrift, minOverlap, overlap, seconds
  real(R64) :: cayleyResidual, cayleyNormDev, cayleyEnergyDev, cayleyRoundTrip
  integer(I64) :: count0
  integer(I32) :: nSteps, iStep, outputEvery, cayleyIterations, io

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call SetupFemDynamics("test/simulationTests/Fem/T_Fem_05_H2FieldFree.json")

  call Say_Fabricate("program")
  dt = Json_Get("program.timeStep", 0.01_R64)
  endTime = Json_Get("program.endTime", 0.5_R64)
  outputEvery = Json_Get("program.outputEvery", 10)
  checkTimeStep = Json_Get("program.checkTimeStep", 0.05_R64)
  nSteps = nint(endTime / dt)

  !==========================================
  call Say_Section("perform action")
  !==========================================

  phi0 = Orbs_orbs(:, 1)

  ! Fresh dipole file (no lines of earlier runs)
  open (newunit=io, file=DIPOLE_FILE, status="replace")
  close (io)
  call FemObservables(energy0, norm0, dipole0, 0.0_R64, DIPOLE_FILE)
  write (*, '(1X, A, I0, A, F18.13, A, F16.13, A, 3ES11.2)') "mesh N = ", Grid_nPoints, ", E_el(0) = ", energy0, &
    ", N(0) = ", norm0, ", dipole(0) = ", dipole0

  ! The Cayley step on the loaded orbital
  call CheckCayleyStep(cayleyResidual, cayleyNormDev, cayleyEnergyDev, cayleyRoundTrip, cayleyIterations, phi0, &
                       checkTimeStep)

  ! Propagation
  maxEnergyDev = 0.0_R64
  maxNormDev = abs(norm0 - N_ELECTRONS)
  maxDipole = maxval(abs(dipole0))
  maxDipoleDrift = 0.0_R64
  minOverlap = 1.0_R64
  Method_Mb_OrbBased_cayleyNSolves = 0
  Method_Mb_OrbBased_cayleyNIterations = 0
  Method_Mb_OrbBased_cayleyMaxIterations = 0

  time = 0.0_R64
  call system_clock(count0)
  do iStep = 1, nSteps
    call Propagator_Propagate(Method_state, time, time + dt)
    time = iStep * dt

    if (mod(iStep, outputEvery) .eq. 0 .or. iStep .eq. nSteps) then
      call FemObservables(energy, norm, dipole, time, DIPOLE_FILE)
      overlap = abs(Grid_InnerProduct(phi0, Orbs_orbs(:, 1)))
      maxEnergyDev = max(maxEnergyDev, abs(energy - energy0))
      maxNormDev = max(maxNormDev, abs(norm - N_ELECTRONS))
      maxDipole = max(maxDipole, maxval(abs(dipole)))
      maxDipoleDrift = max(maxDipoleDrift, maxval(abs(dipole - dipole0)))
      minOverlap = min(minOverlap, overlap)
      write (*, '(1X, A, F6.3, A, F18.13, A, ES10.2, A, ES10.2, A, 3ES10.2, A, ES10.2)') "t = ", time, &
        ": E_el ", energy, ", dE ", energy - energy0, ", N - 2 ", norm - N_ELECTRONS, ", dipole ", dipole, &
        ", 1 - |<phi0|phi>| ", 1.0_R64 - overlap
    end if
  end do
  seconds = ElapsedSeconds(count0)

  !==========================================
  call Say_Section("printout results")
  !==========================================

  print *
  print *, "E_el(0), prototype:         ", energy0, E_EL_PROTOTYPE, " difference ", energy0 - E_EL_PROTOTYPE
  print *, "Cayley step (dt = ", checkTimeStep, "):"
  print *, "  residual of (1+A) phi_new = (1-A) phi:  ", cayleyResidual
  print *, "  | |phi_new|^2 - |phi|^2 |:              ", cayleyNormDev
  print *, "  | <h0>_new - <h0> |:                    ", cayleyEnergyDev
  print *, "  | Cayley(-dt) Cayley(dt) phi - phi |:   ", cayleyRoundTrip
  print *, "  GMRES iterations, bound:              ", cayleyIterations, MAX_ITERATIONS_CHECK
  print *, "max |E(t) - E(0)|:          ", maxEnergyDev
  print *, "max |N(t) - 2|:             ", maxNormDev
  print *, "max |<r_a>|:                ", maxDipole
  print *, "max |<r_a>(t) - <r_a>(0)|:  ", maxDipoleDrift
  print *, "min |<phi(0)|phi(t)>|:      ", minOverlap
  write (*, '(1X, A, I0, A, I0, A, F6.2, A, I0, A, I0)') "Cayley GMRES during the propagation: ", &
    Method_Mb_OrbBased_cayleyNSolves, " solves, ", Method_Mb_OrbBased_cayleyNIterations, " iterations, mean ", &
    real(Method_Mb_OrbBased_cayleyNIterations, R64) / real(max(Method_Mb_OrbBased_cayleyNSolves, 1_I64), R64), &
    ", max ", Method_Mb_OrbBased_cayleyMaxIterations, ", bound ", MAX_ITERATIONS_PROPAGATION
  write (*, '(1X, A, I0, A, F8.2, A, F8.4, A)') "propagation: ", nSteps, " steps in ", seconds, " s, ", &
    seconds / nSteps, " s per step (observables every outputEvery steps included)"

  call check(error, energy0, E_EL_PROTOTYPE, thr=1e-8_R64)
  if (allocated(error)) error stop "T_Fem_05_H2FieldFree failure: loaded state is not the T_Fem_03 ground state"

  call check(error, cayleyResidual < 1e-11_R64)
  if (allocated(error)) error stop "T_Fem_05_H2FieldFree failure: Cayley step does not solve its linear system"

  call check(error, cayleyNormDev < 1e-11_R64 .and. cayleyEnergyDev < 1e-10_R64)
  if (allocated(error)) error stop "T_Fem_05_H2FieldFree failure: Cayley step not unitary"

  call check(error, cayleyRoundTrip < 1e-11_R64)
  if (allocated(error)) error stop "T_Fem_05_H2FieldFree failure: Cayley(-dt) is not the inverse of Cayley(dt)"

  call check(error, cayleyIterations <= MAX_ITERATIONS_CHECK)
  if (allocated(error)) error stop "T_Fem_05_H2FieldFree failure: too many GMRES iterations in the Cayley step at dt = 0.05"

  call check(error, maxEnergyDev < 1e-6_R64)
  if (allocated(error)) error stop "T_Fem_05_H2FieldFree failure: energy not conserved"

  call check(error, maxNormDev < 1e-10_R64)
  if (allocated(error)) error stop "T_Fem_05_H2FieldFree failure: particle number not conserved"

  call check(error, maxDipole < 1e-9_R64 .and. maxDipoleDrift < 1e-11_R64)
  if (allocated(error)) error stop "T_Fem_05_H2FieldFree failure: dipole of the symmetric molecule not zero or not constant"

  call check(error, 1.0_R64 - minOverlap < 1e-6_R64)
  if (allocated(error)) error stop "T_Fem_05_H2FieldFree failure: ground state not stationary"

  call check(error, Method_Mb_OrbBased_cayleyMaxIterations <= MAX_ITERATIONS_PROPAGATION)
  if (allocated(error)) error stop "T_Fem_05_H2FieldFree failure: too many GMRES iterations in the Cayley solves"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  !> One Cayley step of the exact field-free part on phi with dt and the
  !> checks of its contract: relative residual of (1 + A) phi_new =
  !> (1 - A) phi with A = i dt h0/2 (h0 by ApplyH0), changes of the norm and
  !> of <h0>, the relative distance of Cayley(-dt) phi_new to phi and the
  !> number of GMRES iterations of the forward step.
  subroutine CheckCayleyStep(residual, normDev, energyDev, roundTrip, nIterations, phi, dtCheck)
    use M_Utils_Constants

    real(R64), intent(out) :: residual
    real(R64), intent(out) :: normDev
    real(R64), intent(out) :: energyDev
    real(R64), intent(out) :: roundTrip
    integer(I32), intent(out) :: nIterations
    complex(R64), intent(in), contiguous :: phi(:)
    real(R64), intent(in) :: dtCheck

    complex(R64), allocatable :: work(:, :), h0Phi(:), h0New(:)
    integer(I64) :: countStart, iterationsBefore

    allocate (work(size(phi), 1), h0Phi(size(phi)), h0New(size(phi)))

    work(:, 1) = phi(:)
    iterationsBefore = Method_Mb_OrbBased_cayleyNIterations
    call system_clock(countStart)
    call Method_Mb_OrbBased_PropagateFieldFree(work, dtCheck, .false.)
    nIterations = int(Method_Mb_OrbBased_cayleyNIterations - iterationsBefore, I32)
    write (*, '(1X, A, F6.3, A, I0, A, F7.3, A)') "Cayley step dt = ", dtCheck, ": ", nIterations, &
      " GMRES iterations, ", ElapsedSeconds(countStart), " s"

    call ApplyH0(h0Phi, phi)
    call ApplyH0(h0New, work(:, 1))
    residual = MetricNorm(work(:, 1) + 0.5_R64 * IU * dtCheck * h0New - (phi - 0.5_R64 * IU * dtCheck * h0Phi)) &
               / MetricNorm(phi)
    normDev = abs(MetricNorm(work(:, 1))**2 - MetricNorm(phi)**2)
    energyDev = abs(real(Grid_InnerProduct(work(:, 1), h0New) - Grid_InnerProduct(phi, h0Phi), kind=R64))

    call Method_Mb_OrbBased_PropagateFieldFree(work, -dtCheck, .false.)
    roundTrip = MetricNorm(work(:, 1) - phi) / MetricNorm(phi)

  end subroutine

  !> h0phi = (T + V_static) phi for body type 1 (the generator of the exact
  !> field-free part).
  subroutine ApplyH0(h0phi, phi)
    use M_SysKinetic
    use M_SysPotential

    complex(R64), intent(out), contiguous :: h0phi(:)
    complex(R64), intent(in), contiguous :: phi(:)

    complex(R64), allocatable :: potential(:)

    call SysKinetic_MultiplyWithKineticOp(h0phi, phi, 0.0_R64, bt_=1)
    call SysPotential_FillExternalPotential(potential, 0.0_R64)
    h0phi(:) = h0phi(:) + potential(:) * phi(:)

  end subroutine

  !> Norm in the grid metric, sqrt(<f|f>).
  function MetricNorm(f) result(res)
    use M_Grid

    complex(R64), intent(in), contiguous :: f(:)
    real(R64) :: res

    res = sqrt(real(Grid_InnerProduct(f, f), kind=R64))

  end function

  include 'FemDynamics.inc'

end program
