!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> H2 restricted Hartree–Fock ground state (R = 1.4) on the finite-element grid
!> (grid.fem) with groundSolver.scf.stdImpl (method tdhx, restricted orbitals
!> [1, 1], potential mixing with α = 1, i.e. plain SCF) and the preconditioned
!> eigensolver diagonalizerList.lobpcg: ApplyPrecOnVec =>
!> GroundSolver_PreconditionAction ((T + max(−θ, σ_min))⁻¹ by the exact FDM
!> solve SysKinetic_SolveShifted), metricWeights = Grid_Fem_weights (wiring
!> and SCF loop in FemScf.inc, shared with T_Fem_04).
!>
!> Mesh (identical to the numpy prototype of the FEM design). This test writes
!> orb01_01.in (the converged 1σg orbital, Orbs_SaveOrbs); the dynamics tests
!> of this folder load it and must use exactly this grid block:
!>   "grid": { "fem": { "nLocals": 5,
!>     "xmin": -10.0, "xmax": 10.0, "ymin": -10.0, "ymax": 10.0, "zmin": -10.0, "zmax": 10.0,
!>     "centers": { "x": [0.0, 0.0], "y": [0.0, 0.0], "z": [-0.7, 0.7], "hMin": [0.3, 0.3] },
!>     "growth": 1.7, "hMax": 3.0, "multipoleLmax": 6, "multipoleCenter": [0.0, 0.0, 0.0] } }
!>   → n = (47, 47, 59), N = 130331
!>
!> Two SCF runs:
!>   1. from the orthonormalized gerade 1s(A) + 1s(B) start, exchange through
!>      the adaptively compressed exchange operator
!>      (groundSolver.scf.stdImpl.aceQ = true, as in the JSON)
!>   2. exact exchange (aceQ switched to false in the loaded JSON, ground
!>      solver fabricated again; its work arrays from Setup stay valid),
!>      continued from the converged ACE orbitals and LOBPCG block: a
!>      different fixed point of ACE would move the orbitals and the energy.
!>      (A fresh exact-exchange start reaches the same energy to 5e-14, but
!>      LOBPCG stalls at relative residuals 7e-8 … 3e-7 in every SCF step:
!>      on grid.fem the exact exchange operator is Hermitian only up to the
!>      multipole truncation of the Poisson boundary values. Measured before
!>      the stagnation exit: 1853 Fock products, 116 s on one thread; with it
!>      each solve stops 10 iterations into its floor, unconverged. ACE is
!>      Hermitian by construction.)
!>      Run 2 prints one expected warning, the multipole warning of
!>      Grid_Fem_SolvePoisson: the exact exchange solves the pair densities
!>      of the occupied orbital with the guard vectors, box-filling virtual
!>      states, of which a fraction 2.6e-4 > 1e-4 lies outside the convergence
!>      sphere (the bound density of run 1 stays at 4.9e-7). Its LOBPCG solve
!>      converges, so scf.stdImpl prints no aceQ hint (only an unconverged
!>      Fock solve with aceQ = false gets one).
!> Each run stops when |ΔE| < program.convThresh and |Δε| < program.epsThresh
!> (ε converges only like √ΔE).
!>
!> Checks (electronic energies; the nuclear repulsion 1/R is printed apart):
!> - mesh: element bounds of the prototype (printed to 4 decimals there)
!> - GroundSolver_PreconditionAction: (T + 2)·P(v; θ = −2) = v and
!>   (T + σ_min)·P(v; θ = 1) = v
!> - E_el (ACE run) vs the prototype on the same mesh, −1.8479538609 (1e-8),
!>   and vs the HF limit −1.8479152851 (1e-4: Galerkin point nuclei are not
!>   variational, this coarse mesh lies 3.9e-5 below)
!> - ε(1σg) vs the prototype: −0.59467108933 from a re-run converged to
!>   |ΔE| < 1e-13 and |Δε| < 1e-10 (1e-8; this run stops at |Δε| < 1e-9 and
!>   LOBPCG tol 1e-8); the value −0.5946715469 listed with the prototype
!>   comes from a run stopped at |ΔE| < 1e-11, whose ε was still 4.6e-7 off
!>   (checked to 1e-6)
!> - exact exchange and ACE: E_el agree to 1e-8, ε to 1e-8
!> - every LOBPCG solve converged, with few iterations after the first SCF
!>   step (warm start + preconditioner)
program T_Fem_03_H2ScfRhf
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Grid
  use M_Grid_Fem
  use M_SysKinetic
  use M_SysPotential_Fem_Coulomb, only: SysPotential_Fem_Coulomb_nuclearRepulsion
  use M_Orbs
  use M_GroundSolver
  use testdrive, only: check, error_type

  implicit none

  !> Prototype on the same mesh (Galerkin, multipole lmax 6 about the origin)
  real(R64), parameter :: E_EL_PROTOTYPE = -1.8479538609_R64
  real(R64), parameter :: EPS_PROTOTYPE = -0.59467108933_R64
  real(R64), parameter :: EPS_PROTOTYPE_LISTED = -0.5946715469_R64
  !> Hartree–Fock limit of the electronic energy at R = 1.4
  real(R64), parameter :: E_EL_HF_LIMIT = -1.8479152851_R64
  !> Element bounds of the prototype mesh (x = y, and z), rounded to 4 decimals
  real(R64), parameter :: BOUNDS_X(0:12) = [-10.0_R64, -7.3315_R64, -4.663_R64, -2.3465_R64, -1.0603_R64, &
                                            -0.3702_R64, 0.0_R64, 0.3702_R64, 1.0603_R64, 2.3465_R64, &
                                            4.663_R64, 7.3315_R64, 10.0_R64]
  real(R64), parameter :: BOUNDS_Z(0:15) = [-10.0_R64, -7.4482_R64, -4.8963_R64, -2.8288_R64, -1.6814_R64, &
                                            -1.0488_R64, -0.7_R64, -0.3113_R64, 0.3113_R64, 0.7_R64, &
                                            1.0488_R64, 1.6814_R64, 2.8288_R64, 4.8963_R64, 7.4482_R64, 10.0_R64]
  !> Largest accepted LOBPCG iteration count of an SCF step after the first
  !> (warm start from the previous eigenvectors; measured: 7)
  integer(I32), parameter :: MAX_LOBPCG_ITER_WARM = 15

  type(error_type), allocatable :: error
  complex(R64), allocatable :: vec(:), precVec(:), kinVec(:)
  real(R64) :: convThresh, epsThresh, energyAce, energyExact, epsAce(1), epsExact(1), precRes(2), boundsDev
  integer(I32) :: maxScfIterations, nScfAce, nScfExact, maxWarmAce, maxWarmExact
  logical :: convergedAce, convergedExact, foundQ

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call SetupFemScf("test/simulationTests/Fem/T_Fem_03_H2ScfRhf.json")

  call Say_Fabricate("program")
  convThresh = Json_Get("program.convThresh", 1e-12_R64)
  epsThresh = Json_Get("program.epsThresh", 1e-9_R64)
  maxScfIterations = Json_Get("program.maxScfIterations", 40)

  !==========================================
  call Say_Section("perform action")
  !==========================================

  write (*, '(1X, A, 3(I0, 1X), A, I0)') "mesh: n =", Grid_Fem_nx, Grid_Fem_ny, Grid_Fem_nz, ", N = ", Grid_nPoints
  boundsDev = max(maxval(abs(Grid_Fem_axes(1) % bounds - BOUNDS_X)), maxval(abs(Grid_Fem_axes(2) % bounds - BOUNDS_X)), &
                  maxval(abs(Grid_Fem_axes(3) % bounds - BOUNDS_Z)))
  write (*, *) "max deviation of the element bounds from the prototype (4 decimals): ", boundsDev

  ! Preconditioner contract: (T + max(-shift, sigmaMin)) P v = v
  allocate (vec(Grid_nPoints), precVec(Grid_nPoints), kinVec(Grid_nPoints))
  vec(:) = Orbs_orbs(:, 1)
  call GroundSolver_PreconditionAction(precVec, vec, -2.0_R64, 0.0_R64)
  call SysKinetic_MultiplyWithKineticOp(kinVec, precVec, 0.0_R64)
  precRes(1) = MetricNorm(kinVec + 2.0_R64 * precVec - vec) / MetricNorm(vec)
  call GroundSolver_PreconditionAction(precVec, vec, 1.0_R64, 0.0_R64)
  call SysKinetic_MultiplyWithKineticOp(kinVec, precVec, 0.0_R64)
  precRes(2) = MetricNorm(kinVec + GroundSolver_precShiftMin * precVec - vec) / MetricNorm(vec)
  write (*, *) "preconditioner residuals (theta = -2, theta = 1): ", precRes

  ! Run 1: adaptively compressed exchange (aceQ = true in the JSON) from the
  ! start orbitals (no energy yet)
  energyAce = huge(1.0_R64)
  epsAce = huge(1.0_R64)
  call RunScf(nScfAce, maxWarmAce, convergedAce, energyAce, epsAce, "ACE", maxScfIterations, convThresh, epsThresh, &
              firstStepWarmQ=.false.)
  call PrintSummary("ACE", energyAce, epsAce(1))

  ! Converged orbital for the dynamics tests of this folder (orb01_01.in)
  call Orbs_SaveOrbs(Orbs_orbs)

  ! Run 2: exact exchange, continued from the converged ACE orbitals and
  ! LOBPCG block (the ground solver re-reads its keys; its work arrays from
  ! Setup stay valid); ΔE and Δε of the first step are taken against ACE
  call jsonFile % update("groundSolver.scf.stdImpl.aceQ", .false., foundQ)
  if (.not. foundQ) error stop "T_Fem_03_H2ScfRhf: groundSolver.scf.stdImpl.aceQ not found"
  call GroundSolver_Fabricate
  energyExact = energyAce
  epsExact = epsAce
  call RunScf(nScfExact, maxWarmExact, convergedExact, energyExact, epsExact, "exact exchange", maxScfIterations, &
              convThresh, epsThresh, firstStepWarmQ=.true.)
  call PrintSummary("exact exchange", energyExact, epsExact(1))

  !==========================================
  call Say_Section("printout results")
  !==========================================

  print *
  print *, "E_el (ACE):             ", energyAce, " SCF iterations ", nScfAce, ", max warm LOBPCG iterations ", &
    maxWarmAce
  print *, "E_el (exact exchange):  ", energyExact, " SCF iterations ", nScfExact, ", max warm LOBPCG iterations ", &
    maxWarmExact
  print *, "E_el prototype:         ", E_EL_PROTOTYPE, " difference ", energyAce - E_EL_PROTOTYPE
  print *, "E_el HF limit:          ", E_EL_HF_LIMIT, " difference ", energyAce - E_EL_HF_LIMIT
  print *, "eps (ACE, exact):       ", epsAce(1), epsExact(1)
  print *, "eps prototype:          ", EPS_PROTOTYPE, " difference ", epsAce(1) - EPS_PROTOTYPE
  print *, "eps prototype (listed): ", EPS_PROTOTYPE_LISTED, " difference ", epsAce(1) - EPS_PROTOTYPE_LISTED

  call check(error, Grid_nPoints, 130331)
  if (allocated(error)) error stop "T_Fem_03_H2ScfRhf failure: number of nodes differs from the prototype"

  call check(error, boundsDev < 6e-5_R64)
  if (allocated(error)) error stop "T_Fem_03_H2ScfRhf failure: element bounds differ from the prototype"

  call check(error, maxval(precRes) < 1e-10_R64)
  if (allocated(error)) error stop "T_Fem_03_H2ScfRhf failure: GroundSolver_PreconditionAction"

  call check(error, convergedAce .and. convergedExact)
  if (allocated(error)) error stop "T_Fem_03_H2ScfRhf failure: a LOBPCG solve did not converge"

  call check(error, maxWarmAce <= MAX_LOBPCG_ITER_WARM .and. maxWarmExact <= MAX_LOBPCG_ITER_WARM)
  if (allocated(error)) error stop "T_Fem_03_H2ScfRhf failure: too many LOBPCG iterations in an SCF step"

  call check(error, energyAce, E_EL_PROTOTYPE, thr=1e-8_R64)
  if (allocated(error)) error stop "T_Fem_03_H2ScfRhf failure: E_el vs prototype"

  call check(error, energyAce, E_EL_HF_LIMIT, thr=1e-4_R64)
  if (allocated(error)) error stop "T_Fem_03_H2ScfRhf failure: E_el vs HF limit"

  call check(error, epsAce(1), EPS_PROTOTYPE, thr=1e-8_R64)
  if (allocated(error)) error stop "T_Fem_03_H2ScfRhf failure: orbital energy vs prototype"

  call check(error, epsAce(1), EPS_PROTOTYPE_LISTED, thr=1e-6_R64)
  if (allocated(error)) error stop "T_Fem_03_H2ScfRhf failure: orbital energy vs listed prototype value"

  call check(error, energyExact, energyAce, thr=1e-8_R64)
  if (allocated(error)) error stop "T_Fem_03_H2ScfRhf failure: E_el exact exchange vs ACE"

  call check(error, epsExact(1), epsAce(1), thr=1e-8_R64)
  if (allocated(error)) error stop "T_Fem_03_H2ScfRhf failure: orbital energy exact exchange vs ACE"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  !> Electronic, nuclear and total energy and the orbital energy.
  subroutine PrintSummary(label, energy, eps)
    character(len=*), intent(in) :: label
    real(R64), intent(in) :: energy
    real(R64), intent(in) :: eps

    write (*, *)
    write (*, '(1X, 3A)') "summary (", label, ")"
    write (*, '(1X, A, F18.12)') "  electronic energy  ", energy
    write (*, '(1X, A, F18.12)') "  nuclear repulsion  ", SysPotential_Fem_Coulomb_nuclearRepulsion
    write (*, '(1X, A, F18.12)') "  total energy       ", energy + SysPotential_Fem_Coulomb_nuclearRepulsion
    write (*, '(1X, A, F18.12)') "  eps(1sigma_g)      ", eps

  end subroutine

  !> Metric norm sqrt(<f|f>).
  function MetricNorm(f) result(res)
    complex(R64), intent(in), contiguous :: f(:)
    real(R64) :: res

    res = sqrt(real(Grid_InnerProduct(f, f), kind=R64))

  end function

  include 'FemScf.inc'

end program
