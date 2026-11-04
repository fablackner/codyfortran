!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> H2O restricted Hartree–Fock ground state (O at the origin, H at
!> (±1.4305227, 0, −1.1071202), 10 electrons in 5 doubly occupied orbitals)
!> on a small finite-element mesh (grid.fem) with groundSolver.scf.stdImpl
!> (method tdhx, restricted orbitals [5, 5], adaptively compressed exchange,
!> Hartree-potential mixing with α = 0.7) and the preconditioned eigensolver
!> diagonalizerList.lobpcg: ApplyPrecOnVec => GroundSolver_PreconditionAction
!> ((T + max(−θ, σ_min))⁻¹ by the exact FDM solve SysKinetic_SolveShifted),
!> metricWeights = Grid_Fem_weights, no guard columns and tol 1e-7 (nGuard 0
!> instead of 2: ε equal to 2e-10 with a third fewer Fock products; tol 1e-7
!> instead of 1e-8: ε equal to 2e-8 with 20 % fewer products, since the
!> occupied-virtual gap, not the HOMO-1/HOMO gap, sets the accuracy of the
!> occupied space). Wiring and SCF loop (|ΔE| < program.convThresh and
!> max|Δε| < program.epsThresh): FemScf.inc, shared with T_Fem_03.
!> Covers the multi-center mesh (three nuclei, merged center projections on y
!> and z), dipolar multipole boundary values, the STO start guess
!> (orbsInit.fem.sto: O 1s, 2s, 2p plus H 1s combinations) and the
!> five-orbital ACE.
!> SCF path: the plain SCF of the prototype (α = 1) converges as well (23
!> steps, 804 Fock products at tol 1e-8); its leading error mode oscillates
!> with a factor ≈ −0.46, which α ≈ 1/1.46 removes (22 steps, 736 products;
!> slower modes remain). DIIS on the Hartree potential is erratic here (51
!> steps).
!>
!> Mesh (identical to the numpy prototype of the FEM design): box ±8, nLocals
!> 4, hMin 0.1 at O and 0.3 at the H atoms, growth 1.8, hMax 3.0, multipole
!> boundary values about the origin with lmax 6 → n = (47, 35, 44), N = 72380.
!>
!> This is a regression and integration test on a COARSE mesh, never a test of
!> the Hartree–Fock limit (−76.0675): the O core is under-resolved and
!> Galerkin point nuclei are not variational, so this mesh lies 0.34 Ha BELOW
!> the HF limit (the converged prototype mesh with N = 1.1e6 gives −76.0692).
!>
!> Checks against the prototype on the same mesh (independent numpy
!> implementation of the same discretization):
!> - E_tot = E_el + E_nn: −76.4037559605 as listed with the prototype (1e-6)
!>   and −76.403755960801 from a re-run converged to |ΔE| < 1e-13 and
!>   |Δε| < 1e-10 (1e-9)
!> - E_nn = SysPotential_Fem_Coulomb_nuclearRepulsion vs 9.194689103233
!>   (Z_I Z_J/R_IJ of the coordinates above, 1e-10). The prototype listing
!>   gives E_nn = 9.1949450734 and E_el = −85.5987010339: its print script
!>   used a hard-coded R_OH = 1.8088452 instead of the 1.8088975 of these
!>   coordinates; its E_tot is right (computed from the coordinates)
!> - orbital energies vs the re-run (5e-7; this run stops at |Δε| < 1e-7)
!>   and vs the listed values, which are not SCF-converged beyond ~7e-6 (1e-5)
!> - every LOBPCG solve converged
program T_Fem_04_H2oScfRhf
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Grid
  use M_Grid_Fem
  use M_SysPotential_Fem_Coulomb, only: SysPotential_Fem_Coulomb_nuclearRepulsion
  use testdrive, only: check, error_type

  implicit none

  integer(I32), parameter :: N_OCC = 5
  !> Prototype on the same mesh: listed values and the tightly converged re-run
  real(R64), parameter :: E_TOT_PROTOTYPE_LISTED = -76.4037559605_R64
  real(R64), parameter :: E_TOT_PROTOTYPE = -76.403755960801_R64
  real(R64), parameter :: E_NN = 9.194689103233_R64
  real(R64), parameter :: EPS_PROTOTYPE(N_OCC) = [-20.7051853746_R64, -1.35922237664_R64, -0.71676850941_R64, &
                                                  -0.58473765118_R64, -0.50947403454_R64]
  real(R64), parameter :: EPS_PROTOTYPE_LISTED(N_OCC) = [-20.70519246_R64, -1.35922481_R64, -0.71676879_R64, &
                                                         -0.58473888_R64, -0.50947870_R64]
  !> Largest accepted LOBPCG iteration count of an SCF step after the first
  !> (warm start from the previous eigenvectors; measured: 10)
  integer(I32), parameter :: MAX_LOBPCG_ITER_WARM = 25

  type(error_type), allocatable :: error
  real(R64) :: convThresh, epsThresh, energy, eps(N_OCC)
  integer(I32) :: maxScfIterations, nScf, maxWarmIter, i
  logical :: allConvergedQ

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call SetupFemScf("test/simulationTests/Fem/T_Fem_04_H2oScfRhf.json")

  call Say_Fabricate("program")
  convThresh = Json_Get("program.convThresh", 1e-10_R64)
  epsThresh = Json_Get("program.epsThresh", 1e-7_R64)
  maxScfIterations = Json_Get("program.maxScfIterations", 60)

  !==========================================
  call Say_Section("perform action")
  !==========================================

  write (*, '(1X, A, 3(I0, 1X), A, I0)') "mesh: n =", Grid_Fem_nx, Grid_Fem_ny, Grid_Fem_nz, ", N = ", Grid_nPoints

  energy = huge(1.0_R64)
  eps(:) = huge(1.0_R64)
  call RunScf(nScf, maxWarmIter, allConvergedQ, energy, eps, "ACE", maxScfIterations, convThresh, epsThresh, &
              firstStepWarmQ=.false.)

  !==========================================
  call Say_Section("printout results")
  !==========================================

  print *
  write (*, '(1X, A)') "summary"
  write (*, '(1X, A, F18.10)') "  electronic energy  ", energy
  write (*, '(1X, A, F18.10)') "  nuclear repulsion  ", SysPotential_Fem_Coulomb_nuclearRepulsion
  write (*, '(1X, A, F18.10)') "  total energy       ", energy + SysPotential_Fem_Coulomb_nuclearRepulsion
  do i = 1, N_OCC
    write (*, '(1X, A, I0, A, F18.10, A, ES10.2, A, ES10.2)') "  eps(", i, ")          ", eps(i), &
      "   vs prototype ", eps(i) - EPS_PROTOTYPE(i), ", vs listed ", eps(i) - EPS_PROTOTYPE_LISTED(i)
  end do
  print *
  print *, "E_tot - prototype (listed, re-run): ", energy + SysPotential_Fem_Coulomb_nuclearRepulsion - E_TOT_PROTOTYPE_LISTED, &
    energy + SysPotential_Fem_Coulomb_nuclearRepulsion - E_TOT_PROTOTYPE
  print *, "E_nn - Z_I Z_J/R_IJ:                ", SysPotential_Fem_Coulomb_nuclearRepulsion - E_NN
  print *, "max LOBPCG iterations of the warm-started SCF steps: ", maxWarmIter

  call check(error, Grid_Fem_nx .eq. 47 .and. Grid_Fem_ny .eq. 35 .and. Grid_Fem_nz .eq. 44)
  if (allocated(error)) error stop "T_Fem_04_H2oScfRhf failure: mesh differs from the prototype"

  call check(error, allConvergedQ)
  if (allocated(error)) error stop "T_Fem_04_H2oScfRhf failure: a LOBPCG solve did not converge"

  call check(error, maxWarmIter <= MAX_LOBPCG_ITER_WARM)
  if (allocated(error)) error stop "T_Fem_04_H2oScfRhf failure: too many LOBPCG iterations in an SCF step"

  call check(error, SysPotential_Fem_Coulomb_nuclearRepulsion, E_NN, thr=1e-10_R64)
  if (allocated(error)) error stop "T_Fem_04_H2oScfRhf failure: nuclear repulsion"

  call check(error, energy + SysPotential_Fem_Coulomb_nuclearRepulsion, E_TOT_PROTOTYPE_LISTED, thr=1e-6_R64)
  if (allocated(error)) error stop "T_Fem_04_H2oScfRhf failure: E_tot vs listed prototype value"

  call check(error, energy + SysPotential_Fem_Coulomb_nuclearRepulsion, E_TOT_PROTOTYPE, thr=1e-9_R64)
  if (allocated(error)) error stop "T_Fem_04_H2oScfRhf failure: E_tot vs prototype"

  call check(error, maxval(abs(eps - EPS_PROTOTYPE)) < 5e-7_R64)
  if (allocated(error)) error stop "T_Fem_04_H2oScfRhf failure: orbital energies vs prototype"

  call check(error, maxval(abs(eps - EPS_PROTOTYPE_LISTED)) < 1e-5_R64)
  if (allocated(error)) error stop "T_Fem_04_H2oScfRhf failure: orbital energies vs listed prototype values"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  include 'FemScf.inc'

end program
