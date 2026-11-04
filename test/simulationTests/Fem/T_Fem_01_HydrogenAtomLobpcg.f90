!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Hydrogen atom on the finite-element grid (grid.fem): lowest eigenvalue of
!> H = T + V_nuc with the preconditioned block eigensolver
!> diagonalizerList.lobpcg.
!>
!> Mesh (identical to the numpy prototype of the FEM design): box ±12 on every
!> axis, nLocals 6, one center at the origin with hMin 0.4, growth 1.6,
!> hMax 3.0 → n = 59 nodes per axis, N = 205379.
!>
!> - H is applied by a program-local matvec (SysKinetic_MultiplyWithKineticOp
!>   plus the static Galerkin point-nucleus potential of
!>   sysPotential.fem.coulomb, Z = 1)
!> - preconditioner (program-local, as GroundSolver_PreconditionAction):
!>   (T + max(−θ, σ_min))⁻¹ through the exact FDM solve
!>   SysKinetic_SolveShifted, σ_min = program.precShiftMin
!> - metric: the diagonal FEM mass, metricWeights = Grid_Fem_weights
!> - start vector: a normalized 1s Slater function with ζ = 1.3
!>   (method.sb + orbsInit.fem.sto)
!>
!> Checks:
!> - E against the prototype on the same mesh, −0.500012546113 (1e-9): both
!>   implement the same discretization, so they must agree to the solver
!>   tolerance
!> - |E + 0.5| < 2e-5 (exact −0.5; Galerkin point nuclei are not variational,
!>   this coarse mesh lies 1.25e-5 below)
!> - LOBPCG converged within a few tens of iterations (an unpreconditioned
!>   solve of this stiff operator needs hundreds)
!> - ⟨ψ|H|ψ⟩ of method.sb (its own H) at the eigenvector equals E, and the
!>   eigenvector has unit metric norm
program T_Fem_01_HydrogenAtomLobpcg
  use M_Utils_Types
  use M_Utils_Say
  use M_Utils_Json
  use M_Grid
  use M_Grid_Fem
  use M_SysKinetic
  use M_SysPotential
  use M_OrbsInit
  use M_Method
  use M_DiagonalizerList
  use M_DiagonalizerList_Lobpcg, only: T_DiagonalizerList_E_Lobpcg
  use testdrive, only: check, error_type

  implicit none

  !> Exact −0.5 lies 1.25e-5 above the Galerkin value of this mesh
  real(R64), parameter :: E_PROTOTYPE = -0.500012546113_R64
  real(R64), parameter :: E_EXACT = -0.5_R64
  !> Largest accepted LOBPCG iteration count (measured: 7; without the
  !> preconditioner this stiff operator needs hundreds)
  integer(I32), parameter :: MAX_LOBPCG_ITER = 20

  type(error_type), allocatable :: error
  type(T_DiagonalizerList_FabricateInput) :: DiagonalizerListInput(1)
  complex(R64), allocatable :: externalPotential(:), vOrb(:)
  real(R64) :: precShiftMin, energy, energyMethod, normEvec
  integer(I32) :: nIter, nMatVec
  logical :: convergedQ

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call Json_LoadJsonFile(manualFileName_="test/simulationTests/Fem/T_Fem_01_HydrogenAtomLobpcg.json", &
                         relativeToProjectDirQ_=.true.)

  call Grid_Fabricate
  call SysKinetic_Fabricate
  call SysPotential_Fabricate
  call OrbsInit_Fabricate
  call Method_Fabricate
  DiagonalizerListInput(1) % ApplyMatOnVec => ApplyHamiltonian
  DiagonalizerListInput(1) % ApplyPrecOnVec => ApplyPreconditioner
  DiagonalizerListInput(1) % dim = Grid_nPoints
  call DiagonalizerList_Fabricate(DiagonalizerListInput)

  call Say_Fabricate("program")
  precShiftMin = Json_Get("program.precShiftMin", 0.5_R64)

  call Grid_Setup
  call SysKinetic_Setup
  call SysPotential_Setup
  call OrbsInit_Setup
  call Method_Setup
  call DiagonalizerList_Setup

  ! LOBPCG works in the FEM metric (diagonal mass, BLAS-3 Gram matrices)
  DiagonalizerList(1) % e % metricWeights = Grid_Fem_weights

  ! Static nuclear potential for the program-local matvec
  call SysPotential_FillExternalPotential(externalPotential, 0.0_R64)
  allocate (vOrb(Grid_nPoints))

  !==========================================
  call Say_Section("perform action")
  !==========================================

  write (*, '(1X, A, 3(I0, 1X), A, I0)') "mesh: n =", Grid_Fem_nx, Grid_Fem_ny, Grid_Fem_nz, ", N = ", Grid_nPoints

  ! Start block: the normalized 1s Slater function of method.sb (the guard
  ! columns are filled by LOBPCG)
  call DiagonalizerList(1) % e % SetStartVectors(reshape(Method_state, [Grid_nPoints, 1]))

  call DiagonalizerList(1) % e % Diagonalize(0.0_R64, .true.)

  select type (e => DiagonalizerList(1) % e)
  type is (T_DiagonalizerList_E_Lobpcg)
    nIter = e % nIter
    nMatVec = e % nMatVec
    convergedQ = e % convergedQ
  class default
    error stop "T_Fem_01_HydrogenAtomLobpcg: diagonalizerList must be lobpcg"
  end select

  energy = DiagonalizerList(1) % e % evals(1)

  ! The eigenvector as method.sb state: its own H gives <psi|H|psi>
  Method_state(:) = DiagonalizerList(1) % e % evecs(:, 1)
  normEvec = sqrt(real(Grid_InnerProduct(Method_state, Method_state), kind=R64))
  energyMethod = Method_GetEnergy(0.0_R64)

  !==========================================
  call Say_Section("printout results")
  !==========================================

  print *
  print *, "lowest eigenvalue:            ", energy
  print *, "prototype (same mesh):        ", E_PROTOTYPE, " difference ", energy - E_PROTOTYPE
  print *, "exact:                        ", E_EXACT, " difference ", energy - E_EXACT
  print *, "<psi|H|psi> (method.sb):      ", energyMethod
  print *, "metric norm of eigenvector:   ", normEvec
  print *, "LOBPCG iterations, matvecs, converged: ", nIter, nMatVec, convergedQ

  call check(error, Grid_nPoints, 205379)
  if (allocated(error)) error stop "T_Fem_01_HydrogenAtomLobpcg failure: mesh differs from the prototype"

  call check(error, convergedQ)
  if (allocated(error)) error stop "T_Fem_01_HydrogenAtomLobpcg failure: LOBPCG not converged"

  call check(error, nIter <= MAX_LOBPCG_ITER)
  if (allocated(error)) error stop "T_Fem_01_HydrogenAtomLobpcg failure: too many LOBPCG iterations"

  call check(error, energy, E_PROTOTYPE, thr=1e-9_R64)
  if (allocated(error)) error stop "T_Fem_01_HydrogenAtomLobpcg failure: energy vs prototype"

  call check(error, abs(energy - E_EXACT) < 2e-5_R64)
  if (allocated(error)) error stop "T_Fem_01_HydrogenAtomLobpcg failure: energy vs exact"

  call check(error, energyMethod, energy, thr=1e-10_R64)
  if (allocated(error)) error stop "T_Fem_01_HydrogenAtomLobpcg failure: method.sb energy"

  call check(error, normEvec, 1.0_R64, thr=1e-12_R64)
  if (allocated(error)) error stop "T_Fem_01_HydrogenAtomLobpcg failure: eigenvector norm"

  !==========================================
  call Say_Goodbye
  !==========================================

contains

  !> H·state = T·state + V_nuc·state (static potential filled after setup).
  subroutine ApplyHamiltonian(dState, state, time)
    complex(R64), intent(out), contiguous, target :: dState(:)
    complex(R64), intent(in), contiguous, target :: state(:)
    real(R64), intent(in) :: time

    call SysKinetic_MultiplyWithKineticOp(dState, state, time)
    call SysPotential_MultiplyWithExternalPotential(vOrb, externalPotential, state)
    dState(:) = dState(:) + vOrb(:)

  end subroutine

  !> (T + max(−shift, σ_min))⁻¹·state, the kinetic preconditioner of
  !> GroundSolver_PreconditionAction (no ground solver is fabricated here).
  subroutine ApplyPreconditioner(dState, state, shift, time)
    use M_Utils_UnusedVariables

    complex(R64), intent(out), contiguous, target :: dState(:)
    complex(R64), intent(in), contiguous, target :: state(:)
    real(R64), intent(in) :: shift
    real(R64), intent(in) :: time

    if (.false.) call UnusedVariables_Mark(time)

    call SysKinetic_SolveShifted(dState, state, cmplx(max(-shift, precShiftMin), 0.0_R64, kind=R64))

  end subroutine

end program
