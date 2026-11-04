!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> H2+ molecular ion (R = 2) on the finite-element grid (grid.fem): lowest
!> eigenvalue of H = T + V_nuc with the preconditioned block eigensolver
!> diagonalizerList.lobpcg.
!>
!> Mesh (identical to the numpy prototype of the FEM design): box ±12 in x
!> and y, ±13 in z, nLocals 6, centers (0, 0, ∓1) with hMin 0.2 each,
!> growth 1.6, hMax 3.0 → n = (69, 69, 94), N = 447534. The two centers share
!> their x and y coordinates (merged into one graded coordinate per axis);
!> both protons sit at element vertices.
!>
!> - H is applied by a program-local matvec (SysKinetic_MultiplyWithKineticOp
!>   plus the static Galerkin point-nucleus potential of
!>   sysPotential.fem.coulomb, Z = 1 on both centers)
!> - preconditioner (program-local, as GroundSolver_PreconditionAction):
!>   (T + max(−θ, σ_min))⁻¹ through the exact FDM solve
!>   SysKinetic_SolveShifted, σ_min = program.precShiftMin
!> - metric: the diagonal FEM mass, metricWeights = Grid_Fem_weights
!> - start vector: the normalized gerade LCAO 1s(A) + 1s(B) (method.sb +
!>   orbsInit.fem.sto)
!>
!> Checks (electronic energy; the nuclear repulsion 1/R = 0.5 is excluded):
!> - E against the prototype on the same mesh, −1.102636133105 (1e-9): both
!>   implement the same discretization
!> - E against the exact 1sσg energy −1.1026342144949 (3e-6; the Galerkin
!>   point nuclei put this mesh 1.9e-6 below it, they are not variational)
!> - LOBPCG converged within a few tens of iterations
!> - ⟨ψ|H|ψ⟩ of method.sb (its own H) at the eigenvector equals E, unit
!>   metric norm, and the gerade symmetry ψ(x, y, −z) = ψ(x, y, z) of the
!>   eigenvector
program T_Fem_02_H2plusLobpcg
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

  real(R64), parameter :: E_PROTOTYPE = -1.102636133105_R64
  !> Exact electronic 1sσg energy of H2+ at R = 2 (Bates et al.)
  real(R64), parameter :: E_EXACT = -1.1026342144949_R64
  !> Largest accepted LOBPCG iteration count (measured: 10; without the
  !> preconditioner this stiff operator needs hundreds)
  integer(I32), parameter :: MAX_LOBPCG_ITER = 25

  type(error_type), allocatable :: error
  type(T_DiagonalizerList_FabricateInput) :: DiagonalizerListInput(1)
  complex(R64), allocatable :: externalPotential(:), vOrb(:)
  real(R64) :: precShiftMin, energy, energyMethod, normEvec, parityDeviation, meshAsymmetry
  integer(I32) :: nIter, nMatVec, ix, iy, iz
  logical :: convergedQ

  !==========================================
  call Say_Hello
  call Say_Section("start")
  !==========================================

  call Json_LoadJsonFile(manualFileName_="test/simulationTests/Fem/T_Fem_02_H2plusLobpcg.json", &
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

  ! Start block: the normalized gerade LCAO of method.sb (the guard columns
  ! are filled by LOBPCG)
  call DiagonalizerList(1) % e % SetStartVectors(reshape(Method_state, [Grid_nPoints, 1]))

  call DiagonalizerList(1) % e % Diagonalize(0.0_R64, .true.)

  select type (e => DiagonalizerList(1) % e)
  type is (T_DiagonalizerList_E_Lobpcg)
    nIter = e % nIter
    nMatVec = e % nMatVec
    convergedQ = e % convergedQ
  class default
    error stop "T_Fem_02_H2plusLobpcg: diagonalizerList must be lobpcg"
  end select

  energy = DiagonalizerList(1) % e % evals(1)

  ! The eigenvector as method.sb state: its own H gives <psi|H|psi>
  Method_state(:) = DiagonalizerList(1) % e % evecs(:, 1)
  normEvec = sqrt(real(Grid_InnerProduct(Method_state, Method_state), kind=R64))
  energyMethod = Method_GetEnergy(0.0_R64)

  ! Gerade: the z mesh is symmetric about 0 (node iz <-> nz + 1 - iz, up to
  ! the round-off of the grading, meshAsymmetry)
  meshAsymmetry = maxval(abs(Grid_Fem_axes(3) % points(:) + Grid_Fem_axes(3) % points(Grid_Fem_nz:1:-1)))
  parityDeviation = 0.0_R64
  do iz = 1, Grid_Fem_nz
    do iy = 1, Grid_Fem_ny
      do ix = 1, Grid_Fem_nx
        parityDeviation = max(parityDeviation, abs(Method_state(Grid_Fem_Index(ix, iy, iz)) &
                                                   - Method_state(Grid_Fem_Index(ix, iy, Grid_Fem_nz + 1 - iz))))
      end do
    end do
  end do
  parityDeviation = parityDeviation / maxval(abs(Method_state))

  !==========================================
  call Say_Section("printout results")
  !==========================================

  print *
  print *, "lowest eigenvalue:            ", energy
  print *, "prototype (same mesh):        ", E_PROTOTYPE, " difference ", energy - E_PROTOTYPE
  print *, "exact:                        ", E_EXACT, " difference ", energy - E_EXACT
  print *, "<psi|H|psi> (method.sb):      ", energyMethod
  print *, "metric norm of eigenvector:   ", normEvec
  print *, "max |z(iz) + z(nz+1-iz)|:     ", meshAsymmetry
  print *, "max |psi(z) - psi(-z)|/max|psi|: ", parityDeviation
  print *, "LOBPCG iterations, matvecs, converged: ", nIter, nMatVec, convergedQ

  call check(error, Grid_nPoints, 447534)
  if (allocated(error)) error stop "T_Fem_02_H2plusLobpcg failure: mesh differs from the prototype"

  call check(error, convergedQ)
  if (allocated(error)) error stop "T_Fem_02_H2plusLobpcg failure: LOBPCG not converged"

  call check(error, nIter <= MAX_LOBPCG_ITER)
  if (allocated(error)) error stop "T_Fem_02_H2plusLobpcg failure: too many LOBPCG iterations"

  call check(error, energy, E_PROTOTYPE, thr=1e-9_R64)
  if (allocated(error)) error stop "T_Fem_02_H2plusLobpcg failure: energy vs prototype"

  call check(error, energy, E_EXACT, thr=3e-6_R64)
  if (allocated(error)) error stop "T_Fem_02_H2plusLobpcg failure: energy vs exact"

  call check(error, energyMethod, energy, thr=1e-10_R64)
  if (allocated(error)) error stop "T_Fem_02_H2plusLobpcg failure: method.sb energy"

  call check(error, normEvec, 1.0_R64, thr=1e-12_R64)
  if (allocated(error)) error stop "T_Fem_02_H2plusLobpcg failure: eigenvector norm"

  call check(error, meshAsymmetry < 1e-10_R64 .and. parityDeviation < 1e-6_R64)
  if (allocated(error)) error stop "T_Fem_02_H2plusLobpcg failure: eigenvector not gerade"

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
