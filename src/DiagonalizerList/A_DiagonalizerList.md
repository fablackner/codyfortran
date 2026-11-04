# DiagonalizerList Module — Agent Onboarding Guide

## Purpose

The **DiagonalizerList** module provides a unified interface for eigenvalue problem solvers used in quantum many-body simulations. It abstracts away the differences between dense (LAPACK), iterative (ARPACK) and preconditioned block (LOBPCG) eigensolvers, allowing the rest of the codebase to request eigenpairs through a **matrix-free** interface.

**Core problem solved:** Given a linear operator Ĥ (typically a Hamiltonian), find eigenvalues λᵢ and eigenvectors |ψᵢ⟩ satisfying:

```
Ĥ |ψᵢ⟩ = λᵢ |ψᵢ⟩
```

The operator is accessed only through a callback `ApplyMatOnVec(y, x, t)` that computes y = Ĥ(t)·x. Two optional inputs serve the preconditioned solver (LOBPCG) and are ignored by LAPACK and ARPACK: a preconditioner callback `ApplyPrecOnVec(y, x, shift, t)` ≈ (Ĥ − shift)⁻¹·x and the diagonal metric `metricWeights` in which Ĥ is self-adjoint (e.g. the quadrature weights of a grid; Euclidean when unallocated, the convention of IntegratorList).

---

## Architecture Overview

```
src/DiagonalizerList/
├── M_DiagonalizerList.f90        # Public interface (abstract type + procedure pointers)
├── S_DiagonalizerList.f90        # Factory logic (JSON dispatch)
├── A_DiagonalizerList.md         # This file
├── Lapack/
│   ├── M_DiagonalizerList_Lapack.f90   # LAPACK backend interface
│   └── S_DiagonalizerList_Lapack.f90   # LAPACK implementation
├── Arpack/
│   ├── M_DiagonalizerList_Arpack.f90   # ARPACK backend interface
│   └── S_DiagonalizerList_Arpack.f90   # ARPACK implementation
└── Lobpcg/
    ├── M_DiagonalizerList_Lobpcg.f90   # LOBPCG backend interface
    └── S_DiagonalizerList_Lobpcg.f90   # LOBPCG implementation
```

### Design Pattern: Strategy + Factory

| Component | Role |
|-----------|------|
| **M_DiagonalizerList** | Abstract interface `T_DiagonalizerList_E` with deferred methods |
| **S_DiagonalizerList** | Factory that reads JSON and instantiates backends by name |
| **T_DiagonalizerList_E_Lapack** | Dense solver: assembles explicit matrix, calls LAPACK |
| **T_DiagonalizerList_E_Arpack** | Iterative solver: Arnoldi/Lanczos via ARPACK |
| **T_DiagonalizerList_E_Lobpcg** | Preconditioned block solver (LOBPCG) in a metric |

---

## Data Flow

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                           Client Code                                        │
│  1. Define callback: ApplyMatOnVec(dState, state, time)                     │
│     (optional, lobpcg: ApplyPrecOnVec; after fabrication: metricWeights)    │
│  2. Create input: T_DiagonalizerList_FabricateInput                         │
└─────────────────────────────────────────────────────────────────────────────┘
                                   │
                                   ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                    DiagonalizerList_Fabricate(input)                         │
│  1. Read JSON "diagonalizerList" block                                       │
│  2. For each child: allocate backend by name (lapack/arpack/lobpcg)         │
│  3. Bind callbacks, set dim, call backend%Fabricate()                       │
└─────────────────────────────────────────────────────────────────────────────┘
                                   │
                                   ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                       DiagonalizerList_Setup                                 │
│  Calls Setup() on each backend (allocate workspaces)                        │
└─────────────────────────────────────────────────────────────────────────────┘
                                   │
                                   ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│              DiagonalizerList(i)%e%Diagonalize(time, evecsQ)                │
│  (optionally seeded before: DiagonalizerList(i)%e%SetStartVectors(vecs))    │
│  Compute eigenpairs → store in %evals(:), %evecs(:,:), %nFound              │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## JSON Configuration

### Example: Multiple Solvers

```json
{
  "diagonalizerList": {
    "arpack_ground": {
      "nEvals": 5,
      "which": "SR",
      "nKry": 20,
      "checkConvergenceQ": true,
      "printLevel": 0
    },
    "lapack_full": {
      "nEvals": -1,
      "printLevel": 0
    }
  }
}
```

**Backend selection:** The JSON child name must contain `"lapack"`, `"arpack"` or `"lobpcg"` as a substring.

### LAPACK Parameters

| Parameter    | Type | Default | Description                       |
|--------------|------|---------|-----------------------------------|
| `nEvals`     | int  | 1       | Number of eigenvalues (-1 = all)  |
| `printLevel` | int  | 0       | Verbosity (0=silent, 1=progress)  |

### ARPACK Parameters

| Parameter           | Type   | Default      | Description                              |
|---------------------|--------|--------------|------------------------------------------|
| `nEvals`            | int    | 1            | Number of eigenvalues (-1 = all)         |
| `which`             | string | "SR"         | Spectrum region (see below)              |
| `bmat`              | string | "I"          | "I" = standard, "G" = generalized        |
| `nKry`              | int    | 2*nEvals+1   | Krylov subspace dimension                |
| `tol`               | real   | 0 (machine eps) | Ritz value convergence tolerance      |
| `checkConvergenceQ` | bool   | true         | Verify residuals after solve             |
| `warmStartQ`        | bool   | true         | Restart from previous eigenvectors       |
| `printLevel`        | int    | 0            | Verbosity level                          |

**Spectrum selection (`which`):**
- `"LM"` / `"SM"` : Largest / Smallest Magnitude
- `"LR"` / `"SR"` : Largest / Smallest Real part ← **use for ground states**
- `"LI"` / `"SI"` : Largest / Smallest Imaginary part

**Warm start (`warmStartQ`):** on repeated `Diagonalize` calls the Arnoldi
iteration is started from the normalized sum of the previously converged
eigenvectors (ARPACK `info=1` with a deterministic `resid`) instead of a
random vector. This makes repeated solves deterministic and reduces the
matvec count when the operator changes only slightly between calls (SCF
iterations). Note the returned eigenvectors are still only defined up to
gauge (phase / degenerate-subspace rotation); consumers that mix eigenvectors
with previous state must gauge-align first (see
`Orbs_AlignOnReference`).

### LOBPCG Parameters

| Parameter    | Type | Default | Description                                            |
|--------------|------|---------|--------------------------------------------------------|
| `nEvals`     | int  | 1       | Number of wanted (lowest) eigenpairs                   |
| `nGuard`     | int  | 2       | Extra block columns (accelerate convergence)           |
| `maxIter`    | int  | 200     | Iteration limit; one-line warning (no stop) if reached |
| `tol`        | real | 1e-8    | Converged if ‖r_j‖ ≤ tol·max(1, abs(θ_j)), metric norm |
| `warmStartQ` | bool | true    | Start from previous `evecs` or `SetStartVectors` seed  |
| `printLevel` | int  | 0       | 1 = per-iteration residuals                            |

Locally Optimal Block Preconditioned Conjugate Gradient (Knyazev 2001) for
the lowest `nEvals` eigenpairs of an operator that is self-adjoint in the
metric given by `metricWeights` (Euclidean if unallocated), in the robust
form of Hetmaniuk & Lehoucq (2006) / Duersch et al.
(2018):

- block X of m = nEvals + nGuard metric-orthonormal Ritz vectors, A·X stored
- residuals R = AX − XΘ (metric norms); converged columns are soft-locked
  (no new directions, still part of the Rayleigh–Ritz step); stop when all
  j ≤ nEvals satisfy ‖r_j‖ ≤ tol·max(1, |θ_j|)
- W = `ApplyPrecOnVec(r_j, θ_j)` for the active columns (W = R without
  preconditioner), metric-orthonormalized against [X, P] and within (block
  Gram–Schmidt twice, Cholesky-QR with an eigen-decomposition fallback that
  drops near-dependent directions)
- A·W: the only matrix–vector products (one per new column and iteration);
  A·X and A·P are updated linearly
- Rayleigh–Ritz on the orthonormal basis [X, P, W] (ZHEEV of a ≤ 3m × 3m
  Hermitian matrix whose [X, P] block is known from the previous step)
- P = metric-orthonormal basis of the directions in which the active
  columns moved, orthogonal to the new X (built in coefficient space)
- a preconditioned residual without a (numerically) nonzero component
  outside [X, P] — e.g. from a nearly exact shift-and-invert preconditioner,
  (A − θ_j)⁻¹ r_j ≈ x_j — is replaced by the raw residual r_j (orthogonal
  to X) instead of being dropped; only if that vanishes too is the
  direction dropped (count in `nRawResiduals`)

**Early exits.** Besides `maxIter` (one-line warning) the iteration stops
early, with a one-line warning, `convergedQ = .false.` and `stalledQ =
.true.`, when
- it stagnates: over the last 10 iterations (`STALL_WINDOW`) the largest
  relative residual of the wanted columns moved by at most 10 %
  (1 − `STALL_FACTOR`) of its smallest value there, summing all changes up
  and down. A residual that decreases by less than 10 % without fluctuating
  sits on a floor: the Ritz pairs of an operator that is self-adjoint in the
  metric only approximately reach a fixed point (see below). The sum of the
  changes keeps erratic but converging histories going: unpreconditioned
  stiff solves (T_Ne3d_17: 50–482 iterations per solve) have plateaus over
  which the residual drops by less than 10 % while it fluctuates by 40 % or
  more (smallest variation over any 10 of their iterations: 46 %, 4.6
  times the threshold). A floor that fluctuates by more (e.g. round-off
  when `tol` is below the attainable accuracy) is not detected: `maxIter`
  ends it;
- no new search direction is left.

Outputs: `evals(1:nEvals)` ascending and `evecs(:, 1:nEvals)` orthonormal
**in the metric** — the nEvals best approximations, also after an
unconverged call. `nFound` = number of converged wanted columns (nEvals on
success); the generic `convergedQ` = (`nFound .eq. nEvals`). After an
unconverged call the converged columns need not be the leading ones
(`resNorms` tells them apart). `evecs` keeps the whole block (m columns,
also for `evecsQ = .false.`) so the next call reuses the guard vectors. The
start block is the previous `evecs` (warm start) or the vectors passed to
`SetStartVectors` (any column count with `dim` rows, also linearly
dependent; e.g. a client seed before the first call), completed with
deterministic pseudo-random vectors.
Diagnostics of the last call (type `T_DiagonalizerList_E_Lobpcg`): `nIter`,
`nMatVec`, `stalledQ`, `nRawResiduals`, `resNorms`.

**Operators that are self-adjoint only approximately.** LOBPCG Hermitizes
the projected matrices; an operator with a small anti-self-adjoint part
leaves a residual floor of about that size, at which the stagnation exit
stops the solve. The exact exchange of the SCF on `grid.fem` is such an
operator: it is Hermitian only up to the multipole truncation of the
Poisson boundary values (|⟨φ|K̂g⟩ − conj⟨g|K̂φ⟩| ≈ 5e-4 relative for
`multipoleLmax` 6), and LOBPCG stalls at relative residuals 7e-8 … 3e-7.
With LOBPCG use `groundSolver.scf.stdImpl.aceQ = true` (scf.stdImpl prints
a hint once when a Fock solve with `aceQ = false` does not converge): the
ACE operator is Hermitian by construction and costs nUp
inner products instead of nUp Poisson solves per Fock action. Fresh H₂ SCF
of T_Fem_03 (2 threads): 189 Fock products, 5.5 s with ACE; exact exchange
1772 products, 61 s without the stagnation exit (10 of 14 solves at
maxIter) and 604 products, 21 s with it.

**Preconditioning matters:** without a preconditioner the iteration count
grows like the square root of spectral spread over gap, so stiff grid
operators (FEDVR, finite elements) need hundreds to thousands of
iterations; a preconditioner ≈ (A − θ)⁻¹ that captures the stiff part (e.g.
(T + σ)⁻¹ for a Fock operator, `GroundSolver_PreconditionAction`) makes
the count nearly independent of the grid.

**Code structure** (S_DiagonalizerList_Lobpcg.f90). `Diagonalize` holds
the iteration state and the main loop, which shows every exit; its
contained steps follow the bullets above and share that state by host
association, as SIL of IntegratorList does: `BuildStartBlock`,
`UpdateResiduals` (residuals, norms, convergence), `BuildSearchDirections`
(preconditioning, orthonormalization against [X, P], A·W) and
`RayleighRitzStep` (Rayleigh–Ritz, new X and P). The stagnation rule is the
pure function `StagnatedQ(history)`; `ResetKnownProjection` builds the
projected matrix of [X, P] for the next step; the orthonormalization and
the metric and block kernels take the weights, not the solver object.

---

## Key Interfaces

### ApplyMatOnVec Callback

```fortran
subroutine ApplyMatOnVec(dState, state, time)
  complex(R64), intent(out), contiguous :: dState(:)  ! Result: y = H·x
  complex(R64), intent(in), contiguous  :: state(:)   ! Input: x
  real(R64), intent(in)                 :: time       ! For time-dependent H
end subroutine
```

This is the **only** interface LAPACK and ARPACK use to access the operator.

### Preconditioner Callback (LOBPCG)

```fortran
subroutine ApplyPrecOnVec(dState, state, shift, time)
  complex(R64), intent(out), contiguous, target :: dState(:)  ! ≈ (H − shift)⁻¹·x
  complex(R64), intent(in), contiguous, target  :: state(:)   ! residual x
  real(R64), intent(in)                         :: shift      ! current Ritz value
  real(R64), intent(in)                         :: time
end subroutine
```

An unassociated pointer means identity preconditioning.

### Metric (`metricWeights`, LOBPCG)

The metric is diagonal and given by the element component `metricWeights`
(size `dim`), with exactly the contract of IntegratorList's
`metricWeights`: when allocated, inner products are
⟨f|g⟩ = Σ_i conj(f_i) w_i g_i; when unallocated, the Euclidean product.
The metric must make the operator self-adjoint: A = W⁻¹H for a Hermitian H
and diagonal quadrature weights W is self-adjoint in the weighted product,
not in the Euclidean one. Every Hermitian grid metric of the code is
diagonal (`Grid_InnerProduct` is Σ conj(f) g w with the grid's weights;
Lattice is Euclidean).
Assign the weights after fabrication and before the first `Diagonalize`
(grid weights usually exist only after `Grid_Setup`):

```fortran
call DiagonalizerList_Fabricate(inputs)
...
call Grid_Setup
DiagonalizerList(1)%e%metricWeights = Grid_Fem_weights   ! e.g. Grid_Ylm_weights on Ylm grids
```

LOBPCG sets up one weight vector per call (`metricWeights`, or unit weights
when unallocated, as SIL does) and computes all Gram and projected matrices
with ZGEMM on weight-scaled row blocks. A non-diagonal metric (consistent
FEM mass, Gaussian basis) would need a block operator B·V, not per-pair
inner products.

### Start Vectors

```fortran
call DiagonalizerList(1)%e%SetStartVectors(vecs)   ! vecs(dim, k), any k ≥ 1
call DiagonalizerList(1)%e%Diagonalize(0.0_R64, .true.)
```

`SetStartVectors` is an implemented binding of `T_DiagonalizerList_E`: the
columns of `vecs` are offered as start vectors for the next `Diagonalize`,
e.g. the current orbitals before the first SCF solve (the solver holds no
eigenvectors of its own yet). Backends may ignore them; the default is a
no-op: LAPACK is direct, and ARPACK warm-starts only from its own converged
eigenvectors. LOBPCG overrides it: it checks the row count against `dim`
and keeps the vectors as `evecs` with `nFound = 0`, so the next call starts
from them like a warm start — only when `warmStartQ` is true (the default),
the same switch as for its previous eigenvectors; with `warmStartQ: false`
the vectors are stored but ignored. `evecs` is an output; clients do not
write it.

### Fabrication Input

```fortran
type(T_DiagonalizerList_FabricateInput) :: inputs(1)
inputs(1)%ApplyMatOnVec => MyHamiltonianCallback
inputs(1)%dim = hilbert_space_size
! optional (lobpcg only):
inputs(1)%ApplyPrecOnVec => MyPreconditioner
call DiagonalizerList_Fabricate(inputs)
! optional (lobpcg only), after Grid_Setup:
DiagonalizerList(1)%e%metricWeights = myGridWeights
```

### Accessing Results

```fortran
call DiagonalizerList(1)%e%Diagonalize(0.0_R64, .true.)

! Results:
eigenvalues  => DiagonalizerList(1)%e%evals(1:nFound)   ! real(R64)
eigenvectors => DiagonalizerList(1)%e%evecs(:, 1:nFound) ! complex(R64)
num_found    =  DiagonalizerList(1)%e%nFound
converged    =  DiagonalizerList(1)%e%convergedQ  ! all nEvals converged (default .true.;
                                                   ! lobpcg: .false. at maxIter or stagnation)
```

---

## Typical Usage in CodyFortran

```fortran
! 1. Setup fabrication input (callback + dimension)
type(T_DiagonalizerList_FabricateInput) :: DiagonalizerListInput(1)
DiagonalizerListInput(1)%ApplyMatOnVec => Actions_ApplyHamiltonian
DiagonalizerListInput(1)%dim = Coeffs_nCoeffs

! 2. Fabricate (reads JSON, creates backends)
call DiagonalizerList_Fabricate(DiagonalizerListInput)

! 3. Setup (allocate workspaces)
call DiagonalizerList_Setup

! 4. Solve eigenvalue problem
call DiagonalizerList(1)%e%Diagonalize(0.0_R64, .true.)

! 5. Use results
if (DiagonalizerList(1)%e%nFound > 0) then
  ground_state_energy = DiagonalizerList(1)%e%evals(1)
  ground_state_vector = DiagonalizerList(1)%e%evecs(:, 1)
end if
```

---

## Backend Comparison

| Aspect | LAPACK | ARPACK | LOBPCG |
|--------|--------|--------|--------|
| **Matrix access** | Assembles explicit dim×dim matrix | Matrix-free (callback only) | Matrix-free + optional preconditioner |
| **Complexity** | O(dim³) | O(nKry × dim × T_matvec) | O(iterations × m × (T_matvec + T_prec)) |
| **Memory** | O(dim²) | O(nKry × dim) | O(6m × dim), m = nEvals + nGuard |
| **Metric** | Euclidean | Euclidean (non-Hermitian solver) | Diagonal `metricWeights` (Euclidean if unallocated), A self-adjoint in it |
| **Best for** | dim ≲ 1000, full spectrum | dim > 1000, partial spectrum | lowest eigenpairs of large stiff operators with a preconditioner; warm-started repeated solves (SCF) |
| **Convergence** | Deterministic | Iterative (stops if not converged) | Iterative (`convergedQ = .false.` and a warning at maxIter or on stagnation) |

**Rule of thumb:** Use ARPACK for production (large systems), LAPACK for
validation/testing, LOBPCG when a good preconditioner exists (e.g.
finite-element Fock operators, whose spectral radius makes unpreconditioned
Krylov methods hopeless).

---

## Adding a New Backend

1. **Create directory:** `src/DiagonalizerList/MyBackend/`

2. **Interface module** (`M_DiagonalizerList_MyBackend.f90`):
   ```fortran
   module M_DiagonalizerList_MyBackend
     use M_DiagonalizerList, only: T_DiagonalizerList_E

     type, extends(T_DiagonalizerList_E) :: T_DiagonalizerList_E_MyBackend
       ! backend-specific fields
     contains
       procedure :: Fabricate
       procedure :: Setup
       procedure :: Diagonalize
     end type

     interface
       module subroutine DiagonalizerList_MyBackend_Allocate(e, path)
         class(T_DiagonalizerList_E), allocatable, intent(out) :: e
         character(len=*), intent(in) :: path
       end subroutine
     end interface
   end module
   ```

3. **Implementation submodule** (`S_DiagonalizerList_MyBackend.f90`):
   - Implement `Allocate`, `Fabricate`, `Setup`, `Diagonalize`
   - Override the implemented binding `SetStartVectors` (default: no-op)
     only if the backend can start from given vectors (as LOBPCG does)

4. **Register in factory** (`S_DiagonalizerList.f90`):
   ```fortran
   use M_DiagonalizerList_MyBackend, only: DiagonalizerList_MyBackend_Allocate
   ...
   else if (index(childName, "mybackend") .ne. 0) then
     call DiagonalizerList_MyBackend_Allocate(DiagonalizerList(i)%e, "diagonalizerList."//childName)
   ```

5. **Update CMakeLists.txt** to include new source files

---

## Dependencies

| Depends On | Purpose |
|------------|---------|
| `M_Utils_Types` | R64, I32 type kinds |
| `M_Utils_Json` | JSON configuration parsing |
| `M_Utils_Say` | Logging/diagnostics |
| `M_Utils_LapackLib` | LAPACK wrapper (ZHEEVR/DSYEVR; LOBPCG: ZHEEV, ZPOTRF, ZGESVD) |
| `M_Utils_ArpackLib` | ARPACK wrapper (znaupd/zneupd) |
| `M_Utils_BlasLib` | BLAS norms for residual checks; ZGEMM block updates (LOBPCG) |
| `M_Utils_NoOpProcedures` | Default no-op for Setup pointer |

---

## Common Pitfalls

1. **Input array size mismatch:** `size(input)` must equal number of JSON children under `"diagonalizerList"`.

2. **Callback dimension mismatch:** `input(i)%dim` must match the vector size expected by `ApplyMatOnVec`.

3. **ARPACK convergence failure:** If `nFound < nEvals`, increase `nKry` or check operator conditioning.

   **LOBPCG slow or stalling:** the operator must be self-adjoint in the metric given by `metricWeights` (otherwise the residuals stall on a floor and the stagnation exit returns `convergedQ = .false.`; e.g. the exact exchange on grid.fem: use `aceQ`); without a preconditioner stiff operators need many iterations (raise `maxIter` or wire a preconditioner). Clients that continue with unconverged results check `convergedQ` (the ground solvers print a note and continue).

4. **Forgetting Setup:** Always call `DiagonalizerList_Setup` before `Diagonalize`.

5. **Time parameter:** For time-independent Hamiltonians, pass any value (e.g., `0.0_R64`).

6. **Name convention:** JSON child names must contain `"lapack"`, `"arpack"` or `"lobpcg"` substring.

---

## Testing

Relevant test files in `test/`:
- `T_FermiHubbard_01_ExactDiagTdciLapack.f90` — LAPACK on Fermi-Hubbard
- `T_FermiHubbard_02_ExactDiagTdciArpack.f90` — ARPACK on Fermi-Hubbard
- `T_BoseHubbard_01_ExactDiagTdciLapack.f90` — LAPACK on Bose-Hubbard
- `T_BoseHubbard_02_ExactDiagTdciArpack.f90` — ARPACK on Bose-Hubbard
- `T_DiagonalizerList_Lobpcg.f90` (componentTests/DiagonalizerList/Lobpcg) —
  LOBPCG on A = W⁻¹H (`metricWeights` = W) without/with preconditioner:
  dense reference eigenvalues, metric orthonormality, recomputed residuals,
  warm start;
  the Hermitian H without `metricWeights` (the Euclidean metric, unit
  weights; diagonal preconditioner, 13 iterations) against ZHEEV of H;
  on an operator with the prescribed spectrum 1, 2, 2, 2, 3, …: a cluster
  split by nEvals = 3 with nGuard 0 and 2 (22 and 19 iterations), a start
  block seeded with linearly dependent columns (20), the nearly exact
  shift-and-invert preconditioner (AK − θ + 1e-12)⁻¹ (raw residuals
  substituted, 32 iterations; formerly a one-iteration stop with wrong
  Ritz values), and the stagnation exit on the slightly non-self-adjoint
  AK + iεW⁻¹S (ε = 1e-6): stops after 21 of 1000 iterations at the floor
  8.4e-6, `convergedQ` false, `nFound` 0, Ritz values within 2e-6
- `T_Ne3d_17_GroundSolverScfLobpcgAce.f90` — LOBPCG as the SCF Fock
  diagonalizer (`metricWeights` = `Grid_Ylm_weights`, warm starts, ACE
  exchange); asserts that every solve converges (unpreconditioned, 50–482
  iterations, no false stagnation exit) and the SCF step count
- `test/simulationTests/Fem/` — LOBPCG on finite-element grids with the
  exact kinetic preconditioner (T + max(−θ, σ))⁻¹ and `metricWeights`:
  H atom and H₂⁺ (T_Fem_01/02, 7 and 10 iterations from an STO start,
  eigenvalues equal to an independent numpy implementation to 1e-11) and
  the H₂ and H₂O RHF SCF (T_Fem_03/04, `GroundSolver_PreconditionAction`)

Run tests:
```bash
cd build && ctest -R "ExactDiag|Lobpcg|T_Fem"
```

---

## Performance Tuning

### ARPACK
- **Increase `nKry`:** Better convergence, more memory. Rule: `nKry ≈ 2-3 × nEvals`.
- **Decrease `nKry`:** Faster iterations, may fail to converge.
- **`checkConvergenceQ=false`:** Skip residual verification (faster, less safe).

### LAPACK
- **`nEvals=-1`:** Compute all eigenvalues (O(dim³) cost).
- **`nEvals=k`:** Compute only lowest k (uses ZHEEVR index range).

### LOBPCG
- **Preconditioner:** the decisive factor; the shift passed is the current
  Ritz value, so shifted approximations of (A − θ)⁻¹ work best.
- **`nGuard`:** 1–3 extra columns usually minimize the matvec count; more
  columns reduce the iterations but cost matvecs.
- **`tol`:** relative to max(1, |θ|) in the metric norm; SCF energies are
  quadratic in the eigenvector error, so 1e-8 is typically enough.
- Cost per iteration besides the matvecs: O(m²) metric inner products
  (ZGEMM on weight-scaled row blocks) and O(m²·dim) block updates.

---

## Related Modules

- **Method** — Dynamics methods that use DiagonalizerList for ground-state initialization
- **Propagator/EigenExpansion** — Time evolution via diagonalization
- **GroundSolver** — SCF iteration using eigensolvers
  (`GroundSolver_PreconditionAction`, bound by `GroundSolver_Fabricate` for
  every backend, is a ready-made LOBPCG preconditioner for orbital
  eigenproblems: SCF and MCSCF Fock operators)
