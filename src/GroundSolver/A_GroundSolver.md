# GroundSolver Module

## Overview

The GroundSolver module implements **self-consistent field (SCF) ground-state solvers** for quantum many-body systems. It finds the ground state of a Hamiltonian by iteratively diagonalizing effective operators until convergence.

Available methods are **SCF** (Hartree–Fock mean field, one Fock diagonalization
per iteration) and **MCSCF** (multiconfiguration self-consistent field, one CI
diagonalization plus one Fock diagonalization per iteration).

## Architecture

```
M_GroundSolver.f90              ← Interface: procedure pointers + abstract interfaces
└── S_GroundSolver.f90          ← Dispatch: JSON routing to backends
    ├── Scf/
    │   ├── M_GroundSolver_Scf.f90      ← SCF interface
    │   ├── S_GroundSolver_Scf.f90      ← SCF dispatch
    │   ├── StdImpl/
    │   │   ├── M_GroundSolver_Scf_StdImpl.f90  ← Standard impl interface
    │   │   └── S_GroundSolver_Scf_StdImpl.f90  ← Standard impl (general grids)
    │   └── YlmOpt/
    │       ├── M_GroundSolver_Scf_YlmOpt.f90   ← Ylm-optimized interface
    │       └── S_GroundSolver_Scf_YlmOpt.f90   ← Ylm-optimized impl (spherical grids)
    └── Mcscf/
        ├── M_GroundSolver_Mcscf.f90            ← MCSCF interface
        ├── S_GroundSolver_Mcscf.f90            ← MCSCF dispatch
        ├── StdImpl/
        │   ├── M_GroundSolver_Mcscf_StdImpl.f90 ← Standard impl interface
        │   └── S_GroundSolver_Mcscf_StdImpl.f90 ← Standard impl (two diagonalizers)
        └── YlmOpt/
            ├── M_GroundSolver_Mcscf_YlmOpt.f90  ← Ylm-optimized interface
            └── S_GroundSolver_Mcscf_YlmOpt.f90  ← Ylm-optimized impl (spherical grids)
```

The module follows the **CodyFortran fabrication pattern**:
- **M_* files**: Declare procedure pointers and abstract interfaces (contracts)
- **S_* files**: Provide implementations and JSON-based dispatch logic

## Key Concepts

### SCF Iteration

Each `GroundSolver_Approach` call performs one SCF cycle:

1. **Build mean-field potentials** from current orbitals
   - Hartree (direct): Ĵ = ∫ρ(r')/|r−r'| dr' with ρ = 2·∑_j |φ_j|²; stdImpl
     solves for the total source once (the potential is linear in the
     source) instead of once per orbital (`SysInteraction_FillDirectPotential`
     with `factor_ = 2`)
   - External: V̂_ext (e.g., nuclear Coulomb −Z/r)
   - Optionally (stdImpl `aceQ`): the adaptively compressed exchange
     operator of the current orbitals (see "Adaptively Compressed Exchange")

2. **Diagonalize the Fock operator**
   - F̂ = T̂ + V̂_ext + Ĵ − K̂
   - Uses an iterative solver via DiagonalizerList (ARPACK, or the
     preconditioned LOBPCG; see "Preconditioned Eigensolver (lobpcg)")

3. **Gauge-align and mix eigenvectors with old orbitals**
   - Eigenvectors are only defined up to a per-vector phase and, within
     degenerate shells (e.g., Ne 2p), up to an arbitrary unitary rotation of
     the degenerate subspace; iterative solvers re-randomize this gauge on
     every call. `Orbs_AlignOnReference` removes the gauge freedom by
     solving the unitary Procrustes problem: with the occupied overlap matrix
     M_ij = ⟨evec_i|φ_old,j⟩ and SVD M = U·Σ·V†, apply W = U·V† to the
     eigenvectors. Without this step, any α < 1 injects a random O(1−α)
     perturbation per iteration and destroys convergence.
   - φ_new = (1−α)·φ_old + α·evec_aligned (linear default; the update is
     delegated to the Mixing module, so DIIS can be selected via JSON)
   - α < 1 damps oscillations for difficult convergence
   - With `mixTarget: "potential"` (stdImpl) the Hartree potential is mixed
     instead and the orbitals are replaced outright (see "Mixing Strategies")

4. **Orthonormalize** via Gram–Schmidt

The caller controls the convergence loop (checking energy change).

### MCSCF Iteration

Each `GroundSolver_Approach` call performs one MCSCF cycle on the packed
MCTDHX-style state (CI coefficients ⊕ orbitals):

1. **CI step**: Build h1/h2 in the current orbital basis and diagonalize the
   CI Hamiltonian H_{IJ} = ⟨Φ_I|Ĥ|Φ_J⟩ via DiagonalizerList(1); mix the
   phase-aligned ground eigenvector with the old coefficients and normalize.
   With a declared total-M sector (`coeffs.generic.totalM`, Ylm-optimized
   solver only) the mixed vector is first projected on it
   (`Coeffs_ProjectOnSector`): the diagonalizer works on full-length vectors,
   whose eigenvectors lie in the sector only up to round-off.

2. **Orbital step**: Fill the correlated 1-RDM ρ¹ from the updated CI vector,
   build the RDM-weighted Fock operator

   F̂[ρ¹] = T̂ + V̂_ext + Ĵ[ρ¹] − K̂[ρ¹]

   with Ĵ[ρ¹] the direct potential of the correlated density
   ρ(r) = ∑_{pq} ρ¹_{pq} φ_p*(r) φ_q(r) and
   K̂[ρ¹]·φ = ∑_{pq} ρ¹_{pq} W(φ_p, φ) φ_q the RDM-weighted exchange;
   diagonalize it via DiagonalizerList(2), gauge-align the eigenvectors with
   the old orbitals (`Orbs_AlignOnReference`), mix, copy
   spin-up → spin-down, orthonormalize.

For a single determinant (ρ¹ = identity on the occupied space) the orbital
step reduces exactly to the SCF Hartree–Fock iteration.

### Backends

| Backend | Use Case | Grid | Diagonalizers |
|---------|----------|------|---------------|
| **scf.stdImpl** | General grids | Any | 1: Fock, nPoints×nPoints |
| **scf.ylmOpt** | Spherical atoms | Ylm | (lmax+1): per-l, nRadial×nRadial |
| **mcscf.stdImpl** | Multiconfiguration | Any | 2: CI (nCoeffs) + Fock (nPoints) |
| **mcscf.ylmOpt** | Multiconfiguration, spherical atoms | Ylm | (lmax+2): CI (nCoeffs) + per-l, nRadial×nRadial |

**ylmOpt** exploits angular symmetry: instead of one large matrix, it diagonalizes smaller per-l-channel matrices. This requires the user to set up multiple diagonalizers with l-specific wrappers. In the MCSCF variant the per-l operator is the RDM-weighted Fock operator F̂_l[ρ¹] and the direct potential is the monopole (l=0) projection of the correlated density — exact for spherically symmetric (S) states.

## Public Interface

```fortran
! Bind backend from JSON configuration (also binds the orbital preconditioner)
call GroundSolver_Fabricate()

! Allocate working arrays
call GroundSolver_Setup()

! Perform one SCF/MCSCF iteration
call GroundSolver_Approach(state, time)

! Orbital preconditioner (any backend): dOrb ≈ (F̂ − shift)⁻¹·orb
call GroundSolver_PreconditionAction(dOrb, orb, shift, time)

! One-line note if DiagonalizerList(iDiag) did not converge (backends call it)
call GroundSolver_NoteUnconverged(iDiag)
```

### Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `state` | `complex(R64)(:)` | Packed state (inout); SCF: orbitals, MCSCF: coefficients ⊕ orbitals |
| `time` | `real(R64)` | Time for potentials (usually 0) |

The mixing parameter is not an argument: it is configured in the top-level
`mixing` block (e.g. `mixing.linear.alpha`, see "Mixing Strategies").

### Exported Callbacks

| Pointer | Signature | Bound by | Purpose |
|---------|-----------|----------|---------|
| `GroundSolver_PreconditionAction` | `(dOrb, orb, shift, time)` | `GroundSolver_Fabricate`, every backend | Orbital preconditioner (T̂ + max(−shift, σ_min))⁻¹·orb if `SysKinetic_SolveShifted` is associated, else identity; the `ApplyPrecOnVec` of a lobpcg diagonalizer of any orbital eigenproblem (SCF Fock operator, MCSCF orbital step) |
| `GroundSolver_Scf_HartreeFockAction` | `(dOrb, orb, time)` | scf.stdImpl | Fock action F̂·orb, the `ApplyMatOnVec` of the diagonalizer |
| `GroundSolver_Mcscf_FockAction` | `(dOrb, orb, time)` | mcscf.stdImpl | RDM-weighted Fock action, the `ApplyMatOnVec` of the orbital diagonalizer |

### Unconverged Diagonalizers (all backends)

Every backend calls `GroundSolver_NoteUnconverged(i)` after each
`DiagonalizerList(i)` solve: if the generic `convergedQ` of the diagonalizer
is false (lobpcg at `maxIter` or after its stagnation exit; LAPACK is direct
and ARPACK stops instead), it prints one line, e.g.

```
note: groundSolver: diagonalizerList.lobpcg did not converge (0 of 1 wanted eigenpairs); continuing with its best approximations
```

and the SCF/MCSCF step continues with the returned approximations (all
`nEvals` of them; the next solve warm-starts from them). The SCF convergence
test of the program decides whether the ground state is reached. With
`aceQ: false`, scf.stdImpl follows the first such note with a one-time hint
recommending `aceQ: true` (the exact exchange is self-adjoint in the grid
metric only up to its discretization; see "Preconditioned Eigensolver").

### Orbital Preconditioner (all backends)

`GroundSolver_Fabricate` binds `GroundSolver_PreconditionAction` for every
backend and reads its only parameter:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `groundSolver.precShiftMin` | real | 0.5 | Smallest shift σ_min of the preconditioner (T̂ + σ)⁻¹ (> 0) |

With the Ritz value θ of a residual (`shift`) it applies
(T̂ + max(−θ, σ_min))⁻¹ through the optional kinetic hook
`SysKinetic_SolveShifted` (exact FDM solve on grid.fem; body type 1, i.e. the
mass of the spin-restricted orbitals), and the identity when the kinetic
backend does not provide that hook (Ylm, Linear, Prolate, …). For occupied
orbitals (θ < 0) T̂ − θ approximates F̂ − θ in its stiff high-momentum part;
it is Hermitian positive definite in the grid metric for every shift, as
LOBPCG requires. Use it for orbital vectors (`dim = Grid_nPoints`) only, not
for CI vectors. The hook is looked up on every call, so the fabrication
order of SysKinetic and GroundSolver does not matter; wire it after
`GroundSolver_Fabricate` (the pointer assignment captures the target).

## JSON Configuration

### SCF Standard Implementation (scf.stdImpl)

```json
{
    "groundSolver": {
        "scf": {
            "stdImpl": { }
        }
    },
    "diagonalizerList": {
        "arpack": {
            "nEvals": 5,
            "which": "SR",
            "nKry": 15
        }
    }
}
```

| Parameter      | Type   | Default    | Description |
|----------------|--------|------------|-------------|
| `mixTarget`    | string | "orbitals" | "orbitals", "potential" or "density" (see "Mixing Strategies") |
| `aceQ`         | bool   | false      | Exchange via the adaptively compressed exchange operator |

The preconditioner shift `precShiftMin` is a top-level key,
`groundSolver.precShiftMin` (see "Orbital Preconditioner").

### Adaptively Compressed Exchange (stdImpl `aceQ`)

Without ACE every Fock action costs nUp interaction-potential solves (one
per occupied orbital) for the exchange term. With `aceQ` (Lin 2016) Approach
builds, once per iteration and after the Hartree potential:

- ξ_j = K̂φ_j = ∑_k W(φ_k, φ_j) φ_k for the nUp occupied orbitals, with
  nUp(nUp+1)/2 solves when `SysInteraction_conjSymmetricQ` (the swapped pair
  via `SysInteraction_ConjugateInteractionPotential`), else nUp²
- M_ij = ⟨φ_i|ξ_j⟩ (Hermitian positive definite), M = L Lᴴ (ZPOTRF)
- ζ = ξ L⁻ᴴ

and the Fock action applies Ĵ·orb − K̂_ACE·orb = Ĵ·orb − ∑_k ζ_k ⟨ζ_k|orb⟩ in
one pass (`SysInteraction_AddMeanField`: nUp inner products, threaded).
K̂_ACE equals K̂ on span{φ_j} and K̂ − K̂_ACE is positive
semidefinite, so the SCF fixed point (the lowest nUp eigenvectors) and the
converged energy are unchanged; intermediate iterates differ from the exact-
exchange iteration. All mix targets work as before; `Method_GetEnergy`
always uses the exact exchange. Construction and application are the
grid-agnostic mean-field helpers of SysInteraction
(`SysInteraction_FillExchangeAce` with weighted orbitals = the occupied
orbitals, `SysInteraction_AddMeanField`; src/SysInteraction/A_SysInteraction.md),
shared with the frozen mean-field reference of the Cayley exact part of
Method (`exactFieldFreeMeanFieldQ`, src/Method/A_Method.md).

### Preconditioned Eigensolver (lobpcg)

For large grids whose Fock operator has a huge spectral radius (e.g. finite
elements) use the LOBPCG backend with the grid metric and the orbital
preconditioner, warm-started across SCF iterations:

```fortran
call GroundSolver_Fabricate
DiagonalizerListInput(1) % ApplyMatOnVec => GroundSolver_Scf_HartreeFockAction
DiagonalizerListInput(1) % ApplyPrecOnVec => GroundSolver_PreconditionAction
DiagonalizerListInput(1) % dim = Grid_nPoints
call DiagonalizerList_Fabricate(DiagonalizerListInput)
...
call Grid_Setup
DiagonalizerList(1) % e % metricWeights = Grid_Fem_weights   ! the grid metric, e.g. Grid_Ylm_weights on Ylm
```

```json
"diagonalizerList": { "lobpcg": { "nEvals": 5, "nGuard": 2, "tol": 1e-8 } },
"groundSolver": { "precShiftMin": 0.5, "scf": { "stdImpl": { "mixTarget": "potential", "aceQ": true } } }
```

The metric is the diagonal one of the grid (`Grid_InnerProduct` =
Σ conj(f) g w): assign the grid weights to `metricWeights` after
`Grid_Setup`; left unallocated, LOBPCG works in the Euclidean metric, in
which the Fock operator of a grid with quadrature weights is not self-
adjoint (see A_DiagonalizerList.md, "Metric"). Programs orthonormalize the
start orbitals after `Method_Setup` (`Orbs_Orthonormalize(Orbs_orbs)`;
OrbsInit only normalizes).

**Finite-element grids (grid.fem).** With `sysKinetic.fem.laplacian` the
preconditioner is the exact (T̂ + σ)⁻¹, and LOBPCG needs 7–13 iterations
for a cold start and 0–12 per warm-started SCF step (test/simulationTests/Fem).
Use `aceQ: true` (recommended with lobpcg; scf.stdImpl prints a hint once
when a Fock solve with `aceQ: false` does not converge): the exact exchange
operator of grid.fem is Hermitian only up to
the multipole truncation of the Poisson boundary values
(|⟨φ|K̂g⟩ − conj⟨g|K̂φ⟩| ≈ 5e-4 relative for `multipoleLmax` 6), so with
`aceQ: false` LOBPCG stalls at relative residual floors of 7e-8 … 3e-7,
whereas K̂_ACE = ∑ζζᴴ is Hermitian by construction and needs nUp inner
products instead of nUp Poisson solves per Fock action. Fresh H₂ SCF of
T_Fem_03 with the exact exchange (2 threads): without a stagnation exit 10
of 14 LOBPCG solves ran into maxIter 100 (1772 Fock products, 61 s); with
it they stop 10 iterations into their floors (unconverged, noted by
`GroundSolver_NoteUnconverged`): 604 products, 21 s, the same SCF
trajectory and energy (to 1e-13); with ACE: 189 products, 5.5 s. On the coarse
H2O mesh of T_Fem_04 the leading error mode of the plain SCF (α = 1)
oscillates with a factor ≈ −0.46, which Hartree-potential mixing with α ≈ 0.7
removes (23 → 22 SCF steps; slower modes remain); DIIS on the Hartree
potential is erratic there.

On the first Approach a diagonalizer without previous eigenvectors gets the
current occupied orbitals as start vectors (`SetStartVectors`, a binding of
the diagonalizer contract): LOBPCG starts from them; LAPACK and ARPACK keep
the no-op default, so their results are unchanged. The preconditioner
is the identity unless the kinetic backend binds `SysKinetic_SolveShifted`
(then LOBPCG converges in a grid-independent number of iterations).

### SCF Ylm-Optimized Implementation (scf.ylmOpt)

```json
{
    "groundSolver": {
        "scf": {
            "ylmOpt": { }
        }
    },
    "diagonalizerList": {
        "arpack1": { "nEvals": 2, "which": "SR" },
        "arpack2": { "nEvals": 1, "which": "SR" }
    }
}
```

For ylmOpt, you need one diagonalizer entry per l-channel (l=0, l=1, ...).

### MCSCF Standard Implementation (mcscf.stdImpl)

```json
{
    "groundSolver": {
        "mcscf": {
            "stdImpl": { }
        }
    },
    "diagonalizerList": {
        "arpackCi":  { "nEvals": 1, "which": "SR", "nKry": 30, "tol": 1e-12 },
        "lapackOrb": { "nEvals": 9 }
    }
}
```

ARPACK's `tol` defaults to 0 (machine precision); a looser tolerance like
1e-12 substantially reduces the number of restarts of the CI solve without
affecting the converged MCSCF energy at typical `convThresh` values.

The CI space is large and sparse → iterative ARPACK. The orbital space is small
(Grid_nPoints) and the Fock spectrum contains degenerate shells → dense LAPACK
(assembled with Grid_nPoints Fock matvecs) is typically much faster than ARPACK
at machine-precision tolerance.

The program must fabricate two diagonalizers:

```fortran
DiagonalizerListInput(1) % ApplyMatOnVec => ApplyMatOnVecCi   ! wraps GroundSolver_Mcscf_HamiltonianAction
DiagonalizerListInput(1) % dim = Coeffs_nCoeffs
DiagonalizerListInput(2) % ApplyMatOnVec => ApplyMatOnVecOrb  ! wraps GroundSolver_Mcscf_FockAction
DiagonalizerListInput(2) % dim = Grid_nPoints
! optional, for a lobpcg orbital diagonalizer (orbital vectors only):
DiagonalizerListInput(2) % ApplyPrecOnVec => GroundSolver_PreconditionAction
```

`arpackOrb.nEvals` must equal the number of spatial orbitals per spin
(`Orbs_nOrbsInState/2`). Example application:
`app/ne/groundstate/mcscf9orbs/P_GroundSolverStdImpl.f90`
(neon MCSCF ground state with 9 spatial orbitals per spin).

### MCSCF Ylm-Optimized Implementation (mcscf.ylmOpt)

```json
{
    "groundSolver": {
        "mcscf": {
            "ylmOpt": {
                "mixTarget": "potential",
                "alphaCi": 1.0
            }
        }
    },
    "diagonalizerList": {
        "lapack1":  { },
        "arpack2":  { "nEvals": 2, "which": "SR" },
        "arpack3":  { "nEvals": 1, "which": "SR" }
    }
}
```

The diagonalizer layout is: entry 1 for the CI space (`dim = Coeffs_nCoeffs`,
callback wrapping `GroundSolver_Mcscf_HamiltonianAction`), entries l+2 for
l = 0..lmax over the radial space (`dim = Grid_Ylm_nRadial`, l-specific
wrappers around `GroundSolver_Mcscf_YlmOpt_FockAction`):

```fortran
DiagonalizerListInput(1) % ApplyMatOnVec => GroundSolver_Mcscf_HamiltonianAction
DiagonalizerListInput(1) % dim = Coeffs_nCoeffs
DiagonalizerListInput(2) % ApplyMatOnVec => ApplyMatOnVecL0  ! wraps FockAction with l=0
DiagonalizerListInput(2) % dim = Grid_Ylm_nRadial
DiagonalizerListInput(3) % ApplyMatOnVec => ApplyMatOnVecL1  ! wraps FockAction with l=1
DiagonalizerListInput(3) % dim = Grid_Ylm_nRadial
```

Each l-channel's `nEvals` must cover the highest radial quantum number
nr = n - l - 1 among the orbitals of that channel. Use ARPACK for the radial
channels: the radial Fock action (exchange enters through the
`Grid_Ylm_{Set,Get}LmComponent` projection pair) is not symmetric in the
Euclidean metric that the dense LAPACK solver assumes, so `lapack*` entries
converge to a wrong fixed point there.

## Usage Example

```fortran
program scf_example
  use M_Mixing
  use M_GroundSolver
  use M_Method
  
  real(R64) :: energy_old, energy_new, conv = 1e10_R64
  real(R64) :: conv_thresh = 1e-12_R64
  
  ! ... setup Grid, SysKinetic, SysPotential, etc. ...
  
  call Mixing_Fabricate()        ! mixing parameter: JSON mixing.linear.alpha
  call GroundSolver_Fabricate()
  ! ... fabricate the DiagonalizerList (see above) ...
  call GroundSolver_Setup()
  
  ! SCF convergence loop
  energy_new = 1e10_R64
  do while (abs(conv) > conv_thresh)
    call GroundSolver_Approach(Method_state, time=0.0_R64)
    energy_old = energy_new
    energy_new = Method_GetEnergy(0.0_R64)
    conv = energy_new - energy_old
  end do
end program
```

## Dependencies

### Required Modules

- **Grid**: Spatial discretization, `Grid_nPoints`, `Grid_InnerProduct`
- **SysKinetic**: Kinetic energy operator, `MultiplyWithKineticOp`; optional
  `SysKinetic_SolveShifted` (preconditioner)
- **SysPotential**: External potential, `FillExternalPotential`
- **SysInteraction**: Two-body interaction, `FillInteractionSrc/Potential`,
  `ConjugateInteractionPotential` (ACE with `conjSymmetricQ`)
- **Orbs**: Orbital container, `Orbs_Orthonormalize`
- **DiagonalizerList**: Iterative eigensolvers (ARPACK, LOBPCG)

### MCSCF-Specific Dependencies

- **Coeffs**: CI coefficients, `Coeffs_ApplyH1FillRdm1/ApplyH2FillRdm2`,
  `Coeffs_Normalize`, `Coeffs_nCoeffs`
- **Method_Mb_OrbBased**: `FillH1`, `FillH2`

### Ylm-Specific Dependencies

- **Grid_Ylm**: `nRadial`, `lmax`, `SetLmComponent`, `GetLmComponent`
- **SysKinetic_Ylm**: `MultiplyWithRadialKineticOp`
- **SysPotential_Ylm**: `FillExternalPotentialRadial`
- **SysInteraction_Ylm**: `FillInteractionPotentialRadial`
- **OrbsInit_Ylm_HydrogenLike**: Quantum numbers `n(j)`, `l(j)`, `m(j)`

## Extending the Module

### Adding a New Backend

1. Create `M_GroundSolver_NewBackend.f90` with interface + fabricate declaration
2. Create `S_GroundSolver_NewBackend.f90` with implementation
3. Add dispatch branch in `S_GroundSolver.f90`:

```fortran
if (Json_GetExistence("groundSolver.newBackend")) then
  call GroundSolver_NewBackend_Fabricate
```

### Implementation Checklist

Your implementation must bind these procedure pointers:
- `GroundSolver_Setup` → Allocate working arrays
- `GroundSolver_Approach` → One SCF iteration
- (Optional) Backend-specific action pointer(s) for the diagonalizer(s)

`GroundSolver_PreconditionAction` is bound by `GroundSolver_Fabricate`
itself, before the dispatch, so every backend gets it.

## Tests

| Test | Description |
|------|-------------|
| `T_Ne3d_02_GroundSolverScf` | Neon atom (scf.stdImpl), E ≈ −127.74 Ha |
| `T_Ne3d_03_GroundSolverScfYlm` | Neon atom (scf.ylmOpt), E ≈ −128.55 Ha |
| `T_Ne3d_04_GroundSolverMcscf` | Neon atom (mcscf.stdImpl), E ≈ −127.74 Ha |
| `T_Ne3d_05_GroundSolverScfDensity` | Neon atom (scf.stdImpl, `mixTarget: "density"`), same energy as T_Ne3d_02 |
| `T_Ne3d_06_GroundSolverMcscfYlm` | Neon atom (mcscf.ylmOpt), same grid/energy as T_Ne3d_04 |
| `T_Ne3d_17_GroundSolverScfLobpcgAce` | Neon atom (scf.stdImpl with `aceQ`, lobpcg diagonalizer in the grid metric, `GroundSolver_PreconditionAction` = identity on Ylm), same energy as T_Ne3d_02; asserts that every LOBPCG solve converges (no maxIter, no stagnation exit) and at most 30 SCF steps (measured 14; the margin absorbs round-off in the step count: the 1e-12 energy criterion sits at the noise floor of the tol-1e-8 solves) |
| `T_Fem_03_H2ScfRhf` | H₂ RHF on grid.fem (scf.stdImpl, `aceQ`, lobpcg with `GroundSolver_PreconditionAction` = exact (T̂ + σ)⁻¹ and `metricWeights`): E_el = −1.84795386088 (prototype on the same mesh to 2e-11), ACE and exact exchange agree to 1e-13; checks the preconditioner contract; writes `orb01_01.in` |
| `T_Fem_04_H2oScfRhf` | H₂O RHF on a coarse grid.fem mesh (5 orbitals, STO guess): E_tot = −76.40375596 (prototype on the same mesh to 5e-11; 0.34 Ha below the HF limit, Galerkin nuclei are not variational) |

Run with: `ctest -R "T_Ne3d_0[2-6]|T_Ne3d_17|T_Fem_0[34]"`

## Physical Background

The SCF method solves the **Hartree–Fock equations**:

```
F̂ φ_i = ε_i φ_i
```

where the Fock operator is:

```
F̂ = T̂ + V̂_ext + Ĵ − K̂
```

- **T̂**: Kinetic energy (−½∇²)
- **V̂_ext**: External potential (e.g., −Z/r for atoms)
- **Ĵ**: Hartree (direct) potential, electron repulsion mean-field
- **K̂**: Fock (exchange) operator, antisymmetry requirement

The self-consistency arises because Ĵ and K̂ depend on the orbitals {φ_i} being solved for.

The MCSCF method additionally expands the state in configurations |Φ_I⟩ built
from the orbitals and makes both the CI vector and the orbitals self-consistent:
the CI Hamiltonian depends on the orbitals through h1/h2, and the RDM-weighted
Fock operator depends on the CI vector through ρ¹.

## Known Limitations

1. **Restricted calculations only**: Spin-up and spin-down orbitals are assumed identical
2. **ylmOpt assumes m=0 trials**: The Fock action uses m=0 in the diagonalizer callbacks
3. **MCSCF orbital step is Fock-canonical**: The orbital update takes the lowest
   eigenvectors of F̂[ρ¹]; this is a robust heuristic whose fixed point can
   deviate slightly from the fully variational MCSCF orbitals
4. **Declared m** (`method.mb.orbBased.mOfOrbs`, src/Method/A_Method.md): only
   `scf.ylmOpt` and `mcscf.ylmOpt` are accepted (checked in
   `GroundSolver_Fabricate`). They place every radial eigenvector in its
   (l, m) channel and align with the old orbitals per m block
   (`Orbs_AlignOnReference(..., blockOf_=Orbs_blockOfOrb(1:nUp))`), so the
   channels of other m stay exactly zero; the `stdImpl` solvers diagonalize on
   the full grid and mix degenerate orbitals of different m.

## Mixing Strategies

Iterate updates are delegated to the Mixing module (see `src/Mixing/A_Mixing.md`):
the algorithm (`mixing.linear` / `mixing.diis`) is a top-level JSON choice (no
block ⇒ linear), and `groundSolver.scf.stdImpl.mixTarget` selects what gets
mixed:

- `"orbitals"` (default): mix the gauge-aligned eigenvectors into the old
  orbitals, then orthonormalize
- `"potential"`: mix the Hartree potential (≡ density mixing) and replace the
  orbitals outright — often converges in far fewer iterations, since the
  potential is gauge-invariant and the map is smooth in it
- `"density"`: purify the mixed 1-RDM (1−α)·D_old + α·D_new to its top
  natural orbitals (`Orbs_MixOccupiedSpace`, α from `mixing.linear.alpha`)

Rule of thumb: `potential` + `linear` with α ≈ 0.5–0.9 is a strong default;
add DIIS for systems where linear mixing stalls. Avoid DIIS with very small α
(the exchange operator lags the mixed potential; see A_Mixing.md).

## Troubleshooting

| Issue | Solution |
|-------|----------|
| Slow/no convergence | Try `mixTarget: "potential"`; reduce α (e.g., 0.3–0.5); try `mixing.diis` |
| ARPACK convergence failure | Increase `nKry`, check `nEvals` ≤ dim/2 |
| ARPACK hopeless on a huge grid (stiff Fock operator) | Use `lobpcg` with `GroundSolver_PreconditionAction`, `metricWeights` (the grid weights) and `aceQ` |
| LOBPCG hits maxIter | Without `SysKinetic_SolveShifted` the preconditioner is the identity: raise `maxIter` / relax `tol` |
| "LOBPCG ... stagnates" warnings and "did not converge" notes in every SCF step | The Fock operator is not self-adjoint in the grid metric to the requested `tol`; on grid.fem with `aceQ: false` (scf.stdImpl adds a one-time hint): set `aceQ: true` |
| Wrong energy (ylmOpt) | Verify DiagonalizerList order matches l=0,1,2,... |
| Degenerate Fock levels (MCSCF) | Use a `lapack*` orbital diagonalizer; ensure `nEvals` covers full shells |
