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

The module follows the **CodyFortranRDM fabrication pattern**:
- **M_* files**: Declare procedure pointers and abstract interfaces (contracts)
- **S_* files**: Provide implementations and JSON-based dispatch logic

## Key Concepts

### SCF Iteration

Each `GroundSolver_Approach` call performs one SCF cycle:

1. **Build mean-field potentials** from current orbitals
   - Hartree (direct): Ĵ = ∑_j ∫|φ_j|²/|r−r'| dr'
   - External: V̂_ext (e.g., nuclear Coulomb −Z/r)

2. **Diagonalize the Fock operator**
   - F̂ = T̂ + V̂_ext + Ĵ − K̂
   - Uses iterative ARPACK solver via DiagonalizerList

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
! Bind backend from JSON configuration
call GroundSolver_Fabricate()

! Allocate working arrays
call GroundSolver_Setup()

! Perform one SCF/MCSCF iteration
call GroundSolver_Approach(state, alpha, time)
```

### Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `state` | `complex(R64)(:)` | Packed state (inout); SCF: orbitals, MCSCF: coefficients ⊕ orbitals |
| `alpha` | `real(R64)` | Mixing parameter (0 < α ≤ 1) |
| `time` | `real(R64)` | Time for potentials (usually 0) |

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
  
  real(R64) :: energy_old, energy_new, conv
  real(R64) :: alpha = 0.5_R64
  real(R64) :: conv_thresh = 1e-12_R64
  
  ! ... setup Grid, SysKinetic, SysPotential, etc. ...
  
  call Mixing_Fabricate()
  call GroundSolver_Fabricate()
  call GroundSolver_Setup()
  
  ! SCF convergence loop
  energy_new = 1e10_R64
  do while (abs(conv) > conv_thresh)
    call GroundSolver_Approach(Method_state, alpha, time=0.0_R64)
    energy_old = energy_new
    energy_new = Method_GetEnergy(0.0_R64)
    conv = energy_new - energy_old
  end do
end program
```

## Dependencies

### Required Modules

- **Grid**: Spatial discretization, `Grid_nPoints`
- **SysKinetic**: Kinetic energy operator, `MultiplyWithKineticOp`
- **SysPotential**: External potential, `FillExternalPotential`
- **SysInteraction**: Two-body interaction, `FillInteractionSrc/Potential`
- **Orbs**: Orbital container, `Orbs_Orthonormalize`
- **DiagonalizerList**: Iterative eigensolvers (ARPACK)

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

## Tests

| Test | Description |
|------|-------------|
| `T_Ne3d_02_GroundSolverScf` | Neon atom (scf.stdImpl), E ≈ −127.74 Ha |
| `T_Ne3d_03_GroundSolverScfYlm` | Neon atom (scf.ylmOpt), E ≈ −128.55 Ha |
| `T_Ne3d_04_GroundSolverMcscf` | Neon atom (mcscf.stdImpl), E ≈ −127.74 Ha |
| `T_Ne3d_06_GroundSolverMcscfYlm` | Neon atom (mcscf.ylmOpt), same grid/energy as T_Ne3d_04 |

Run with: `ctest -R T_Ne3d_0[2346]`

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

Rule of thumb: `potential` + `linear` with α ≈ 0.5–0.9 is a strong default;
add DIIS for systems where linear mixing stalls. Avoid DIIS with very small α
(the exchange operator lags the mixed potential; see A_Mixing.md).

## Troubleshooting

| Issue | Solution |
|-------|----------|
| Slow/no convergence | Try `mixTarget: "potential"`; reduce α (e.g., 0.3–0.5); try `mixing.diis` |
| ARPACK convergence failure | Increase `nKry`, check `nEvals` ≤ dim/2 |
| Wrong energy (ylmOpt) | Verify DiagonalizerList order matches l=0,1,2,... |
| Degenerate Fock levels (MCSCF) | Use a `lapack*` orbital diagonalizer; ensure `nEvals` covers full shells |
