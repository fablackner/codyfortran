# SysKinetic Module — Onboarding Guide

## Overview

The **SysKinetic** module provides the kinetic energy operator T̂ for quantum many-body simulations. It implements:

```
T̂ = −(ℏ²/2m) ∇²
```

In atomic units (ℏ = m_e = 1), this simplifies to T̂ = −(1/2m)∇² for a particle of mass m.

## Architecture

The module follows the **CodyFortran fabrication pattern**:

```
┌─────────────────────────────────────────────────────────────────┐
│                    M_SysKinetic (Facade)                        │
│  - SysKinetic_Fabricate          → runtime configuration        │
│  - SysKinetic_Setup              → precompute data (pointer)    │
│  - SysKinetic_MultiplyWithKineticOp → apply T̂·ψ (pointer)       │
│  - SysKinetic_SolveShifted       → (T̂ + s)⁻¹ (optional hook)    │
├─────────────────────────────────────────────────────────────────┤
│                          │                                      │
│    ┌────────────────┬────┴────────────┬───────────────┐         │
│    ▼                ▼                 ▼               ▼         │
│  Lattice         Linear             Ylm             Fem         │
│ (hopping)    (Cartesian ∇²)   (spherical radial) (3D spectral   │
│    │                │                 │           elements)     │
│    ▼                ▼                 ▼               ▼         │
│ NearestNeighbor  Laplacian        Laplacian       Laplacian     │
│                  ├─FinDiff        ├─FinDiff       (+ exact      │
│                  └─Fourier        └─Fedvr          (T̂+s)⁻¹)     │
└─────────────────────────────────────────────────────────────────┘
```

### Key Design Principles

1. **Interface/Implementation Separation**: `M_*.f90` files define contracts (procedure pointers), `S_*.f90` files provide implementations.

2. **Runtime Binding**: Procedure pointers are assigned by `*_Fabricate()` based on JSON configuration. This allows swapping implementations without recompilation.

3. **Hierarchical Dispatch**: The top-level `SysKinetic_Fabricate` reads JSON and delegates to the appropriate backend (Lattice/Linear/Ylm/Prolate/Fem), which further delegates to specific implementations.

4. **Optional hooks**: `SysKinetic_SolveShifted` (exact solve of (T̂ + s) x = b) is `null()` unless the backend has a fast exact inverse (Fem). `SysKinetic_Fabricate` nullifies it before it dispatches, so a backend fabricated earlier in the same process leaves no stale binding. Callers test `associated(SysKinetic_SolveShifted)` and fall back otherwise.

## File Structure

```
SysKinetic/
├── M_SysKinetic.f90          # Public facade: exports pointers and Fabricate
├── S_SysKinetic.f90          # Dispatch to Lattice/Linear/Ylm/Prolate/Fem
├── A_SysKinetic.md           # This file
│
├── Fem/                      # 3D Cartesian spectral elements (grid.fem)
│   ├── M_SysKinetic_Fem.f90
│   ├── S_SysKinetic_Fem.f90  # checks grid.fem, dispatches to laplacian
│   └── Laplacian/            # T = −lap/(2m) and the exact (T + s)⁻¹ (FDM)
│       ├── M_SysKinetic_Fem_Laplacian.f90
│       └── S_SysKinetic_Fem_Laplacian.f90
│
├── Lattice/                  # Tight-binding on discrete lattice
│   ├── M_SysKinetic_Lattice.f90
│   ├── S_SysKinetic_Lattice.f90
│   └── NearestNeighbor/      # −t Σ (|i⟩⟨j| + h.c.) hopping
│       ├── M_SysKinetic_Linear_NearestNeighbor.f90
│       └── S_SysKinetic_Linear_NearestNeighbor.f90
│
├── Linear/                   # 1D Cartesian uniform grid
│   ├── M_SysKinetic_Linear.f90
│   ├── S_SysKinetic_Linear.f90
│   └── Laplacian/
│       ├── M_SysKinetic_Linear_Laplacian.f90
│       ├── S_SysKinetic_Linear_Laplacian.f90
│       ├── FinDiff/          # Central finite-difference stencil
│       │   ├── M_SysKinetic_Linear_Laplacian_FinDiff.f90
│       │   └── S_SysKinetic_Linear_Laplacian_FinDiff.f90
│       └── Fourier/          # FFT-based spectral method
│           ├── M_SysKinetic_Linear_Laplacian_Fourier.f90
│           └── S_SysKinetic_Linear_Laplacian_Fourier.f90
│
└── Ylm/                      # Spherical harmonics (radial + angular)
    ├── M_SysKinetic_Ylm.f90
    ├── S_SysKinetic_Ylm.f90
    └── Laplacian/
        ├── M_SysKinetic_Ylm_Laplacian.f90
        ├── S_SysKinetic_Ylm_Laplacian.f90
        ├── FinDiff/          # FD on uniform radial grid
        │   ├── M_SysKinetic_Ylm_Laplacian_FinDiff.f90
        │   └── S_SysKinetic_Ylm_Laplacian_FinDiff.f90
        ├── Fedvr/            # Finite-element DVR on nonuniform grid
        │   ├── M_SysKinetic_Ylm_Laplacian_Fedvr.f90
        │   └── S_SysKinetic_Ylm_Laplacian_Fedvr.f90
        └── FedvrEcs/         # FEDVR on exterior-complex-scaled contour
            ├── M_SysKinetic_Ylm_Laplacian_FedvrEcs.f90
            └── S_SysKinetic_Ylm_Laplacian_FedvrEcs.f90
```

## Backend Summary

| Backend | Grid Type | Operator | Boundary Conditions |
|---------|-----------|----------|---------------------|
| **Lattice/NearestNeighbor** | 3D discrete sites | Tight-binding hopping | Periodic or hard-wall |
| **Linear/Laplacian/FinDiff** | 1D uniform Cartesian | −(1/2m) d²/dx² | Dirichlet (zero ends) |
| **Linear/Laplacian/Fourier** | 1D uniform Cartesian | −(1/2m) d²/dx² via FFT | Periodic |
| **Ylm/Laplacian/FinDiff** | Radial uniform | −(1/2m)[∇²_r − l(l+1)/r²] | f(0)=0, f(r_max)=0 |
| **Ylm/Laplacian/Fedvr** | Radial nonuniform | −(1/2m)[∇²_r − l(l+1)/r²] | f(0)=0, f(r_max)=0 |
| **Ylm/Laplacian/FedvrEcs** | Radial nonuniform, ECS contour | −(1/2m)[∇²_z − l(l+1)/z²] on contour z(r) | f(0)=0, absorbing beyond ECS radius |
| **Fem/Laplacian** | 3D Cartesian tensor mesh (GLL spectral elements, graded) | −(1/2m) lap, lap = −W⁻¹K3 (nodal strong form); exact (T̂ + s)⁻¹ | Dirichlet on the box surface |

## Usage

### Initialization Sequence

```fortran
! 1. Read JSON configuration and bind procedure pointers
call SysKinetic_Fabricate

! 2. Precompute grid-dependent data (stencils, FFT plans, etc.)
call SysKinetic_Setup

! 3. Apply the kinetic operator during time evolution
call SysKinetic_MultiplyWithKineticOp(dOrb, orb, time, bt_)
```

### JSON Configuration Examples

**Linear grid with finite-difference Laplacian:**
```json
{
  "sysKinetic": {
    "linear": {
      "laplacian": {
        "bodyMass": [1.0],
        "finDiff": {}
      }
    }
  }
}
```

**Lattice with nearest-neighbor hopping:**
```json
{
  "sysKinetic": {
    "lattice": {
      "nearestNeighbor": {
        "hoppX": 1.0,
        "hoppY": 1.0,
        "hoppZ": 0.5
      }
    }
  }
}
```

**Ylm (atomic) with FEDVR radial grid:**
```json
{
  "sysKinetic": {
    "ylm": {
      "laplacian": {
        "bodyMass": [1.0],
        "fedvr": {}
      }
    }
  }
}
```

**Fem (molecules, requires `grid.fem`):**
```json
{
  "sysKinetic": {
    "fem": {
      "laplacian": {
        "bodyMass": [1.0, 1.0]
      }
    }
  }
}
```

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `sysKinetic.fem.laplacian.bodyMass` | real(:) | [1.0] | Masses, every entry > 0: a single entry is the mass of every body type (body-type independent); several entries give one mass per body type (e.g. `[1.0, 1836.0]`), and a body type beyond them stops. An absent `bt_` uses `bodyMass(1)`. |

## Exported Symbols

### Module Data (M_SysKinetic)

| Symbol | Type | Description |
|--------|------|-------------|
| `SysKinetic_timeIndependentQ` | `logical` | `.true.` if T̂ is time-independent (all current backends) |
| `SysKinetic_bodyTypeIndependentQ` | `logical` | `.true.` if same mass for all body types |

### Procedure Pointers (M_SysKinetic)

| Pointer | Signature | Description |
|---------|-----------|-------------|
| `SysKinetic_Setup` | `subroutine()` | Initialize grid-dependent data |
| `SysKinetic_MultiplyWithKineticOp` | `subroutine(dOrb, orb, time, bt_)` | Apply T̂·ψ → dψ |
| `SysKinetic_SolveShifted` | `subroutine(dOrb, orb, shift, bt_)`, `complex(R64) :: shift` | Optional (default `null()`): solve (T̂ + shift) dOrb = orb exactly; only grids with a fast exact kinetic inverse bind it (Fem). Real shift: preconditioner (T̂ + σ)⁻¹; shift = 2/(i dt) = −2i/dt: Crank–Nicolson factor 1 + i dt T̂/2 = (i dt/2)(T̂ − 2i/dt) |

### Fem Laplacian (M_SysKinetic_Fem_Laplacian)

| Symbol | Type | Description |
|--------|------|-------------|
| `SysKinetic_Fem_Laplacian_bodyMass` | `real(R64)(:)` | Mass per body type |

### Radial Operator (M_SysKinetic_Ylm)

| Pointer | Signature | Description |
|---------|-----------|-------------|
| `SysKinetic_Ylm_MultiplyWithRadialKineticOp` | `subroutine(dOrbLm, orbLm, l, m, time, bt_)` | Apply T̂_{lm} per channel |

## Physics Notes

### Ylm Radial Kinetic Operator

The spherical Laplacian separates as:

```
∇² = (1/r²) ∂/∂r (r² ∂/∂r) − L̂²/r²
```

For a single (l,m) channel, L̂² → l(l+1), giving:

```
T̂_{lm} f(r) = −(1/2m) [ d²f/dr² + (2/r) df/dr − l(l+1)/r² f ]
```

**Implementation trick**: Transform g(r) = r·f(r) to convert the radial Laplacian to a simple second derivative:

```
[ d²/dr² + (2/r) d/dr ] f = (1/r) d²g/dr²
```

This avoids the 2/r singularity in finite-difference stencils.

### Fem Spectral-Element Laplacian

On `grid.fem` orbitals are nodal values on the interior GLL nodes of the
tensor mesh (see `src/Grid/A_Grid.md`). With the lumped (diagonal) mass W and
the stiffness K3 the Galerkin Laplacian in nodal (strong) form is
lap = −W⁻¹K3 = −Σ_a W_a⁻¹K_a, applied as three banded stencils
(`Grid_Fem_ApplyLaplacian`), and

```
T̂ f = −lap f / (2 m_bt)
```

is self-adjoint in the metric of `Grid_InnerProduct` (⟨f|T̂g⟩ = conj⟨g|T̂f⟩ to
1e-14), so h1 = ⟨φ_i|T̂φ_j⟩ is the Galerkin matrix element.

**Exact shifted solve.** The fast diagonalization of the separable Laplacian
(−lap = S3 D S3ᵀ W3, D(i,j,k) = λx(i) + λy(j) + λz(k)) gives

```
(T̂ + s)⁻¹ f = 2m (−lap + 2m s)⁻¹ f = 2m · S3 (D + 2m s)⁻¹ S3ᵀ W3 f
```

for a real or complex shift s (`Grid_Fem_SolveShiftedLaplacian`, generic in
the type of the shift). The mode products stay real (real and imaginary parts
transformed separately); only the diagonal (D + 2m s)⁻¹ is complex for a
complex shift, where the two parts mix. A real shift must exceed −(smallest
eigenvalue of T̂) (stop otherwise); a shift with nonzero imaginary part is
always invertible. Uses: the preconditioner (T̂ + σ)⁻¹ of the orbital
eigenproblems (GroundSolver) and the preconditioner of the Crank–Nicolson /
Cayley step of the exact field-free part h0 = T̂ + V_static (shift −2i/dt,
`S_Method_Mb_OrbBased_FieldFreeCayley`; its presence selects the Cayley
exact part on non-Ylm grids).

**Flags and masses.** `SysKinetic_timeIndependentQ = .true.`;
`SysKinetic_bodyTypeIndependentQ = .true.` only if all masses are equal,
always for a single `bodyMass` entry, which then serves every body type
(the default `[1.0]` works for any number of body types). With several
entries a body type beyond them stops, also when they are all equal.

**Threading.** Both operators are reentrant (local work arrays only) and are
called by Method inside OpenMP loops over orbitals; their own OpenMP loops
carry `if (.not. omp_in_parallel())`.

**Cost** (2 threads, N = 328509, complex fields): MultiplyWithKineticOp
≈ 8 ms, SolveShifted ≈ 38 ms (real shift) / 41 ms (complex shift).

**Test:** `test/componentTests/SysKinetic/Fem/Laplacian/T_SysKinetic_Fem_Laplacian`
(box ±8, nLocals 6, one center, hMin 0.2, growth 1.6, N = 328509):
Hermiticity 4e-14; ⟨g|T̂|g⟩/⟨g|g⟩ for g = e^(−a r²): 1.4999990984 (a = 1,
exact 1.5), 0.7500001450 (a = 0.5); mass scaling 1.5e-16; (T̂ + s)x − b
relative residual 3.2e-12 (s = 0.7), 1.3e-13 (s = −2i/0.05), 7.5e-14
(mass 2.5); refabricated with the single mass `[1.7]`: body-type independent,
T̂ for bt_ = 3 and without bt_ equal to bt_ = 1 (difference 0),
⟨g|T̂|g⟩/⟨g|g⟩ = 0.8823524108 (exact 1.5/1.7), (T̂ + 0.7)x − b for
bt_ = 3: 2.8e-12; `[1.0, 1.0]`: body-type independent.
`T_SysKinetic_Fem_Laplacian_MissingMass` (expected error) checks that T̂ for
body type 3 with `bodyMass: [1.0, 2.5]` stops with "no mass for body type 3".

### Lattice Tight-Binding

The tight-binding kinetic operator:

```
T̂ = −Σ_{⟨i,j⟩} t_{ij} (ĉ†_i ĉ_j + h.c.)
```

approximates the continuum Laplacian for hopping t = 1/(2ma²) where a is the lattice spacing.

## Adding a New Backend

1. Create a new directory under the appropriate branch (e.g., `Linear/Laplacian/MyMethod/`)

2. Create `M_SysKinetic_*_MyMethod.f90` with:
   - Interface for `*_Fabricate` subroutine
   - Any backend-specific module data

3. Create `S_SysKinetic_*_MyMethod.f90` with:
   - `*_Fabricate`: read JSON, bind pointers, validate requirements
   - `Setup`: precompute any grid-dependent data
   - `MultiplyWithKineticOp`: the actual operator implementation

4. Add `use M_SysKinetic_*_MyMethod` to parent fabricator

5. Add JSON branch detection in parent's `*_Fabricate`

6. Add corresponding test in `test/` directory

## Dependencies

- `M_Utils_Types`: Type definitions (I32, R64)
- `M_Utils_Json`: JSON configuration reader
- `M_Utils_Say`: Logging/output
- `M_Grid`, `M_Grid_*`: Grid definitions and utilities
- `M_Utils_DerivativeFinDiff`: FD stencil implementation
- `M_Utils_DerivativeFftw`: FFT-based derivatives
- `M_Utils_DerivativeFedvr`: FEDVR derivatives
- `M_Grid_Fem`: nodal Laplacian and fast-diagonalization shifted solves (Fem)
