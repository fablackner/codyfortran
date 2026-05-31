# SysKinetic Module — Onboarding Guide

## Overview

The **SysKinetic** module provides the kinetic energy operator T̂ for quantum many-body simulations. It implements:

```
T̂ = −(ℏ²/2m) ∇²
```

In atomic units (ℏ = m_e = 1), this simplifies to T̂ = −(1/2m)∇² for a particle of mass m.

## Architecture

The module follows the **CodyFortranRDM fabrication pattern**:

```
┌─────────────────────────────────────────────────────────────────┐
│                    M_SysKinetic (Facade)                        │
│  - SysKinetic_Fabricate          → runtime configuration        │
│  - SysKinetic_Setup              → precompute data (pointer)    │
│  - SysKinetic_MultiplyWithKineticOp → apply T̂·ψ (pointer)       │
├─────────────────────────────────────────────────────────────────┤
│                          │                                      │
│    ┌────────────────┬────┴────────────┬───────────────┐         │
│    ▼                ▼                 ▼               │         │
│  Lattice         Linear             Ylm              │         │
│ (hopping)    (Cartesian ∇²)   (spherical radial)    │         │
│    │                │                 │              │         │
│    ▼                ▼                 ▼              │         │
│ NearestNeighbor  Laplacian        Laplacian         │         │
│                  ├─FinDiff        ├─FinDiff         │         │
│                  └─Fourier        └─Fedvr           │         │
└─────────────────────────────────────────────────────────────────┘
```

### Key Design Principles

1. **Interface/Implementation Separation**: `M_*.f90` files define contracts (procedure pointers), `S_*.f90` files provide implementations.

2. **Runtime Binding**: Procedure pointers are assigned by `*_Fabricate()` based on JSON configuration. This allows swapping implementations without recompilation.

3. **Hierarchical Dispatch**: The top-level `SysKinetic_Fabricate` reads JSON and delegates to the appropriate backend (Lattice/Linear/Ylm), which further delegates to specific implementations.

## File Structure

```
SysKinetic/
├── M_SysKinetic.f90          # Public facade: exports pointers and Fabricate
├── S_SysKinetic.f90          # Dispatch to Lattice/Linear/Ylm
├── AGENTS.md                 # This file
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
        └── Fedvr/            # Finite-element DVR on nonuniform grid
            ├── M_SysKinetic_Ylm_Laplacian_Fedvr.f90
            └── S_SysKinetic_Ylm_Laplacian_Fedvr.f90
```

## Backend Summary

| Backend | Grid Type | Operator | Boundary Conditions |
|---------|-----------|----------|---------------------|
| **Lattice/NearestNeighbor** | 3D discrete sites | Tight-binding hopping | Periodic or hard-wall |
| **Linear/Laplacian/FinDiff** | 1D uniform Cartesian | −(1/2m) d²/dx² | Dirichlet (zero ends) |
| **Linear/Laplacian/Fourier** | 1D uniform Cartesian | −(1/2m) d²/dx² via FFT | Periodic |
| **Ylm/Laplacian/FinDiff** | Radial uniform | −(1/2m)[∇²_r − l(l+1)/r²] | f(0)=0, f(r_max)=0 |
| **Ylm/Laplacian/Fedvr** | Radial nonuniform | −(1/2m)[∇²_r − l(l+1)/r²] | f(0)=0, f(r_max)=0 |

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
