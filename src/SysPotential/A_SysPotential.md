# SysPotential Module — Agent Onboarding Guide

## Overview

The **SysPotential** module implements external (one-body) potentials for the
CodyFortranRDM quantum many-body simulation framework. It provides the operator
**V̂_ext** in the many-body Hamiltonian:

```
Ĥ = T̂ + V̂_ext + Ŵ
```

where `V̂_ext = Σᵢ V(rᵢ,t)` describes trapping potentials, nuclear attraction,
disorder fields, or time-dependent driving.

---

## Architecture

### Module/Submodule Pattern

This module follows the framework's **Interface Module + Implementation Submodule**
paradigm:

```
src/SysPotential/
├── M_SysPotential.f90              ← Public interface (abstract contracts)
├── S_SysPotential.f90              ← Top-level fabrication (runtime dispatch)
├── AGENTS.md                       ← This file
│
├── Lattice/                        ← 3D discrete lattice backends
│   ├── M_SysPotential_Lattice.f90
│   ├── S_SysPotential_Lattice.f90
│   ├── Harmonic/                   ← Parabolic trap V = ½ω²r²
│   ├── RandomUniform/              ← Uniform disorder V ∈ [a,b]
│   ├── RandomGauss/                ← Gaussian disorder V ~ N(μ,σ²)
│   ├── Seesaw/                     ← Time-dependent linear tilt
│   └── Manual/                     ← User-specified site values
│
├── Linear/                         ← 1D real-space backends
│   ├── M_SysPotential_Linear.f90
│   ├── S_SysPotential_Linear.f90
│   ├── Harmonic/                   ← 1D harmonic trap
│   └── SoftYukawa/                 ← Multi-center softened Yukawa
│
└── Ylm/                            ← Spherical harmonics backends
    ├── M_SysPotential_Ylm.f90
    ├── S_SysPotential_Ylm.f90
    └── Coulomb/                    ← Central Coulomb -Z/r
```

**Key rule:** Interface modules (`M_*.f90`) only import other `M_*` modules.
Implementation submodules (`S_*.f90`) contain the actual algorithms.

### Procedure Pointer Binding

At runtime, `SysPotential_Fabricate` reads the JSON configuration and binds
procedure pointers to concrete implementations. This enables:

- Swapping potentials without recompilation
- JSON-driven physics configuration
- Clean separation of interface from implementation

---

## Core API

### Public Procedure Pointers

| Pointer | Signature | Purpose |
|---------|-----------|---------|
| `SysPotential_Setup` | `()` | One-time initialization (precompute caches) |
| `SysPotential_FillExternalPotential` | `(pot, t, bt)` | Compute V(r,t) on the grid |
| `SysPotential_MultiplyWithExternalPotential` | `(dOrb, pot, orb)` | Apply dOrb = V · orb |

### Public Flags

| Variable | Type | Meaning |
|----------|------|---------|
| `SysPotential_timeIndependentQ` | `logical` | True if V̂ is time-independent (enables caching) |
| `SysPotential_bodyTypeIndependentQ` | `logical` | True if V̂ treats all species identically |

### Ylm-Specific API

| Pointer/Variable | Purpose |
|------------------|---------|
| `SysPotential_Ylm_FillExternalPotentialRadial` | Get V_{lm}(r) component |
| `SysPotential_Ylm_lmax` | Angular expansion cutoff |
| `SysPotential_Ylm_mIndependentQ` | True if potential is azimuthally symmetric |

---

## Grid Back-Ends

### 1. Lattice (3D Discrete)

**Use case:** Hubbard models, optical lattices, disordered systems

**Implementations:**

| Type | Formula | Time-dep? | Parameters |
|------|---------|-----------|------------|
| Harmonic | V = ½(ω_x(x-x₀)² + ω_y(y-y₀)² + ω_z(z-z₀)²) | No | omegaX/Y/Z, positionX/Y/Z |
| RandomUniform | V(i) ~ U[a, b] | No | minValue, maxValue, seed, sites |
| RandomGauss | V(i) ~ N(μ, σ²) | No | meanValue, stdValue, seed |
| Seesaw | V = slope(t) · (r - r_center) | **Yes** | slopeMin/Max, frequency per axis |
| Manual | V(i) = user values | No | sites[], values[] |

**JSON path:** `sysPotential.lattice.<type>.stdImpl`

### 2. Linear (1D Real-Space)

**Use case:** 1D model systems, trapped cold atoms, molecular chains

**Implementations:**

| Type | Formula | Time-dep? | Parameters |
|------|---------|-----------|------------|
| Harmonic | V = ½ Σₙ ωₙ²(x - xₙ)² | No | position[], omega[] |
| SoftYukawa | V = -Σₙ qₙ e^{-αr}/(√(r²+s₁²)+s₂) | No | position[], charge[], softening1/2[], dampening[] |

**Soft-Yukawa kernel:**
```
V(x) = -q × exp(-α|x-x₀|) / (√((x-x₀)² + s₁²) + s₂)
```
where `s₁` regularizes the Coulomb singularity and `s₂` provides additional smoothing.

**JSON path:** `sysPotential.linear.<type>.stdImpl`

### 3. Ylm (Spherical Harmonics)

**Use case:** Atoms, ions, spherically symmetric systems

**Implementations:**

| Type | Formula | Time-dep? | Parameters |
|------|---------|-----------|------------|
| Coulomb | V(r) = -Z/r (l=0, m=0 only) | No | charge (Z) |

**Expansion:**
The potential is stored as V_{lm}(r) coefficients where:
```
V(r,Ω) = Σ_{l,m} V_{lm}(r) Y_{lm}(Ω)
```

For a central Coulomb potential, only the (l=0, m=0) term is nonzero:
```
V_{00}(r) = -Z × √(4π) / r
```

**JSON path:** `sysPotential.ylm.coulomb.stdImpl`

---

## Workflow

### Initialization Sequence

```fortran
call SysPotential_Fabricate    ! Parse JSON, bind pointers
call SysPotential_Setup        ! Precompute any cached arrays
```

### Typical Usage (in Method/Hamiltonian)

```fortran
! Compute external potential (once if time-independent)
call SysPotential_FillExternalPotential(extPot, time)

! Apply to orbital:  dOrb = V · orb
call SysPotential_MultiplyWithExternalPotential(dOrb, extPot, orb)
```

### Caching Pattern

```fortran
if (SysPotential_timeIndependentQ .and. .not. allocated(cachedPot)) then
  call SysPotential_FillExternalPotential(cachedPot, 0.0_R64)
end if

! Use cachedPot on all subsequent time steps
call SysPotential_MultiplyWithExternalPotential(dOrb, cachedPot, orb)
```

---

## JSON Configuration Examples

### 1D Harmonic Trap (Linear)

```json
{
  "sysPotential": {
    "linear": {
      "harmonic": {
        "position": [0.0],
        "omega": [1.0],
        "stdImpl": {}
      }
    }
  }
}
```

### Multi-center Soft-Coulomb (Linear)

```json
{
  "sysPotential": {
    "linear": {
      "softYukawa": {
        "position": [-2.0, 2.0],
        "charge": [1.0, 1.0],
        "softening1": [0.5, 0.5],
        "softening2": [0.0, 0.0],
        "dampening": [0.0, 0.0],
        "stdImpl": {}
      }
    }
  }
}
```

### 3D Harmonic Trap (Lattice)

```json
{
  "sysPotential": {
    "lattice": {
      "harmonic": {
        "omegaX": 0.01,
        "omegaY": 0.01,
        "omegaZ": 0.0,
        "positionX": 2.5,
        "positionY": 2.5,
        "positionZ": 0.5,
        "stdImpl": {}
      }
    }
  }
}
```

### Anderson Disorder (Lattice)

```json
{
  "sysPotential": {
    "lattice": {
      "randomUniform": {
        "minValue": -2.0,
        "maxValue": 2.0,
        "seed": 12345,
        "stdImpl": {}
      }
    }
  }
}
```

### Time-Dependent Tilt (Lattice)

```json
{
  "sysPotential": {
    "lattice": {
      "seesaw": {
        "slopeMinX": -0.1,
        "slopeMaxX": 0.1,
        "frequencyX": 0.5,
        "slopeMinY": 0.0,
        "slopeMaxY": 0.0,
        "frequencyY": 0.0,
        "slopeMinZ": 0.0,
        "slopeMaxZ": 0.0,
        "frequencyZ": 0.0,
        "stdImpl": {}
      }
    }
  }
}
```

### Hydrogen Atom (Ylm)

```json
{
  "sysPotential": {
    "ylm": {
      "coulomb": {
        "charge": 1.0,
        "stdImpl": {}
      }
    }
  }
}
```

### Helium Ion (Ylm)

```json
{
  "sysPotential": {
    "ylm": {
      "coulomb": {
        "charge": 2.0,
        "stdImpl": {}
      }
    }
  }
}
```

---

## Adding a New Potential

### Step 1: Create Directory Structure

```
src/SysPotential/<Grid>/<NewPotential>/
├── M_SysPotential_<Grid>_<NewPotential>.f90
├── S_SysPotential_<Grid>_<NewPotential>.f90
└── StdImpl/
    ├── M_SysPotential_<Grid>_<NewPotential>_StdImpl.f90
    └── S_SysPotential_<Grid>_<NewPotential>_StdImpl.f90
```

### Step 2: Interface Module (`M_*.f90`)

```fortran
module M_SysPotential_Linear_NewPotential
  use M_Utils_Types
  implicit none

  interface
    module subroutine SysPotential_Linear_NewPotential_Fabricate
    end subroutine
  end interface

  ! Model parameters
  real(R64) :: SysPotential_Linear_NewPotential_param1
  ! ...

end module
```

### Step 3: Implementation Submodule (`S_*.f90`)

```fortran
submodule(M_SysPotential_Linear_NewPotential) S_SysPotential_Linear_NewPotential
  implicit none
contains

  module subroutine SysPotential_Linear_NewPotential_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysPotential
    use M_SysPotential_Linear_NewPotential_StdImpl

    call Say_Fabricate("sysPotential.linear.newPotential")

    ! Read parameters
    SysPotential_Linear_NewPotential_param1 = &
      Json_Get("sysPotential.linear.newPotential.param1", 1.0_R64)

    ! Set optimization flags
    SysPotential_timeIndependentQ = .true.
    SysPotential_bodyTypeIndependentQ = .true.

    ! Dispatch to implementation
    if (Json_GetExistence("sysPotential.linear.newPotential.stdImpl")) then
      call SysPotential_Linear_NewPotential_StdImpl_Fabricate
    else
      error stop "sysPotential.linear.newPotential is missing: stdImpl"
    end if

  end subroutine

end submodule
```

### Step 4: Wire Into Parent

In `S_SysPotential_Linear.f90`, add:
```fortran
use M_SysPotential_Linear_NewPotential
! ...
else if (Json_GetExistence("sysPotential.linear.newPotential")) then
  call SysPotential_Linear_NewPotential_Fabricate
```

### Step 5: Add Tests

Create `test/T_SysPotential_NewPotential.f90` and corresponding JSON fixture.

---

## Mathematical Details

### Point-wise Multiplication

For potentials diagonal in position space (most cases):
```
(V̂ψ)(r) = V(r) × ψ(r)
```
This is implemented as element-wise array multiplication.

### Ylm Spatial Product

For spherical harmonics, applying a potential involves Gaunt coefficients:
```
(Vψ)_{l,m} = Σ_{l',m'; l'',m''} G(l',m'; l'',m''; l,m) V_{l',m'}(r) ψ_{l'',m''}(r)
```

This is handled by `Grid_Ylm_SpatialProduct` and is O(lmax⁴ × N_radial).

### Seesaw Time-Dependence

The seesaw potential has slopes that oscillate sinusoidally:
```
slope(t) = (slope_max + slope_min)/2 + (slope_max - slope_min)/2 × sin(2π × freq × t)
```

This creates an effective "tilting" or "rocking" of the lattice potential.

---

## Dependencies

| Module | Usage |
|--------|-------|
| `M_Grid*` | Grid_nPoints, coordinate arrays |
| `M_Utils_Json` | Configuration parsing |
| `M_Utils_Say` | Logging during fabrication |
| `M_Utils_Types` | R64, I32 kinds |
| `M_Utils_Constants` | PI (for Seesaw, Coulomb) |

---

## Performance Considerations

| Grid | Notes |
|------|-------|
| Lattice | All O(N) with minimal overhead |
| Linear | O(N) for evaluation, O(N×M) for M centers |
| Ylm | O(lmax⁴ × N_rad) for spatial product application |

### Caching Recommendations

| Flag Setting | Recommendation |
|--------------|----------------|
| `timeIndependentQ = .true.` | Compute potential once, reuse |
| `bodyTypeIndependentQ = .true.` | Same potential for all particle types |
| Both true | Single cached array for entire simulation |

---

## Common Pitfalls

1. **Forgetting `SysPotential_Setup`:** Some backends may require initialization.

2. **Not checking flags:** Recomputing time-independent potentials every step
   wastes cycles.

3. **Wrong array sizes:** `FillExternalPotential` allocates based on `Grid_nPoints`;
   ensure grid is set up first.

4. **Ylm lmax mismatch:** The potential `lmax` (usually 0 for Coulomb) differs
   from the orbital `lmax`; the spatial product handles the coupling.

5. **Seesaw frequency = 0:** Results in constant slope = (max+min)/2.

6. **Random seed = -1:** Uses system clock, giving different disorder each run.
   Set a fixed seed for reproducibility.

---

## File Summary

| File | Purpose |
|------|---------|
| `M_SysPotential.f90` | Public interface, procedure pointers, flags |
| `S_SysPotential.f90` | Top-level fabrication dispatch |
| `*/M_*.f90` | Grid/model interface modules |
| `*/S_*.f90` | Implementation submodules |
| `*/StdImpl/*` | Default numerical implementations |

Total: ~30 files across all backends.
