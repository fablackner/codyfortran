# SysPotential Module — Agent Onboarding Guide

## Overview

The **SysPotential** module implements external (one-body) potentials for the
CodyFortran quantum many-body simulation framework. It provides the operator
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
├── Ylm/                            ← Spherical harmonics backends
│   ├── M_SysPotential_Ylm.f90
│   ├── S_SysPotential_Ylm.f90
│   └── Coulomb/                    ← Central Coulomb -Z/r
│
└── Fem/                            ← 3D spectral elements (grid.fem)
    ├── M_SysPotential_Fem.f90      ← family interface (Fabricate), no data
    ├── S_SysPotential_Fem.f90      ← point-wise Multiply, dispatch
    └── Coulomb/                    ← nuclei of the grid.fem centers (Galerkin or erf-smeared):
                                       charges, smearing, nuclear repulsion; private static cache
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
| `SysPotential_Ylm_mmax` | Maximum \|m\| of the potential expansion; must be set by the model (Coulomb/CoulombLaser: 0); Fabricate refuses models with mmax > 0 when `grid.ylm.mmax < grid.ylm.lmax` |
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

### 4. Fem (3D Cartesian Spectral Elements)

**Use case:** molecules (all-electron, several nuclei) on `grid.fem`.

Potentials are nodal arrays (`Grid_nPoints`); the family binds the point-wise
`MultiplyWithExternalPotential` (dOrb = V·orb; reentrant: serial inside the
OpenMP loops of Method, parallel over the nodes when called outside a
parallel region, e.g. by the Cayley exact part of Method), the model binds
`Setup` and `FillExternalPotential`. Family modules bind through
`use M_SysPotential, only: SysPotential_Setup, ...` and never redeclare
`SysPotential_Setup` (a redeclared pointer would be a dead copy).

**Implementations:**

| Type | Formula | Time-dep? | Parameters |
|------|---------|-----------|------------|
| Coulomb | V(r) = Σ_I V_I(r) over the `grid.fem` centers | No | charges[], smearing[] |

with, per center I (position `Grid_Fem_centers(:, I)`, charge Z_I, radius rc_I):

- **rc_I = 0 — Galerkin point nucleus.** All point nuclei together are one
  `Grid_Fem_SolvePoissonPointCharges` solve with q_I = −Z_I: the load is
  q_I φ_i(R_I) (q_I δ_iJ for a nucleus on node J) and the box values are the
  exact Σ q_I/|r_B − R_I|. Away from the nucleus V = −Z/r up to the mesh error;
  the value **at the nucleus node is finite**, roughly −30 … −50 Z/hMin for
  nLocals 5–6, depending on the mesh (measured: −27.6 Z/hMin for nLocals 5,
  growth 1.8; −40.9 … −49.6 Z/hMin for nLocals 6, growth 1.8 … 1.6; about
  −85 Z/hMin for nLocals 8), hMin the local element size.
  The attraction Σ_i w_i ρ_i V_i equals −Z × (Galerkin Hartree potential of
  ρ)(R_I), i.e. it is exactly as accurate as the FEM Hartree potential (the
  hydrogen ⟨1s|V|1s⟩ is −0.99999874 at hMin = 0.2). Galerkin nuclei are not
  variational: total energies may lie slightly below the basis-set limit.
- **rc_I > 0 — smeared nucleus.** Nodal values −Z erf(|r − R|/rc)/|r − R|
  with the limit −2Z/(√π rc) on the nucleus node (no division by r there;
  `SysPotential_Fem_Coulomb_SmearedPotential(charge, smearing, r)`).
- **Z_I = 0 — ghost center** (mesh refinement only).

**Why smearing helps explicit time steps.** The spike of a Galerkin point
nucleus (roughly −30 … −50 Z/hMin: about −1600 Ha for O at hMin = 0.2, and
−6500 Ha (−40.9 Z/hMin, nLocals 6) at the hMin = 0.05 that resolves its cusp
on the accurate H₂O mesh of app/fem_molecules) and the fine core elements
stretch the spectrum of h = T + V, so explicit integrators of the full h need
dt ≲ 2.8/‖h‖, about 4e-4 a.u. for O. A smeared nucleus is bounded by
2Z/(√π rc) and has no cusp to resolve, so the core mesh can be coarse and
explicit steps of about 0.01–0.05 a.u. become possible; only the core
physics changes. With the exact field-free part
(`method.mb.orbBased.exactFieldFreeQ`, h0 = T + V_static propagated
implicitly by a Crank–Nicolson/Cayley step on Fem, A_Method.md "Exact
Field-Free Part on grid.fem") the spike does not limit the step, and point
nuclei are the accurate choice (T_Fem_05–07: dt = 0.01–0.02 with the
Galerkin point nuclei of H2).

**Static by design.** Method refills the external potential on every
right-hand-side evaluation, so Setup builds the static potential once (a
private cache of `S_SysPotential_Fem_Coulomb`) and `FillExternalPotential`
only copies it (any time, any `bt_`; read it through
`SysPotential_FillExternalPotential`). The laser is not a SysPotential model on
Fem: the length-gauge dipole coupling is `sysGauge.fem.length`, which keeps
SysPotential static so that h0 = T + V_static can be treated exactly.

**Exported data (M_SysPotential_Fem_Coulomb):** the model owns its parameters,
as `SysPotential_Ylm_Coulomb_charge` does; the family module
`M_SysPotential_Fem` holds no data, so another Fem model never leaves stale
Coulomb data behind. `SysPotential_Fem_Coulomb_charges(nCenters)`,
`SysPotential_Fem_Coulomb_smearing(nCenters)`,
`SysPotential_Fem_Coulomb_nuclearRepulsion` (Σ_{I<J} Z_I Z_J/R_IJ of the point
charges, set at fabrication, printed at Setup; not included in Method
energies). Used by app/fem_molecules (`P_FemScf`: nuclear dipole and E_tot;
`P_FemKickTdhf`: E_tot) and T_Fem_03/04.

**JSON path:** `sysPotential.fem.coulomb` (no implementation leaf)

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `charges` | real(:) | required | Z_I, one per `grid.fem` center (size = nCenters), 0 = ghost |
| `smearing` | real(:) | 0 each | rc_I ≥ 0, one per center: 0 = Galerkin point nucleus, > 0 = erf smearing |

Errors stop with clear messages for a missing `charges`, sizes ≠ nCenters,
a negative radius and two charged centers on the same position.

**Test:** `test/componentTests/SysPotential/Fem/Coulomb/T_SysPotential_Fem_Coulomb`
(box ±8, point proton at the origin with hMin 0.2, smeared Z = 2 and Z = 3
nuclei, one ghost; N = 485139): static = Galerkin + erf part to 7e-17
(relative), r = 0 limits exact; point proton max|V + 1/r| = 3.9e-6 (r > 1.5),
5.1e-7 (r > 3), −247.94 = −49.6/hMin at the nucleus node, ⟨1s|V|1s⟩ =
−0.99999874; nuclear repulsion exact to 1e-14.

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

### Water (Fem, grid.fem centers O, H, H)

```json
{
  "sysPotential": {
    "fem": {
      "coulomb": {
        "charges": [8.0, 1.0, 1.0],
        "smearing": [0.0, 0.0, 0.0]
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
| `M_Grid_Fem` | Centers, nodes, `Grid_Fem_SolvePoissonPointCharges` (Fem) |

---

## Performance Considerations

| Grid | Notes |
|------|-------|
| Lattice | All O(N) with minimal overhead |
| Linear | O(N) for evaluation, O(N×M) for M centers |
| Ylm | O(lmax⁴ × N_rad) for spatial product application |
| Fem | Setup: one FDM point-charge solve (~25 ms at N = 4e5) + O(N) per smeared nucleus; Fill: one copy; Multiply: O(N) |

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

7. **Redeclared `SysPotential_Setup`:** the Linear, Lattice and Ylm family
   modules declare their own pointer of that name; a model that binds it
   through such a module sets a dead copy. Bind through
   `use M_SysPotential, only: SysPotential_Setup` (as Fem does).

8. **Nuclei on nodes (Fem):** every nucleus sits on a grid node, so never
   evaluate 1/r (or z/r³) there; Galerkin point nuclei and the erf limit are
   finite by construction.

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
