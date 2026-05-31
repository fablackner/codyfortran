# SysInteraction Module — Agent Onboarding Guide

## Overview

The **SysInteraction** module implements particle-particle (two-body) interactions
for the CodyFortranRDM quantum many-body simulation framework. It provides the
operator **Ŵ** that appears in the many-body Hamiltonian:

```
Ĥ = T̂ + V̂_ext + Ŵ
```

where `Ŵ = ½ Σᵢⱼ w(rᵢ,rⱼ)` describes electron-electron repulsion, Hubbard-U
interactions, or other pairwise potentials.

---

## Architecture

### Module/Submodule Pattern

This module follows the framework's **Interface Module + Implementation Submodule**
paradigm:

```
M_SysInteraction.f90          ← Public interface (abstract contracts)
S_SysInteraction.f90          ← Fabrication logic (runtime dispatch)
├── Lattice/
│   ├── M_SysInteraction_Lattice.f90
│   ├── S_SysInteraction_Lattice.f90
│   └── OnSite/...
├── Linear/
│   ├── M_SysInteraction_Linear.f90
│   ├── S_SysInteraction_Linear.f90
│   └── SoftYukawa/...
└── Ylm/
    ├── M_SysInteraction_Ylm.f90
    ├── S_SysInteraction_Ylm.f90
    └── Coulomb/...
```

**Key rule:** Interface modules (`M_*.f90`) only import other `M_*` modules.
Implementation submodules (`S_*.f90`) contain the actual algorithms.

### Procedure Pointer Binding

At runtime, `SysInteraction_Fabricate` reads the JSON configuration and binds
procedure pointers to concrete implementations. This enables:

- Swapping algorithms without recompilation
- JSON-driven physics configuration
- Clean separation of interface from implementation

---

## Core API

### Public Procedure Pointers

| Pointer | Signature | Purpose |
|---------|-----------|---------|
| `SysInteraction_Setup` | `()` | One-time initialization (precompute kernels, FFT plans) |
| `SysInteraction_FillInteractionSrc` | `(src, orbConjg, orb)` | Build source ρ = ψ*·ψ |
| `SysInteraction_FillInteractionPotential` | `(pot, src, t, bt1, bt2)` | Solve for V from ρ |
| `SysInteraction_MultiplyWithInteractionPotential` | `(dOrb, pot, orb)` | Apply V·ψ |

### Public Data

| Variable | Type | Meaning |
|----------|------|---------|
| `SysInteraction_timeIndependentQ` | `logical` | True if Ŵ is time-independent (enables caching) |
| `SysInteraction_bodyTypeIndependentQ` | `logical` | True if Ŵ treats all species identically |

---

## Grid Back-Ends

### 1. Lattice (Discrete)

**Use case:** Hubbard models, tight-binding systems

**Physics:** On-site density-density interaction
```
Ŵ = U Σᵢ n↑ᵢ n↓ᵢ
```

**Implementations:**
- `OnSite/StdImpl`: Direct V(i) = U·ρ(i)

**JSON path:** `sysInteraction.lattice.onSite`

### 2. Linear (1D Real-Space)

**Use case:** 1D model systems, trapped gases

**Physics:** Convolution with interaction kernel
```
V(x) = ∫ w(|x-x'|) ρ(x') dx'
```

**Implementations:**
- `SoftYukawa/StdImpl`: O(N²) direct convolution
- `SoftYukawa/Fftw3`: O(N log N) FFT convolution

**Kernel form:**
```
w(r) = strength × exp(-dampening×r) / (√(r² + softening1²) + softening2)
```

**JSON path:** `sysInteraction.linear.softYukawa`

### 3. Ylm (Spherical Harmonics)

**Use case:** Atoms, molecules with spherical symmetry

**Physics:** Coulomb 1/r expanded in multipoles
```
1/|r₁-r₂| = Σₗₘ (4π/(2l+1)) × (r<ˡ/r>ˡ⁺¹) × Yₗₘ*(Ω₁)Yₗₘ(Ω₂)
```

**Implementations:**
| Variant | Complexity | Notes |
|---------|------------|-------|
| `StdImpl` | O(N²) | Reference, Green's function |
| `TwoScan` | O(N) | Forward/backward prefix sums |
| `FullEq` | O(N³) | FEDVR + LU per call |
| `BlockEq` | O(nE² + nE×nLoc) | FEDVR + precomputed Schur |

**JSON path:** `sysInteraction.ylm.coulomb`

---

## Workflow

### Initialization Sequence

```fortran
call SysInteraction_Fabricate    ! Parse JSON, bind pointers
call SysInteraction_Setup        ! Precompute kernels/plans
```

### Typical Usage (in Method/Hamiltonian)

```fortran
! Build source from orbital product
call SysInteraction_FillInteractionSrc(src, conjg(orb_i), orb_j)

! Solve for interaction potential
call SysInteraction_FillInteractionPotential(pot, src, time)

! Apply to orbital
call SysInteraction_MultiplyWithInteractionPotential(dOrb, pot, orb_k)
```

---

## JSON Configuration Examples

### Hubbard Model (Lattice)
```json
{
  "sysInteraction": {
    "lattice": {
      "onSite": {
        "strength": 4.0,
        "stdImpl": {}
      }
    }
  }
}
```

### 1D Soft-Coulomb (Linear)
```json
{
  "sysInteraction": {
    "linear": {
      "softYukawa": {
        "strength": 1.0,
        "softening1": 0.5,
        "fftw": {}
      }
    }
  }
}
```

### Helium Atom (Ylm)
```json
{
  "sysInteraction": {
    "ylm": {
      "lmax": 6,
      "coulomb": {
        "strength": 1.0,
        "twoScan": {}
      }
    }
  }
}
```

---

## Adding a New Interaction

1. **Create directory structure:**
   ```
   src/SysInteraction/<Grid>/<NewModel>/
   ├── M_SysInteraction_<Grid>_<NewModel>.f90
   └── S_SysInteraction_<Grid>_<NewModel>.f90
   ```

2. **In the interface module (`M_*.f90`):**
   - Declare any model-specific parameters
   - Declare the `_Fabricate` subroutine interface

3. **In the implementation submodule (`S_*.f90`):**
   - Implement `_Fabricate` to read JSON and bind pointers
   - Implement the actual physics routines

4. **Wire into parent:**
   - Add `use M_SysInteraction_<Grid>_<NewModel>` in parent fabricate
   - Add JSON existence check and dispatch call

5. **Add tests:**
   - Create `test/T_SysInteraction_<NewModel>.f90`
   - Create corresponding JSON fixture

---

## Performance Considerations

| Grid | Recommendation |
|------|----------------|
| Lattice | Always fast (local interaction) |
| Linear | Use `fftw` for N > 500 |
| Ylm | Use `twoScan` for large grids; `blockEq` if many repeated solves |

### Caching Opportunities

When `timeIndependentQ = .true.` and `bodyTypeIndependentQ = .true.`:
- The interaction potential can be computed once and reused
- Higher-level code (Method) should exploit this

---

## Dependencies

| Module | Usage |
|--------|-------|
| `M_Grid*` | Grid points, weights, nPoints |
| `M_Utils_Json` | Configuration parsing |
| `M_Utils_Fftw3Lib` | FFT for convolution |
| `M_Utils_LapackLib` | LU factorization for FEDVR solvers |
| `M_Utils_Fedvr` | Finite element DVR infrastructure |
| `M_Utils_ConvolutionIntegral` | Direct convolution utility |
| `M_Utils_ConvolutionFftw` | FFT convolution utility |

---

## Common Pitfalls

1. **Forgetting `SysInteraction_Setup`:** Some implementations (FFT, BlockEq)
   require setup to precompute plans/factorizations.

2. **Wrong `lmax` for Ylm:** The potential expansion needs `lmax_pot ≥ 2×lmax_orb`
   to capture the full density product.

3. **Missing quadrature weights:** Linear/Ylm source terms include weights;
   Lattice does not (discrete sum).

4. **Body type arguments:** Pass `bt1_`, `bt2_` if simulating multiple species
   with different interactions.

---

## File Summary

| File | Lines | Purpose |
|------|-------|---------|
| `M_SysInteraction.f90` | ~120 | Public interface, procedure pointers |
| `S_SysInteraction.f90` | ~60 | Top-level fabrication dispatch |
| `*/M_*.f90` | ~30-50 | Grid/model interface modules |
| `*/S_*.f90` | ~50-450 | Implementation submodules |

Total: ~25 files, ~2000 lines of physics code.
