# Method Module — Architecture & Onboarding Guide

## Overview

The **Method** module defines the unified dynamics interface for quantum many-body
simulations in CodyFortranRDM. It provides the time-derivative contract
`Method_TimeDerivative(dState, state, t)` that all propagators consume, enabling
the framework to evolve arbitrary quantum states through time.

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                            Propagator / Integrator                          │
│                      calls Method_TimeDerivative(dState, state, t)          │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                              M_Method (interface)                            │
│         Method_state(:)  |  Method_Setup  |  Method_TimeDerivative          │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
              ┌───────────────────────┴───────────────────────┐
              ▼                                               ▼
       ┌─────────────┐                               ┌───────────────┐
       │  Sb (single │                               │ Mb (many-body)│
       │   particle) │                               │               │
       └─────────────┘                               └───────────────┘
                                                              │
                          ┌───────────────────────────────────┼───────────────┐
                          ▼                                   ▼               ▼
                   ┌─────────────┐                    ┌─────────────┐  ┌─────────────┐
                   │  OrbBased   │                    │  GridBased  │  │  GemBased   │
                   │ (orbital    │                    │ (direct     │  │ (Gaussian   │
                   │  expansion) │                    │  product)   │  │  expansion) │
                   └─────────────┘                    └─────────────┘  └─────────────┘
                          │
       ┌──────────┬───────┴───────┬────────────┬─────────────────┐
       ▼          ▼               ▼            ▼                 ▼
    ┌──────┐  ┌──────┐       ┌────────┐   ┌────────┐    ┌─────────────────┐
    │ TDCI │  │ TDHX │       │ MCTDHX │   │TD-2RDM │    │TD-2RDM-Static   │
    │      │  │      │       │        │   │        │    │     Orbs        │
    └──────┘  └──────┘       └────────┘   └────────┘    └─────────────────┘
```

---

## Key Concepts

### 1. Packed State Vector

All methods operate on a **single contiguous 1D complex array** `Method_state(:)`.
The internal layout is method-specific:

| Method              | State Layout                                      |
|---------------------|---------------------------------------------------|
| Sb                  | Single wavefunction ψ(r) on the grid              |
| TDCI                | CI coefficients c_I only (static orbitals)        |
| TDHX                | Orbitals φ_i(r) only (single determinant)         |
| MCTDHX              | Orbitals φ_i(r) ⊕ CI coefficients c_I             |
| TD-2RDM             | Orbitals φ_i(r) ⊕ 2-RDM D²_{ijkl}                 |
| TD-2RDM-StaticOrbs  | 2-RDM D²_{ijkl} only (fixed orbitals)             |
| GridBased/Full      | Full many-body wavefunction Ψ(r₁,r₂,...) on grid  |

### 2. Procedure Pointer Binding (Fabrication Pattern)

At runtime, the JSON configuration determines which implementation is bound:

```fortran
! In Method_Fabricate:
if (Json_GetExistence("method.sb")) then
  call Method_Sb_Fabricate     ! binds Sb procedures
else if (Json_GetExistence("method.mb")) then
  call Method_Mb_Fabricate     ! descends into many-body branch
end if
```

Each `*_Fabricate` routine:
1. Reads method-specific JSON parameters
2. Binds the module's procedure pointers (`Method_Setup`, `Method_TimeDerivative`, etc.)
3. Calls the next level's fabricate if further branching is needed

### 3. Time Derivative Contract

The central interface every method must implement:

```fortran
subroutine Method_TimeDerivative(dState, state, time)
  complex(R64), intent(out) :: dState(:)   ! output: -i·H·|Ψ⟩
  complex(R64), intent(in)  :: state(:)    ! input:  |Ψ⟩
  real(R64), intent(in)     :: time        ! current time
end subroutine
```

**Semantics:**
- Computes dΨ/dt = −i·Ĥ(t)·Ψ (ħ=1 units)
- Overwrites `dState`; no aliasing with `state`
- May use cached operators or module-level data

---

## Module Hierarchy

### Top Level: `M_Method` / `S_Method`

| Export                    | Purpose                                          |
|---------------------------|--------------------------------------------------|
| `Method_state(:)`         | The packed quantum state (allocated by Setup)    |
| `Method_Setup`            | Procedure pointer → initialize state & caches    |
| `Method_TimeDerivative`   | Procedure pointer → compute −i·H·|Ψ⟩             |
| `Method_GetEnergy`        | Procedure pointer → compute ⟨Ψ|H|Ψ⟩              |
| `Method_Fabricate`        | Entry point; branches to Sb or Mb                |

### Single-Body Branch: `M_Method_Sb`

Solves the single-particle TDSE under an external potential V(r,t):
```
i ∂ψ/∂t = (T̂ + V̂) ψ
```
No mean-field or exchange terms. Useful for testing grids and basic physics.

### Many-Body Branch: `M_Method_Mb`

Shared metadata for all many-body methods:

| Data                           | Purpose                                    |
|--------------------------------|--------------------------------------------|
| `Method_Mb_nBodyTypes`         | Number of particle species/spin types      |
| `Method_Mb_nBodies(:)`         | Particle count per body type               |
| `Method_Mb_bodyStatistics(:)`  | 'f' (fermion) or 'b' (boson) per type      |
| `Method_Mb_nBodiesStart/End(:)`| Index ranges for packed body-type storage  |

Branches to: **OrbBased**, **GridBased**, or **GemBased**.

---

## OrbBased Methods (Most Feature-Rich)

### Shared Infrastructure: `M_Method_Mb_OrbBased`

Provides orbital indexing, RDM storage, and operator-application interfaces:

| Procedure Pointer                              | Signature / Purpose                        |
|------------------------------------------------|--------------------------------------------|
| `Method_Mb_OrbBased_FillH1`                    | Build 1-body Hamiltonian h¹_{ij}           |
| `Method_Mb_OrbBased_FillH2`                    | Build 2-body Hamiltonian h²_{ijkl}         |
| `Method_Mb_OrbBased_FillRdm1`                  | Extract 1-RDM from state                   |
| `Method_Mb_OrbBased_FillRdm2`                  | Extract 2-RDM from state                   |
| `Method_Mb_OrbBased_ApplyKineticOp`            | T̂ · orbitals (accumulate)                  |
| `Method_Mb_OrbBased_ApplyPotentialOp`          | V̂_ext · orbitals (accumulate)              |
| `Method_Mb_OrbBased_ApplyInteractionOp`        | Ŵ · orbitals for index pair (i₂,j₂)        |
| `Method_Mb_OrbBased_ApplyHartreeFockOp`        | Mean-field + exchange (Hartree–Fock)       |
| `Method_Mb_OrbBased_ApplyCorrelationOp`        | Beyond-mean-field correlation contribution |
| `Method_Mb_OrbBased_ApplySingleBodyOp`         | Combined T̂ + V̂_ext                         |
| `Method_Mb_OrbBased_TimeDerivativeOrbsLin`     | Linear (static-orbital) time derivative    |
| `...TimeDerivativeCoeffsPlusOrbsNonLin`        | Non-linear part for MCTDHX-like schemes    |

### Concrete Methods

| Method           | Orbitals   | Coefficients | Description                                |
|------------------|------------|--------------|-------------------------------------------|
| **TDCI**         | Static     | Dynamic      | Full CI; exact within orbital basis        |
| **TDHX**         | Dynamic    | Implicit     | Single determinant; mean-field + exchange  |
| **MCTDHX**       | Dynamic    | Dynamic      | Multiconfiguration; Dirac–Frenkel VP       |
| **TD-2RDM**      | Dynamic    | 2-RDM        | Propagates 2-RDM with adaptive orbitals    |
| **TD-2RDM-Static**| Static    | 2-RDM        | Propagates 2-RDM only; fixed basis         |

---

## GridBased Methods

Represents the full many-body wavefunction Ψ(r₁,r₂,...,rₙ) directly on a tensor
product grid. Exponentially expensive but numerically exact for small N.

| Procedure Pointer                              | Purpose                            |
|------------------------------------------------|------------------------------------|
| `Method_Mb_GridBased_ApplyKineticOp`           | Laplacian on 1D grid slice         |
| `Method_Mb_GridBased_ApplyPotentialOp`         | V(r) on 1D grid slice              |
| `Method_Mb_GridBased_ApplyInteractionOp`       | W(r₁,r₂) on 2D grid slice          |

Currently supports: **Full** (complete tensor-product representation).

---

## GemBased Methods (Placeholder)

Intended for Gaussian-Expansion Methods (GEM) where the wavefunction is expanded
in explicitly correlated Gaussian basis functions. Infrastructure is in place
but concrete implementations are not yet finalized.

---

## JSON Configuration Examples

### Single-Body (Sb)
```json
{
  "grid": { "linear": { "const": { "nPoints": 128, "xMin": -10, "xMax": 10 } } },
  "sysPotential": { "linear": { "harmonic": { "omega": 1.0 } } },
  "method": { "sb": {} }
}
```

### Many-Body MCTDHX (2 fermions, 4 orbitals)
```json
{
  "method": {
    "mb": {
      "nBodyTypes": 1,
      "nBodies": [2],
      "bodyStatistics": ["f"],
      "orbBased": {
        "nOrbs": [4],
        "mctdhx": {}
      }
    }
  }
}
```

### Many-Body TD-2RDM with Static Orbitals
```json
{
  "method": {
    "mb": {
      "nBodyTypes": 2,
      "nBodies": [2, 2],
      "bodyStatistics": ["f", "f"],
      "orbBased": {
        "nOrbs": [4, 4],
        "td2rdmStaticOrbs": {}
      }
    }
  }
}
```

---

## Adding a New Method

1. **Create directory structure:**
   ```
   src/Method/Mb/OrbBased/MyMethod/
       M_Method_Mb_OrbBased_MyMethod.f90   (interface module)
       S_Method_Mb_OrbBased_MyMethod.f90   (submodule implementation)
   ```

2. **Define the interface module** with a `*_Fabricate` subroutine.

3. **Implement the submodule:**
   - Bind `Method_Setup` → allocate `Method_state`, initialize
   - Bind `Method_TimeDerivative` → compute −i·H·|Ψ⟩
   - Optionally bind `Method_GetEnergy`

4. **Register in parent fabricate:**
   Add a branch in `S_Method_Mb_OrbBased.f90`:
   ```fortran
   else if (Json_GetExistence("method.mb.orbBased.myMethod")) then
     call Method_Mb_OrbBased_MyMethod_Fabricate
   ```

5. **Update CMakeLists.txt** to include the new source files.

---

## Testing

Tests live in `test/` and follow the naming convention `T_Method_*.f90` with
corresponding `T_Method_*.json` configuration files.

Run tests:
```bash
cd build && ctest -R T_Method
```

---

## Common Pitfalls

1. **Forgetting to initialize `dState` to zero** — `TimeDerivative` overwrites,
   but operator-application routines *accumulate*. Ensure proper initialization.

2. **Body-type index mismatch** — When looping over orbitals, always check
   `bodyTypeOfOrb(i) == bodyTypeOfOrb(j)` before computing matrix elements
   between body-type-restricted quantities.

3. **Restricted vs. unrestricted orbitals** — In restricted calculations,
   `Orbs_nOrbsInState < Method_Mb_OrbBased_nOrbsSum`. The code duplicates
   blocks in h1/h2 matrices; be careful with index arithmetic.

4. **Regularization** — `Method_Mb_OrbBased_regularizationParameter` prevents
   division by zero in natural-orbital transformations. Keep it small (1e-10)
   but not zero.

---

## Key Equations

### TDSE (all methods)
```
i ∂|Ψ⟩/∂t = Ĥ(t) |Ψ⟩
```

### Energy Expectation
```
E(t) = ⟨Ψ(t)| Ĥ(t) |Ψ(t)⟩ = Tr[ρ¹ h¹] + ½ Tr[ρ² h²]
```

### MCTDHX Equations of Motion
```
i ∂c_I/∂t = Σ_J ⟨Φ_I| Ĥ |Φ_J⟩ c_J

i ∂|φ_j⟩/∂t = (1 - P̂) Σ_k (ρ¹)⁻¹_{jk} ĥ^{mf}_k |φ_k⟩
```
where P̂ is the projector onto the occupied orbital space.

### TD-2RDM (BBGKY Hierarchy)
```
i ∂D²/∂t = [h¹, D²] + [Ŵ, D³]   (D³ approximated via reconstruction)
```

---

## File Summary

| File                                  | Role                                      |
|---------------------------------------|-------------------------------------------|
| `M_Method.f90`                        | Top-level interface; state & pointers     |
| `S_Method.f90`                        | Fabricate dispatcher (Sb vs Mb)           |
| `Sb/M_Method_Sb.f90`                  | Single-body interface                     |
| `Sb/S_Method_Sb.f90`                  | Single-body implementation                |
| `Mb/M_Method_Mb.f90`                  | Many-body shared data (body types, etc.)  |
| `Mb/S_Method_Mb.f90`                  | Many-body fabricate dispatcher            |
| `Mb/OrbBased/M_Method_Mb_OrbBased.f90`| Orbital-based shared data & interfaces    |
| `Mb/OrbBased/S_Method_Mb_OrbBased.f90`| Orbital-based operator implementations    |
| `Mb/OrbBased/Tdci/*`                  | Full CI with static orbitals              |
| `Mb/OrbBased/Tdhx/*`                  | Time-dependent Hartree–Fock/Exchange      |
| `Mb/OrbBased/Mctdhx/*`                | Multiconfiguration TDHX                   |
| `Mb/OrbBased/Td2rdm/*`                | TD-2RDM with dynamic orbitals             |
| `Mb/OrbBased/Td2rdmStaticOrbs/*`      | TD-2RDM with static orbitals              |
| `Mb/GridBased/M_Method_Mb_GridBased.f90` | Grid-based interfaces                  |
| `Mb/GridBased/S_Method_Mb_GridBased.f90` | Grid-based operator implementations    |
| `Mb/GridBased/Full/*`                 | Full tensor-product grid method           |
| `Mb/GemBased/*`                       | Gaussian expansion (placeholder)          |

---

## Contact & Further Reading

- Main documentation: `README.md` in repository root
- Build instructions: `CMakeLists.txt`, `setBuildVars.sh`
- Related modules: `M_Propagator`, `M_Orbs`, `M_Coeffs`, `M_TwoRdm`
