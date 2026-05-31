# ConfigList Module — Agent Onboarding Guide

## Overview

The **ConfigList** module provides Fock-basis enumeration and second-quantized excitation machinery for quantum many-body simulations. It is the foundation for configuration-interaction (CI), MCTDH, and TD-2RDM methods.

**In plain terms:** Given N particles distributed among M orbitals, this module:
1. Enumerates all valid configurations (Fock states)
2. Precomputes which configurations connect via 1- and 2-body excitations
3. Applies creation/annihilation operators with proper statistics (fermionic signs, bosonic ladder factors)

---

## Architecture

```
M_ConfigList                          ← Abstract interface + global container
    │
    └── S_ConfigList                  ← Factory: JSON → concrete elements
            │
            └── M_ConfigList_AllActive    ← "All-active" truncation base
                    │
                    ├── Fermions/
                    │   ├── M_ConfigList_AllActive_Fermionic  ← Type definition
                    │   └── S_ConfigList_AllActive_Fermionic  ← Bit-pattern encoding, ±1 signs
                    │
                    └── Bosons/
                        ├── M_ConfigList_AllActive_Bosonic    ← Type definition
                        └── S_ConfigList_AllActive_Bosonic    ← Base-(N+1) encoding, √n factors
```

**Design Pattern:** Interface Module (`M_*.f90`) + Implementation Submodule (`S_*.f90`)

---

## Key Concepts

### Configuration Encoding

| Statistics | Encoding | Example (3 particles, 5 orbitals) |
|------------|----------|-----------------------------------|
| Fermionic  | Bit-pattern: bit i = 1 ⟺ orbital i occupied | `|1,3,5⟩` → `0b10101` = 21 |
| Bosonic    | Base-(N+1): digit i = occupation of orbital i | `|1,1,3⟩` → `2·4⁰ + 0·4¹ + 1·4²` / 3 |

### Excitation Connectivity

During `Setup`, the module precomputes **sparse connectivity graphs**:

```
singles.excitedC(k, iC) = target config from iC via k-th single excitation
singles.factor(k, iC)   = scalar factor (±1 for fermions, √n for bosons)
singles.orbCode(k, iC)  = packed orbital indices (i,j) for a†_i a_j
```

Same structure for `doubles` with 4 orbital indices.

This enables O(1) Hamiltonian matrix construction in the CI/RDM modules.

---

## Data Flow

```
JSON Input                    Runtime
─────────────────────────────────────────────────────────
{                             
  "configList": {             ConfigList_Fabricate()
    "allActive1": {               │
      "fermionic": {              ├─→ Allocate concrete type
        "bodyTarget": 1,          ├─→ Read params, compute nConfigurations
        "nExcitations": 2         │
      }                       ConfigList_Setup()
    }                             │
  }                               ├─→ Build codeFromConfig mapping
}                                 ├─→ Enumerate all singles/doubles
                                  └─→ Store connectivity in sparse arrays
                              
                              ExciteConfiguration(iCNew, factor, [i], [j], iC)
                                  │
                                  └─→ Apply a†_i a_j |iC⟩, return |iCNew⟩ and factor
```

---

## Key Procedures

| Procedure | Purpose |
|-----------|---------|
| `ConfigList_Fabricate` | Parse JSON, allocate elements, call `Fabricate()` |
| `ConfigList_Setup` | Iterate all elements, call `Setup()` |
| `ExciteConfiguration(iCNew, factor, creates, destroys, iC)` | Apply operators to config |

### ExciteConfiguration Interface

```fortran
call element%ExciteConfiguration(iCNew, factor, [i1, i2], [j1, j2], iC)
!                                  │      │       │         │       └─ input config
!                                  │      │       │         └─ annihilate j1, j2
!                                  │      │       └─ create i1, i2
!                                  │      └─ scalar factor (sign or √n)
!                                  └─ output config (0 if forbidden)
```

**Operator ordering:** Right-to-left (physics convention). `a†_1 a†_2` applies `a†_2` first.

---

## Extension Guide

### Adding a New Truncation Scheme (e.g., "restricted")

1. Create directory: `src/ConfigList/Restricted/`
2. Create interface module `M_ConfigList_Restricted.f90`:
   - Define `T_ConfigList_E_Restricted` extending `T_ConfigList_E`
   - Declare `ConfigList_Restricted_Allocate`
3. Create submodule `S_ConfigList_Restricted.f90`:
   - Implement `Setup`, `Fabricate`, `ExciteConfiguration`
4. Register in `S_ConfigList.f90`:
   ```fortran
   else if (index(childName, "restricted") .ne. 0) then
     call ConfigList_Restricted_Allocate(configList(i)%e, "configList."//childName)
   ```

### Adding New Statistics (e.g., anyons)

Follow the pattern in `AllActive/Fermions/` or `AllActive/Bosons/`:
1. Create `M_ConfigList_AllActive_Anyonic.f90` with type extending `T_ConfigList_E_AllActive`
2. Implement statistics-specific `ExciteConfiguration` with fractional phase
3. Register in `S_ConfigList_AllActive.f90` allocator dispatch

---

## Performance Notes

### OpenMP Parallelization

Connectivity graph construction is parallelized over configurations:
```fortran
!$omp parallel do private(iC, ...)
do iC = 1, nConfigurations
  ! Each thread processes independent columns
end do
```

### Memory Layout

Arrays are column-major (Fortran default). `excitedC(k, iC)` stores connections for config `iC` contiguously in the first index, enabling cache-friendly iteration when processing one config at a time.

### Complexity

| Operation | Complexity |
|-----------|------------|
| Setup (singles) | O(nConfigs × nOrbs²) |
| Setup (doubles) | O(nConfigs × nOrbs⁴) |
| ExciteConfiguration | O(nOrbs) per call |
| Hamiltonian construction | O(nConfigs × maxConnections) |

---

## Common Issues

### "index calc failed" during Setup

The round-trip verification `indexOfCombi...(i1(:,iC)) == iC` failed. This indicates a bug in the combinatorial utilities or incorrect parameters (nOrbs, nBodies).

### "bodyTarget not fermionic/bosonic"

The `bodyStatistics` in the Method module doesn't match the ConfigList element type. Check JSON consistency.

### Zero factor returned unexpectedly

For fermions: Pauli exclusion (creating on occupied orbital).
For bosons: Annihilating an empty orbital.

---

## Dependencies

```
M_ConfigList
    ├── M_Utils_Types          (I32, I64, R64)
    ├── M_Utils_NoOpProcedures (NoOpProcedures_Setup)
    ├── M_Utils_Json           (Json_Get, Json_GetExistence, ...)
    ├── M_Utils_Say            (Say_Fabricate, Say_Setup)
    ├── M_Utils_Combinatorics  (indexOfCombiNoRepeat, indexOfCombiWithRepeat)
    ├── M_Utils_CombinationGslLib  (CombiNoRepeat)
    ├── M_Utils_MultisetGslLib     (CombiWithRepeat)
    ├── M_Utils_SfGslLib           (Binomial)
    ├── M_Method_Mb            (nBodies, bodyStatistics)
    └── M_Method_Mb_OrbBased   (nOrbs, nOrbsSum)
```

---

## Testing

Tests are located in `test/` with corresponding JSON input files:
```bash
cd build && ctest -R T_ConfigList
```

Typical test pattern:
1. Fabricate from test JSON
2. Setup
3. Verify nConfigurations matches analytical formula
4. Test specific excitations and verify factors

---

## References

- Helgaker, Jørgensen, Olsen: *Molecular Electronic-Structure Theory*, Ch. 11 (CI methods)
- Second quantization formalism: Fetter & Walecka, *Quantum Theory of Many-Particle Systems*
