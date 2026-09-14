# Purification Module — Developer Onboarding Guide

## Overview

The **Purification** module enforces N-representability conditions on the two-particle reduced density matrix (2-RDM) during time-dependent 2-RDM (TD-2RDM) propagation. During dynamics, numerical errors and approximate BBGKY closures can cause the 2-RDM to acquire unphysical negative eigenvalues. This module provides iterative algorithms to restore positivity.

### Physical Motivation

The 2-RDM must satisfy ensemble N-representability constraints. The two most commonly enforced are:

| Condition | Meaning                                  | Constraint         |
|-----------|------------------------------------------|--------------------|
| **D-condition** | Particle-particle RDM is positive | D²(ijkl) ≥ 0       |
| **Q-condition** | Hole-hole RDM is positive         | Q²(ijkl) ≥ 0       |

When either D² or Q² develops negative eigenvalues, the state is no longer N-representable. The purification algorithm iteratively projects out the defective (negative) eigenspace from both matrices simultaneously.

---

## Architecture

The module follows the **CodyFortran architectural paradigm**:

```
M_Purification.f90          ← Interface module (procedure pointers)
S_Purification.f90          ← Submodule (fabrication dispatcher)
    │
    ├── Singlet/            ← Singlet spin-symmetry branch
    │   ├── M_Purification_Singlet.f90
    │   ├── S_Purification_Singlet.f90
    │   └── Reduce/
    │       ├── M_Purification_Singlet_Reduce.f90
    │       └── S_Purification_Singlet_Reduce.f90  ← IMPLEMENTATION
    │
    ├── SingletTr0/         ← Singlet with trace-free representation
    │   ├── M_Purification_SingletTr0.f90
    │   ├── S_Purification_SingletTr0.f90
    │   └── Reduce/
    │       ├── M_Purification_SingletTr0_Reduce.f90
    │       └── S_Purification_SingletTr0_Reduce.f90  ← IMPLEMENTATION
    │
    └── Generic/            ← Spin-agnostic branch
        ├── M_Purification_Generic.f90
        ├── S_Purification_Generic.f90
        └── Alcoba/
            ├── M_Purification_Generic_Alcoba.f90
            └── S_Purification_Generic_Alcoba.f90  ← STUB (not implemented)
```

### Key Patterns

1. **Interface Modules (`M_*.f90`)**: Define procedure pointer contracts. Only `M_*` files are imported by external code.

2. **Submodules (`S_*.f90`)**: Contain implementations. Each `*_Fabricate` routine reads JSON and binds procedure pointers.

3. **Hierarchical Dispatch**: `Purification_Fabricate` → `Singlet_Fabricate` → `Singlet_Reduce_Fabricate`

---

## JSON Configuration

```json
{
  "purification": {
    "singlet": {
      "reduce": {
        "nSteps": 40,
        "threshold": -1e-5,
        "conserveEnergyQ": false
      }
    }
  }
}
```

### Configuration Paths

| JSON Path                                        | Type    | Default  | Description                                |
|-------------------------------------------------|---------|----------|--------------------------------------------|
| `purification.singlet.reduce.nSteps`            | int     | 40       | Max iterations                             |
| `purification.singlet.reduce.threshold`         | real    | -1e-5    | Eigenvalue cutoff (negative = some allowed)|
| `purification.singlet.reduce.conserveEnergyQ`   | bool    | false    | Project corrections ⊥ to h² (energy cons.) |
| `purification.singletTr0.reduce.*`              | —       | (same)   | Same parameters for SingletTr0 backend     |
| `purification.generic.alcoba`                   | —       | —        | Enable Alcoba backend (stub)               |

---

## Algorithm: Iterative Defect Reduction

The "reduce" algorithm (used by both `Singlet` and `SingletTr0`) operates as follows:

```
1. UNPACK rdm2ab from twoRdmState

2. LOOP until min(eigenvalues) > threshold OR iteration > nSteps:

   a. DECOMPOSE rdm2ab into spin-adapted blocks:
      • rdm2singlet (symmetric under exchange)
      • rdm2triplet (antisymmetric under exchange)

   b. For each block:
      • Eigendecompose to find negative eigenvalues
      • Extract defective part (sum of λ_i |ψ_i⟩⟨ψ_i| for λ_i < threshold)
      • Make defect trace-free (preserves particle number)

   c. SUBTRACT defects from rdm2ab

   d. COMPUTE Q-RDM from updated rdm2ab:
      qrdm2ab = hole-RDM derived from rdm2ab

   e. REPEAT steps (a-c) for Q-RDM blocks

   f. OPTIONAL: Project combined correction onto energy-conserving subspace
      (orthogonal complement of h²)

3. REPACK rdm2ab into twoRdmState
```

### Why Both D² and Q²?

Enforcing D² ≥ 0 alone is not sufficient—the correction to D² may cause Q² to become negative. By iteratively correcting both, the algorithm converges to a state satisfying both conditions.

---

## Key Dependencies

| Module                              | Purpose                                      |
|-------------------------------------|----------------------------------------------|
| `M_TwoRdm_Fermionic_Singlet`        | Pack/unpack singlet 2-RDM state              |
| `M_TwoRdm_Fermionic_SingletTr0`     | Pack/unpack trace-free singlet 2-RDM         |
| `M_Utils_RdmHole`                   | Compute Q-RDM from D²                        |
| `M_Utils_RdmAlgebra`                | Spin-block decomposition (singlet/triplet)   |
| `M_Utils_RdmDefective`              | Extract negative eigenspace                  |
| `M_Utils_RdmEnforceTraceM2sym`      | Make symmetric 2-RDM trace-free              |
| `M_Utils_RdmEnforceTraceM2anti`     | Make antisymmetric 2-RDM trace-free          |
| `M_Utils_RdmEnforceTraceM2ab`       | Make α-β 2-RDM trace-free (SingletTr0)       |
| `M_Utils_RdmProjection`             | Project onto energy-conserving subspace      |

---

## Usage Example

```fortran
use M_Purification

! During initialization (after JSON loaded)
call Purification_Fabricate
call Purification_Setup

! During TD-2RDM propagation
call Propagator_Propagate(state, t0, t1)
call Purification_PurifyTwoRdmState(twoRdmState, h2)
```

---

## Extending the Module

### Adding a New Backend

1. Create directory structure: `Purification/<Family>/<Method>/`

2. Create interface module `M_Purification_<Family>_<Method>.f90`:
   ```fortran
   module M_Purification_<Family>_<Method>
     interface
       module subroutine Purification_<Family>_<Method>_Fabricate
       end subroutine
     end interface
   end module
   ```

3. Create submodule `S_Purification_<Family>_<Method>.f90`:
   - Read parameters from JSON
   - Bind `Purification_PurifyTwoRdmState` pointer
   - Implement local `PurifyTwoRdmState` subroutine

4. Register in parent dispatcher (`S_Purification_<Family>.f90`):
   ```fortran
   if (Json_GetExistence("purification.<family>.<method>")) then
     call Purification_<Family>_<Method>_Fabricate
   end if
   ```

### Adding a New Family (e.g., Unrestricted)

1. Create `Purification/Unrestricted/` directory tree

2. Add dispatch branch in `S_Purification.f90`:
   ```fortran
   else if (Json_GetExistence("purification.unrestricted")) then
     call Purification_Unrestricted_Fabricate
   ```

---

## Diagnostics

The algorithm prints progress to stdout:

```
 PURIFICATION start with min eval -0.123456E-03
 PURIFICATION ended with min eval -0.987654E-06 at step     12
```

- **Start message**: Only printed if initial eigenvalues violate threshold
- **End message**: Always printed with final min eigenvalue and iteration count

---

## Known Limitations

1. **Generic/Alcoba**: Currently a stub—does not perform any purification.

2. **Thread Safety**: Module-level procedure pointers are not thread-safe. Do not call `Purification_Fabricate` from multiple threads.

3. **Energy Conservation**: The `conserveEnergyQ` option projects corrections orthogonal to h², which approximately conserves energy but may slow convergence.

---

## References

- Mazziotti, D.A., "Contracted Schrödinger equation", Phys. Rev. A 57, 4219 (1998)
- Alcoba, D.R. et al., "Purification of correlated reduced density matrices", J. Chem. Phys. 114, 2047 (2001)
