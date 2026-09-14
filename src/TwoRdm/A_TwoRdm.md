# TwoRdm Module - Agent Onboarding Guide

## Overview

The **TwoRdm** module implements the two-body reduced density matrix (2-RDM) subsystem for quantum many-body simulations. The 2-RDM captures all two-particle correlations in an N-particle system and is central to time-dependent 2-RDM (TD-2RDM) propagation methods.

### What is a 2-RDM?

Given an N-particle wavefunction |Ψ⟩, the 2-RDM is defined as:

```
D²(i₁,i₂;j₁,j₂) = ⟨Ψ| a†_{j₁} a†_{j₂} a_{i₂} a_{i₁} |Ψ⟩
```

Key properties:
- **Hermitian**: D²(i₁,i₂;j₁,j₂) = [D²(j₁,j₂;i₁,i₂)]*
- **Antisymmetric** (fermions): D²(i₁,i₂;j₁,j₂) = -D²(i₂,i₁;j₁,j₂)
- **Normalization**: Tr(D²) = N(N-1)
- **Trace relation**: The 1-RDM is obtained by partial trace: D¹ ∝ Tr₂(D²)

---

## Architecture

### Module Hierarchy

```
M_TwoRdm (Interface Module)
│   └── S_TwoRdm (Submodule: Fabrication)
│
└── Fermionic/
    ├── M_TwoRdm_Fermionic (Backend Selector)
    │   └── S_TwoRdm_Fermionic (Submodule: Normalization)
    │
    ├── Generic/
    │   ├── M_TwoRdm_Fermionic_Generic
    │   └── S_TwoRdm_Fermionic_Generic
    │
    ├── Singlet/
    │   ├── M_TwoRdm_Fermionic_Singlet
    │   └── S_TwoRdm_Fermionic_Singlet
    │
    └── SingletTr0/
        ├── M_TwoRdm_Fermionic_SingletTr0
        └── S_TwoRdm_Fermionic_SingletTr0
```

### Design Pattern

The module follows CodyFortran's **Interface/Implementation** pattern:

1. **Interface modules** (`M_*.f90`) declare procedure pointer interfaces and shared state
2. **Implementation submodules** (`S_*.f90`) provide concrete implementations
3. **Fabrication** reads JSON config and binds procedure pointers at runtime

This allows swapping backends (e.g., Generic ↔ Singlet) without recompiling dependent code.

---

## Available Backends

| Backend | JSON Path | Storage Size | Use Case |
|---------|-----------|--------------|----------|
| **Generic** | `twoRdm.fermionic.generic` | nOrbs⁴ | Unrestricted, no spin symmetry |
| **Singlet** | `twoRdm.fermionic.singlet` | nα² × nβ² | Spin-singlet (S=0), restricted |
| **SingletTr0** | `twoRdm.fermionic.singletTr0` | nα² + nα²×nβ² | Singlet with explicit 1-RDM |

### When to Use Each Backend

- **Generic**: Broken-symmetry states, benchmarking, or when spin symmetry cannot be assumed
- **Singlet**: Closed-shell systems, large orbital spaces where memory matters
- **SingletTr0**: TD-2RDM propagation where maintaining trace consistency is critical

---

## Key Data Structures

### Module Variables (M_TwoRdm)

```fortran
integer(I32) :: TwoRdm_nElements
! Total number of elements in linearized 2-RDM

complex(R64), contiguous, pointer :: TwoRdm_twoRdmState(:)
! Linearized 2-RDM storage (layout depends on backend)
```

### Linearized Storage Layouts

**Generic Backend:**
```
twoRdmState(count) = rdm2(i₁, i₂, j₁, j₂)
count = i₁ + nO*(i₂-1) + nO²*(j₁-1) + nO³*(j₂-1)
```

**Singlet Backend:**
```
twoRdmState(count) = rdm2ab(i₁, i₂, j₁, j₂)
count = i₁ + nα*(i₂-1) + nα*nβ*(j₁-1) + nα*nβ*nα*(j₂-1)
```
(stores only αβ block; same-spin blocks reconstructed via singlet relations)

**SingletTr0 Backend:**
```
twoRdmState = [ rdm1a(nα²) | rdm2abTr0(nα²×nβ²) ]
```
(explicit 1-RDM + trace-free 2-RDM)

---

## Core Procedures

### Fabrication & Setup

```fortran
call TwoRdm_Fabricate()  ! Wire procedure pointers from JSON
call TwoRdm_Setup()      ! Allocate storage (backend-specific)
```

### Commutator Operations (for time propagation)

```fortran
call TwoRdm_ApplyCommutatorH1(dTwoRdmState, twoRdmState, h1)
! Accumulates [ĥ₁, D²] into dTwoRdmState

call TwoRdm_ApplyCommutatorH2(dTwoRdmState, twoRdmState, h2)
! Accumulates [Ŵ₁₂, D²] into dTwoRdmState
```

These are central to the BBGKY hierarchy equation of motion:
```
iℏ ∂D²/∂t = [Ĥ₁+Ĥ₂, D²] + [Ŵ, D³]
                           ↑ handled by Reconstruction module
```

### Conversion Routines

```fortran
call TwoRdm_FillRdm1FromTwoRdmState(rdm1, twoRdmState)
! Compute 1-RDM from 2-RDM via partial trace

call TwoRdm_FillRdm2FromTwoRdmState(rdm2, twoRdmState)
! Unpack linearized → full 4-index tensor

call TwoRdm_FillTwoRdmStateFromRdm2(twoRdmState, rdm2)
! Pack full tensor → linearized storage
```

### Normalization

```fortran
call TwoRdm_Normalize(twoRdmState)
! Rescale to enforce Tr(D²) = N(N-1)
```

### I/O

```fortran
call TwoRdm_SaveTwoRdmState(twoRdmState)
! Write to 'twoRdmState.dat' in binary format
```

---

## JSON Configuration

```json
{
  "twoRdm": {
    "fermionic": {
      "generic": {}
    }
  }
}
```

Or for singlet systems:
```json
{
  "twoRdm": {
    "fermionic": {
      "singlet": {}
    }
  }
}
```

Or for trace-free singlet:
```json
{
  "twoRdm": {
    "fermionic": {
      "singletTr0": {}
    }
  }
}
```

---

## Mathematical Details

### One-Body Commutator

```
[ĥ₁, D²]_{i₁i₂;j₁j₂} = Σ_k ( h_{i₁k}D²_{ki₂;j₁j₂} + h_{i₂k}D²_{i₁k;j₁j₂}
                            - h_{kj₁}D²_{i₁i₂;kj₂} - h_{kj₂}D²_{i₁i₂;j₁k} )
```

### Two-Body Commutator

```
[Ŵ₁₂, D²]_{i₁i₂;j₁j₂} = 2·Σ_{k₁k₂} ( W_{i₁i₂;k₁k₂}D²_{k₁k₂;j₁j₂}
                                     - W_{k₁k₂;j₁j₂}D²_{i₁i₂;k₁k₂} )
```

### Partial Trace (1-RDM from 2-RDM)

```
D¹(i,j) = (1/(N-1)) · Σ_k D²(i,k;j,k)
```

### Singlet Relations

For spin-singlet states:
- D²_αα = D²_ββ (same-spin blocks equal)
- D²_αβ = D²_βα (mixed-spin blocks equal)
- D¹_α = D¹_β (spin densities equal)

---

## Dependencies

### Required Modules

- `M_Utils_Types` - Type definitions (R64, I32, etc.)
- `M_Utils_Json` - JSON configuration parsing
- `M_Utils_Say` - Logging/output
- `M_Utils_DataStorage` - Binary I/O
- `M_Method_Mb` - Many-body method parameters
- `M_Method_Mb_OrbBased` - Orbital counts

### Singlet-Specific Dependencies

- `M_Utils_RdmTrace` - Partial trace operations
- `M_Utils_RdmSingletRelations` - Singlet symmetry reconstruction
- `M_Utils_RdmEnforceTraceM2ab` - Trace-free enforcement (SingletTr0)

---

## Extending the Module

### Adding a New Backend

1. Create directory: `src/TwoRdm/Fermionic/NewBackend/`
2. Create interface module `M_TwoRdm_Fermionic_NewBackend.f90`
3. Create submodule `S_TwoRdm_Fermionic_NewBackend.f90`
4. Implement all procedures from the `M_TwoRdm` interface
5. Add fabrication call in `S_TwoRdm_Fermionic.f90`
6. Add JSON key check for `twoRdm.fermionic.newBackend`

### Implementation Checklist

Your backend must implement:
- [ ] `TwoRdm_nElements` calculation
- [ ] `ApplyCommutatorH1`
- [ ] `ApplyCommutatorH2`
- [ ] `FillRdm1FromTwoRdmState`
- [ ] `FillRdm2FromTwoRdmState`
- [ ] `FillTwoRdmStateFromRdm2`

---

## Testing

Run module tests with:
```bash
cd build && ctest -R T_TwoRdm
```

Or use the full test suite:
```bash
bash runTestsAndCov.sh
```

---

## Common Pitfalls

1. **Commutators accumulate**: Always pre-initialize `dTwoRdmState` before calling commutator routines

2. **Conversion allocates**: `FillRdm1FromTwoRdmState` and `FillRdm2FromTwoRdmState` allocate their output arrays

3. **Backend mismatch**: The linearized storage layout is backend-specific; don't mix states from different backends

4. **Index convention**: This code uses physicist's convention where creation operators precede annihilation in RDM definition

5. **Singlet assumption**: Singlet/SingletTr0 backends assume nα = nβ (equal α and β orbital counts)
