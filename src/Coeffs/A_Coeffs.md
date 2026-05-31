# Coeffs Module — Agent Onboarding Guide

## Overview

The **Coeffs** module manages Configuration Interaction (CI) coefficient vectors—the
expansion amplitudes of a many-body wave function in a Fock-space basis:

```
|Ψ⟩ = Σ_I c_I |I⟩
```

where `|I⟩` is a Slater determinant (fermions) or permanent (bosons).

This module is central to all many-body dynamics (TDCI, MCTDHX, TD-2RDM) as it:
1. Stores and manipulates CI coefficient vectors
2. Applies one-body (H1) and two-body (H2) Hamiltonians
3. Computes reduced density matrices (RDMs) for observables
4. Handles index mapping between configurations and linear arrays

---

## Architecture

```
┌─────────────────────────────────────────────────────────────────────────┐
│                         M_Coeffs (Interface)                            │
│  • Procedure pointer declarations                                       │
│  • Module data: Coeffs_nCoeffs, Coeffs_coeffs(:)                       │
└────────────────────────────────┬────────────────────────────────────────┘
                                 │
                    S_Coeffs (Fabricate + shared routines)
                                 │
         ┌───────────────────────┴───────────────────────┐
         │                                               │
         ▼                                               ▼
┌─────────────────────┐                     ┌─────────────────────────┐
│  M_Coeffs_Generic   │                     │   M_Coeffs_Hubbard      │
│  (tensor-product)   │                     │   (bit-encoded)         │
└─────────┬───────────┘                     └───────────┬─────────────┘
          │                                             │
          ▼                                   ┌─────────┴─────────┐
┌─────────────────────┐                       │                   │
│ S_Coeffs_Generic    │           ┌───────────┴───┐   ┌───────────┴───┐
│ • ApplyH1FillRdm1   │           │  NoSpinSym    │   │  PlusSpinSym  │
│ • ApplyH2FillRdm2   │           │  n_up × n_dn  │   │  n(n+1)/2     │
│ • ApplyExcitation   │           └───────────────┘   └───────────────┘
│ • Index mappings    │                     │
└─────────────────────┘           ┌─────────┴─────────┐
                                  │   MinusSpinSym    │
                                  │   n(n-1)/2        │
                                  └───────────────────┘
```

### Design Pattern

The module follows CodyFortran's **M_/S_ interface-submodule pattern**:

| File | Role |
|------|------|
| `M_Coeffs.f90` | Interface module: declares procedure pointers, abstract interfaces, module data |
| `S_Coeffs.f90` | Submodule: implements `Fabricate` and backend-independent routines |
| `M_Coeffs_*.f90` | Backend interface modules |
| `S_Coeffs_*.f90` | Backend implementation submodules |

**Runtime binding:** `Coeffs_Fabricate` reads JSON configuration and binds procedure
pointers to the appropriate backend. This enables swapping implementations without
recompiling dependent code.

---

## Backends

### 1. Generic (`coeffs.generic`)

**Use case:** General many-body systems (atoms, molecules, custom models)

**Basis construction:**
- Tensor product of per-body-type ConfigLists
- Total dimension: `∏_{bt} nConfigurations(bt)`
- Supports fermions, bosons, or mixtures

**Index mapping:**
```
iCoeff = 1 + Σ_bt (C_bt - 1) × stride_bt
stride_bt = ∏_{bt' < bt} nConfigurations(bt')
```

**Key features:**
- Uses precomputed `singles`/`doubles` connectivity graphs from ConfigList
- OpenMP parallelization over CI indices
- Linearized H1/H2 storage for cache efficiency

### 2. Hubbard (`coeffs.hubbard`)

**Use case:** Fermi-Hubbard model on lattice grids

**Basis construction:**
- Bit-encoded Fock states: bit j=1 means site j occupied
- Separate encoding for spin-up and spin-down
- Three spin-symmetry variants reduce state space

**Spin symmetry variants:**

| Variant | Dimension | Physical meaning |
|---------|-----------|------------------|
| `noSpinSym` | n_up × n_dn | Full product space |
| `plusSpinSym` | n(n+1)/2 | Symmetric under spin exchange (triplet-like) |
| `minusSpinSym` | n(n-1)/2 | Antisymmetric under spin exchange (singlet-like) |

**Key features:**
- Precomputed hopping graphs for O(1) kinetic term lookup
- Bit operations (`popcnt`, `btest`) for O(1) sign computation
- Diagonal interaction: `interactionValues` indexed by `bitcodesUP AND bitcodesDN`

---

## Key Procedures

### Hamiltonian Application

| Procedure | Description |
|-----------|-------------|
| `Coeffs_ApplyH1FillRdm1` | Apply Ĥ₁ = Σ h¹_{pq} a†_p a_q; optionally compute 1-RDM |
| `Coeffs_ApplyH2FillRdm2` | Apply Ĥ₂ = ½Σ h²_{pqrs} a†_p a†_q a_s a_r; optionally compute 2-RDM |

These are the performance-critical routines. They simultaneously:
1. Compute `dCoeffs = H|coeffs⟩` for time propagation
2. Accumulate RDM elements for observables

### RDM Computation

| Procedure | Output | Use case |
|-----------|--------|----------|
| `Coeffs_FillRdm1Bt` | ρ_{pq} = ⟨a†_p a_q⟩ | Natural orbitals, occupations |
| `Coeffs_FillRdm2Bt` | ρ_{pqrs} = ⟨a†_p a†_q a_s a_r⟩ | Pair correlations |
| `Coeffs_FillRdm3Bt` | ρ_{pqrstu} = ⟨a†_p a†_q a†_r a_u a_t a_s⟩ | BBGKY closure |

### Ladder Operators

```fortran
call Coeffs_ApplyExcitation(coeffs, creates, destroys, bt)
```

Applies `a†_{creates(1)} a†_{creates(2)} ... a_{destroys(n)} ... a_{destroys(1)}`
to body type `bt`. Used for:
- Computing RDM elements
- Excited state preparation
- Transition dipole moments

### Index Mapping

```fortran
call Coeffs_ConfigurationsFromIndex(configurations, iCoeff)  ! Linear → tuple
call Coeffs_IndexFromConfigurations(iCoeff, configurations)  ! Tuple → linear
```

---

## JSON Configuration

### Generic backend
```json
{
  "coeffs": {
    "generic": { }
  }
}
```

### Hubbard backend
```json
{
  "coeffs": {
    "hubbard": {
      "noSpinSym": { }
    }
  }
}
```

Replace `noSpinSym` with `plusSpinSym` or `minusSpinSym` for symmetry-adapted bases.

---

## Initialization Sequence

```fortran
! 1. Fabricate: read JSON, bind procedure pointers, compute nCoeffs
call Coeffs_Fabricate

! 2. Setup: build configuration tables, hopping graphs, etc.
call Coeffs_Setup

! 3. Initialize: set initial values (via CoeffsInit module)
call CoeffsInit_Initialize(Coeffs_coeffs)

! 4. Use
call Coeffs_ApplyH1FillRdm1(coeffs, rdm1_=rdm1, dCoeffs_=dCoeffs, h1_=h1)
call Coeffs_Normalize(coeffs)
```

---

## Data Structures

### Module-level (M_Coeffs)

| Variable | Type | Description |
|----------|------|-------------|
| `Coeffs_nCoeffs` | `integer(I32)` | Total number of CI coefficients |
| `Coeffs_coeffs` | `complex(R64), pointer` | Pointer to coefficient segment in state vector |

### Hubbard-specific (M_Coeffs_Hubbard)

| Variable | Description |
|----------|-------------|
| `bitcodesUP(:)` | Bit patterns for spin-up configurations |
| `bitcodesDN(:)` | Bit patterns for spin-down configurations |
| `hoppUP(k,i)` | Target config index for k-th hop from config i |
| `weightUP(k,i)` | Matrix element (including fermionic sign) |
| `nConnectedUP(i)` | Number of valid hops from config i |
| `interactionValues(:)` | Precomputed U·n↑n↓ for each (up,dn) pair |

---

## Performance Considerations

1. **OpenMP parallelization:** Loops over CI indices use `!$omp parallel do` with
   thread-local accumulators for RDM construction.

2. **Linearized storage:** H1/H2 tensors are packed into 1D arrays (`h1Lin`, `h2Lin`)
   for cache-friendly access during Hamiltonian application.

3. **Symmetry exploitation:** 
   - Hubbard backends reduce state space by 2× with spin symmetry
   - Generic backend uses precomputed sparse connectivity graphs

4. **Buffer reuse:** Generic `ApplyExcitation` uses a module-level `coeffsBuffer`
   to avoid repeated allocation.

---

## File Index

| File | Purpose |
|------|---------|
| `M_Coeffs.f90` | Interface: procedure pointers, abstract interfaces |
| `S_Coeffs.f90` | Fabricate + shared routines (normalize, project, RDM fills) |
| `Generic/M_Coeffs_Generic.f90` | Generic backend interface |
| `Generic/S_Coeffs_Generic.f90` | Generic backend implementation |
| `Hubbard/M_Coeffs_Hubbard.f90` | Hubbard backend interface + module data |
| `Hubbard/S_Coeffs_Hubbard.f90` | Hubbard setup, bit operations, excitations |
| `Hubbard/NoSpinSym/*.f90` | Full product space (n_up × n_dn) |
| `Hubbard/PlusSpinSym/*.f90` | Symmetric sector n(n+1)/2 |
| `Hubbard/MinusSpinSym/*.f90` | Antisymmetric sector n(n-1)/2 |

---

## Common Tasks

### Add a new backend

1. Create `M_Coeffs_MyBackend.f90` with `Coeffs_MyBackend_Fabricate` interface
2. Create `S_Coeffs_MyBackend.f90` implementing all required procedures
3. Add JSON branch in `S_Coeffs.f90::Coeffs_Fabricate`

### Debug RDM computation

Use the body-type-resolved routines (`FillRdm1Bt`, etc.) which are simple
brute-force loops. Compare against the optimized `ApplyH1FillRdm1` output.

### Checkpoint/restart

```fortran
call Coeffs_SaveCoeffs(coeffs)  ! Writes coeffs.in
! Later...
! Use CoeffsInit with "load": {} to reload
```

---

## Dependencies

| Module | Purpose |
|--------|---------|
| `M_ConfigList` | Fock-space basis enumeration (configurations) |
| `M_Method_Mb` | Many-body method parameters (nBodyTypes, nBodies, statistics) |
| `M_Method_Mb_OrbBased` | Orbital indexing (nOrbs, orbital ranges) |
| `M_Utils_BlasLib` | BLAS wrappers for normalization |
| `M_Utils_Json` | JSON configuration parsing |
| `M_Grid` / `M_Grid_Lattice` | Spatial grid (Hubbard backend) |
