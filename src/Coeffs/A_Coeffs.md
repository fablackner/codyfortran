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
- Optional restriction to a total-M sector (below)

**Total-M sector (`coeffs.generic.totalM`, optional):** with declared orbital
m (`method.mb.orbBased.mOfOrbs`, required) every configuration of a body type
has a definite M = Σ m_i n_i, and h1, h2 conserve the total M (z-polarized
field). `Coeffs_Setup` (bound by the Generic backend) then builds
- `sector`: the indices of the configurations of total M = `totalM`,
- `sectorBt(bt)`: the bucket bounds of the ConfigList connection tables of each
  body type, which Setup sorts in place (stably) by the change d of M they
  cause: singles `singlesFirst(d, iC)` .. `singlesFirst(d+1, iC)-1`, doubles
  that conserve M `doublesFirst(0, iC)` .. `doublesFirst(1, iC)-1` (the order
  of the connections is not part of the ConfigList contract; no copy, so the
  connectivity memory is that of the full space).

`ApplyH1FillRdm1` and `ApplyH2FillRdm2` run `!$omp do` over the sector, take
the singles of bucket 0 and the M-conserving doubles, and pair a single of
change d on one body type with a single of change −d on the other, so no
connection leaves the sector. Storage stays full length (`Coeffs_nCoeffs`,
`coeffs.in` and the integrators unchanged); rows outside the sector are never
touched, and the result is exact for vectors in the sector,
which `Coeffs_ProjectOnSector` enforces once at initialization (called by
`CoeffsInit_Initialize`, and on the CI eigenvectors of the Ylm-optimized MCSCF
ground solver). It stops with an error if more than 1e-8 of the norm lies
outside the sector. On full-length vectors the Hamiltonian is thus
H_sector ⊕ 0: a diagonalizer, the eigen-expansion propagator or a CI spectrum
of length `Coeffs_nCoeffs` sees `Coeffs_nCoeffs` − |sector| spurious
eigenvalues 0 (11760 for Ne below), harmless only for a ground state of
negative energy (MCSCF of atoms); excited or positive-energy spectra need the
key removed. Without the key the sector holds all configurations and one
bucket all connections in ConfigList order, i.e. the loops of the full
space. Ne 9 orbitals (5+5 electrons): the M = 0 sector has 4116 of 15876
configurations; 16 steps of the production setup take 41.2 s instead of
45.4 s (−0.26 s per step). Tests: `test/componentTests/Coeffs/Generic/T_Coeffs_Generic_Sector`
(sector of total M = 1 against the full space for two fermionic and one
bosonic body type: rdm1, rdm2 and H|c⟩ of random h1, h2, exact zeros outside
the sector, `Coeffs_ProjectOnSector`), `T_Coeffs_Generic_SectorOutside` (a
unary configuration outside the sector stops with the error), and
`test/simulationTests/Ne3d/T_Ne3d_15/16` (MCSCF and MCTDHF laser propagation
of Be-like Ne6+, 34 of 100 configurations, against the full space).

### 2. Hubbard (`coeffs.hubbard`)

CI of two fermion species (spin up and spin down, body types 1 and 2) in the
same n orbitals, for a two-body Hamiltonian that is an on-site interaction
Σ_s U_s n_{s↑} n_{s↓} in the orbital basis: the Fermi-Hubbard model with
the lattice sites as orbitals (`orbsInit.lattice.onSite`) and static
orbitals (`method.mb.orbBased.tdci`). Bosons and other interactions use
`coeffs.generic`. No `configList` is needed.

**Basis:** a configuration of one spin is a bit pattern (bit o − 1 for
orbital o, at most 62 orbitals). The configurations are numbered in
increasing order of the pattern, so the index of a pattern is one plus its
colex rank Σ_t binom(position of the t-th set bit, t), without a lookup
table. The CI vector is the matrix c(iUp, iDn), up index fastest. An
operator on orbital o picks up (−1)^(occupied orbitals below o).

**Setup:** `Coeffs_Setup` tabulates the single excitations
a†_c a_d |i⟩ = ±|t⟩ of every configuration i of each spin (d occupied, c
empty or c = d), N(n − N) + N entries per configuration. The tables do not
depend on the Hamiltonian.

**H1** comes from the `h1_` of each call and acts on the rows (up) and the
columns (down) of c(iUp, iDn). Each call sums the diagonal entries of a
configuration into one element and keeps only the nonzero off-diagonal ones
(the hops of the lattice). A time-dependent
h1, the mean-field reference of MCTDHX or the basis change of
`Method_Mb_OrbBased_OrthonormalizeKeepingState` act like any other h1.

**H2** is diagonal, with U_s read from `h2_` on every call (adiabatic
switching scales it). U_s is the combination of h2 elements that
`coeffs.generic` uses for two body types, at a = b = g = d = s. The call
stops if that combination does not vanish for every other (a, b, g, d), or
if the antisymmetrized interaction within one spin does not vanish:
orbitals that are not the sites, moving orbitals or a nearest-neighbour
interaction.

**RDMs** come from chains of single excitations, in the conventions of
`coeffs.generic`, including the antisymmetry under the exchange of an up and
a down index (`FillRdm2FromRdm2Lin`).

**Spin symmetry variants:**

| Variant | Coefficients | Stored vector |
|---------|--------------|---------------|
| `noSpinSym` | n_up n_dn | c(iUp, iDn) |
| `plusSpinSym` | n(n+1)/2 | (\|i,k⟩ + \|k,i⟩)/√2 for i < k, and \|i,i⟩ |
| `minusSpinSym` | n(n−1)/2 | (\|i,k⟩ − \|k,i⟩)/√2 for i < k |

The sectors are the eigenspaces of the exchange of the two spins. They need
equal particle numbers and an h1 that is the same for both spins (checked at
every call). Their bases are orthonormal, so norms, inner products,
diagonalizers and propagators act on the stored vectors as on the states.
H and the RDMs of a sector vector are computed on the product space and
projected back. A sector has no excitations of one spin:
`Coeffs_ApplyExcitation`, `coeffsInit.excited` and the transition RDMs
stop. For 1 + 1 particles the plus sector holds the singlets (symmetric
spatial part).

**Cost** (open chain at half filling, release build, one thread of an
i7-8565U; H1 + H2 is one TDCI time derivative):

| Sites | Coefficients | Backend | Setup | H1 + H2 | 2-RDM |
|-------|--------------|---------|-------|---------|-------|
| 8 | 4 900 | generic | 15 ms | 11.6 ms | 9.8 ms |
| 8 | 4 900 | hubbard | 0.07 ms | 0.13 ms | 9.9 ms |
| 10 | 63 504 | generic | 137 ms | 323 ms | 308 ms |
| 10 | 63 504 | hubbard | 0.2 ms | 1.4 ms | 274 ms |
| 12 | 853 776 | generic | 1.1 s | 10.3 s | 10.2 s |
| 12 | 853 776 | hubbard | 1 ms | 23 ms | 7.7 s |

The 2-RDM runs over pairs of single excitations in both backends. A sector
stores half the vector, and one application costs about twice that of the
product space (expansion and projection; 48 ms at 12 sites).

**Tests:** `test/componentTests/Coeffs/Hubbard` compares noSpinSym with
coeffs.generic on random Hamiltonians (spectrum, ground-state RDMs,
excitation-based 2-RDM blocks), the sectors with the product space, and
checks the refusal of an off-site interaction.
`test/simulationTests/FermiHubbard/T_FermiHubbard_07` to `_12` run exact
diagonalization, imaginary time and real time on lattices against the
coeffs.generic values, real time also in the plus sector.

---

## Key Procedures

### Hamiltonian Application

| Procedure | Description |
|-----------|-------------|
| `Coeffs_ApplyH1FillRdm1` | Apply Ĥ₁ = Σ h¹_{pq} a†_p a_q; optionally compute 1-RDM |
| `Coeffs_ApplyH2FillRdm2` | Apply Ĥ₂ = Σ h²_{pqrs} a†_p a†_q a_s a_r with h²_{pqrs} = ½⟨pq\|W\|rs⟩; optionally compute 2-RDM |
| `Coeffs_ProjectOnSector` | Zero the coefficients outside the declared total-M sector (NoOp without `coeffs.generic.totalM`) |

Index convention (shared with `Method_Mb_OrbBased_FillH1/H2`, TwoRdm and the
orbital equations of motion): h¹_{pq} = ⟨p|h|q⟩, h²_{pqrs} = ½⟨pq|W|rs⟩ (the ½
is part of h2, as filled by `Method_Mb_OrbBased_FillH2`), rdm1(p,q) = ⟨a†_q a_p⟩ and
rdm2(p,q,r,s) = ⟨a†_r a†_s a_q a_p⟩. The Generic backend's connectivity loops
act with the transposed index order internally; the packing helpers
(`FillH1Lin`, `FillH2Lin`, `FillRdm{1,2}FromRdm{1,2}Lin`) translate. The
distinction only matters for complex orbitals or Hamiltonians (e.g., the
velocity-gauge laser coupling), where the transposed matrices are the complex
conjugates.

These are the performance-critical routines. They simultaneously:
1. Compute `dCoeffs = H|coeffs⟩` for time propagation
2. Accumulate RDM elements for observables

### RDM Computation

The 1- and 2-RDM come with `Coeffs_ApplyH1FillRdm1` and
`Coeffs_ApplyH2FillRdm2`; body-type blocks of 2- and 3-RDMs and transition
RDMs from `Coeffs_FillTransitionRdm{2,3}Bt` (below).

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
    "generic": { "totalM": 0 }
  }
}
```

`totalM` is optional (requires `method.mb.orbBased.mOfOrbs`); `"generic": { }`
uses the full space.

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

Replace `noSpinSym` with `plusSpinSym` or `minusSpinSym` for the spin-exchange
sectors. No `configList` block.

---

## Initialization Sequence

```fortran
! 1. Fabricate: read JSON, bind procedure pointers, compute nCoeffs
call Coeffs_Fabricate

! 2. Setup: build configuration tables, hopping graphs, the Generic sector, etc.
call Coeffs_Setup

! 3. Initialize: set initial values (via CoeffsInit module; ends with
!    Coeffs_ProjectOnSector)
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
| `Coeffs_Hubbard_nSites` | Orbitals of each spin |
| `Coeffs_Hubbard_nCoeffsUP`, `_nCoeffsDN` | Configurations of each spin, binom(n, N) |
| `Coeffs_Hubbard_parity` | 0 product space, +1 / −1 spin-exchange sector (set by the variant) |
| `Coeffs_Hubbard_bitcodesUP(:)`, `_bitcodesDN(:)` | Bit patterns of the configurations, increasing |
| `Coeffs_Hubbard_singlesUP`, `_singlesDN` | Single excitations by source configuration: `first(i)`, `target`, `create`, `destroy`, `factor` (±1) |

---

## Performance Considerations

1. **OpenMP parallelization:** Loops over CI indices use `!$omp parallel do` with
   thread-local accumulators for RDM construction.

2. **Linearized storage:** H1/H2 tensors are packed into 1D arrays (`h1Lin`, `h2Lin`)
   for cache-friendly access during Hamiltonian application.

3. **Symmetry exploitation:** 
   - The Hubbard spin-exchange sectors halve the CI dimension
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
| `Hubbard/S_Coeffs_Hubbard.f90` | Hubbard setup, Hamiltonian application, RDMs, sector maps, bit operations |
| `Hubbard/NoSpinSym/*.f90` | Full product space (n_up × n_dn) |
| `Hubbard/PlusSpinSym/*.f90` | Symmetric sector n(n+1)/2 |
| `Hubbard/MinusSpinSym/*.f90` | Antisymmetric sector n(n-1)/2 |

---

## Common Tasks

### Add a new backend

1. Create `M_Coeffs_MyBackend.f90` with `Coeffs_MyBackend_Fabricate` interface
2. Create `S_Coeffs_MyBackend.f90` implementing all required procedures
3. Add JSON branch in `S_Coeffs.f90::Coeffs_Fabricate`

### Excitation-based RDM blocks

One private worker backs every excitation-based block: one excitation of
the bra per element, one matrix-vector product with the kets, OpenMP over
the elements (`Coeffs_ApplyExcitation` must be thread-safe; the buffer of the
generic backend is threadprivate). The cost is n^(2R) excitations of the CI
vector plus the nKets products per element (N2, 10 orbitals, 14 400
coefficients: the reference and 12 transition 3-RDMs in 80 s on three
threads).

| Routine | Convention |
|---------|------------|
| `Coeffs_FillTransitionRdm2Bt(rdm2, bra, kets, bt1, bt2)` | ket-first, rdm2(i1,i2,j1,j2,k) = <bra\| a†_j1 a†_j2 a_i2 a_i1 \|kets(:,k)> (the layout of TwoRdm and `Coeffs_ApplyH2FillRdm2`) |
| `Coeffs_FillTransitionRdm3Bt(rdm3, bra, kets, bt1, bt2, bt3)` | ket-first, as above with three indices each |

To debug, compare a block with bra = ket against `Coeffs_ApplyH2FillRdm2`
(`T_Coeffs_Generic_Sector` does it for a complex vector with fermions and
bosons). The Hubbard spin-exchange sectors (plusSpinSym, minusSpinSym) have
no excitations: an excitation of one spin species leaves the sector.

### Checkpoint/restart

```fortran
call Coeffs_SaveCoeffs(coeffs)  ! Writes coeffs.in
! Later...
! Use CoeffsInit with "load": {} to reload
```

`Coeffs_SaveTwoRdm(coeffs)` writes the 2-RDM blocks `rdm2BtXX_YY.in` (initial
state for `twoRdmInit.load`) via `Method_Mb_OrbBased_SaveRdm2Blocks`.

---

## Dependencies

| Module | Purpose |
|--------|---------|
| `M_ConfigList` | Fock-space basis enumeration (configurations) |
| `M_Method_Mb` | Many-body method parameters (nBodyTypes, nBodies, statistics) |
| `M_Method_Mb_OrbBased` | Orbital indexing (nOrbs, orbital ranges) |
| `M_Utils_BlasLib` | BLAS wrappers for normalization |
| `M_Utils_Json` | JSON configuration parsing |
