# OrbsInit Module — Agent Onboarding Guide

## Overview

The **OrbsInit** module provides orbital initialization for quantum many-body simulations in CodyFortran. It constructs initial single-particle wavefunctions (orbitals) on various spatial discretizations, serving as the starting point for time-dependent and ground-state calculations.

## Architecture

```
M_OrbsInit (Interface)
    │
    └── S_OrbsInit (Fabrication + Generic Initialize)
            │
            ├── Linear/          ─── 1D continuous grids
            │   └── Harmonic/    ─── Quantum HO eigenstates ψ_n(x)
            │
            ├── Lattice/         ─── 3D discrete lattice sites
            │   └── OnSite/      ─── Kronecker δ orbitals (Hubbard models)
            │
            ├── GridPoint/       ─── Generic δ-basis on any grid
            │
            ├── Load/            ─── File-based orbital loading
            │
            ├── Ylm/             ─── Spherical harmonics R(r)·Y_lm
            │   └── HydrogenLike/─── Hydrogenic radial functions R_nl(r)
            │
            ├── Prolate/         ─── Two-center prolate grid
            │   └── Lcao/        ─── gerade/ungerade hydrogenic LCAO
            │
            └── Fem/             ─── 3D spectral elements (grid.fem)
                └── Sto/         ─── Slater-type orbitals with real solid
                                     harmonics on the centers (LCAO)
```

### Design Patterns

1. **Interface/Implementation Split**: Each backend has:
   - `M_*.f90` — Module interface with procedure pointers and abstract interfaces
   - `S_*.f90` — Submodule containing concrete implementations

2. **Hierarchical Procedure Pointers**: Runtime JSON config drives which implementations get bound:
   ```
   OrbsInit_Fabricate → binds OrbsInit_Initialize
                      → calls OrbsInit_Linear_Fabricate
                                   → binds OrbsInit_InitializeOrb
                                   → calls OrbsInit_Linear_Harmonic_Fabricate
                                                → binds OrbsInit_Linear_InitFunction
   ```

3. **Two-Level Normalization**:
   - Backend `InitFunction` returns unnormalized values
   - `InitializeOrb` samples all grid points, then normalizes via `Grid_InnerProduct`

## Key Interfaces

### Public API (M_OrbsInit)

| Procedure | Signature | Purpose |
|-----------|-----------|---------|
| `OrbsInit_Fabricate` | `()` | JSON-driven backend selection |
| `OrbsInit_Setup` | `()` | Allocate buffers, precompute constants |
| `OrbsInit_Initialize` | `(orbs(:,:))` | Fill all orbitals (nGrid × nOrbs) |
| `OrbsInit_InitializeOrb` | `(orb(:), ind, bt_)` | Initialize single orbital |

### Backend Init Functions

| Backend | Signature | Returns |
|---------|-----------|---------|
| Linear | `InitFunction(x, index, bt_) → real` | Amplitude at position x |
| Lattice | `InitFunction(ix,iy,iz, index, bt_) → real` | Amplitude at site (ix,iy,iz) |
| Ylm | `InitFunction(r, l, m, index, bt_) → complex` | Radial part at (r,l,m) |
| Fem | `OrbsInit_Fem_InitFunction(x, y, z, ind, bt_) → real` | Amplitude at the node (x,y,z); must be thread-safe (OpenMP node loop) |

## JSON Configuration

### Linear / Harmonic Oscillator
```json
{
  "orbsInit": {
    "linear": {
      "harmonic": {
        "position": 0.0,
        "omega": 1.0
      }
    }
  }
}
```
Produces eigenstates: `ψ_n(x) = H_n(ξ)·exp(-ξ²/2)` where `ξ = (x-x₀)/a`, `a = 1/√ω`, `n = index-1`.

### Lattice / On-Site (Hubbard Models)
```json
{
  "orbsInit": {
    "lattice": {
      "onSite": { }
    }
  }
}
```
Produces δ-orbitals: `φ_i(j) = δ_{ij}` with index→site mapping in row-major order.

### Ylm / Hydrogen-Like (Atoms)
```json
{
  "orbsInit": {
    "ylm": {
      "hydrogenLike": {
        "charge": 1.0,
        "n": [1, 2, 2, 2],
        "l": [0, 0, 1, 1],
        "m": [0, 0, 0, 1]
      }
    }
  }
}
```
Produces `R_nl(r) = ρ^l · L_{n-l-1}^{2l+1}(ρ) · exp(-ρ/2)` where `ρ = 2Zr/(na₀)`.
The `m` list only sets the initial orbitals; it is independent of the
declared m of the method (`method.mb.orbBased.mOfOrbs`, optional), which the
initialization then checks (see below).

### Load from File
```json
{
  "orbsInit": {
    "load": { }
  }
}
```
Reads binary files `orb{bt}_{ind}.in` via DataStorage routines. With declared
m, round-off admixtures of other m in the files are removed by the final
projection of `OrbsInit_Initialize` (below).

### GridPoint (δ-basis)
```json
{
  "orbsInit": {
    "gridPoint": { }
  }
}
```
Produces `φ_i(j) = δ_{ij}` normalized via grid inner product.

### Fem / Slater-Type Orbitals (Molecules)

Requires `grid.fem`. The family-level `InitializeOrb` (template:
OrbsInit_Linear) samples `OrbsInit_Fem_InitFunction(x, y, z, ind, bt_)` at
every interior node (OpenMP over z planes outside parallel regions) and
normalizes the orbital with `Grid_InnerProduct`; it does **not**
orthonormalize (programs orthonormalize after `Method_Setup`). An absent
`bt_` (method.sb) is fine; the same orbitals are used for every body type.

The `sto` model takes the terms as parallel arrays (one entry per term);
orbital `ind` is the sum over the terms t with `orbital(t) = ind`:

```
φ_ind(r) = Σ_t coeff_t · ρ^(n_t − l_t − 1) · S_{l_t m_t}(r − R_t) · exp(−ζ_t ρ),
ρ = |r − R_t|,   R_t = Grid_Fem_centers(:, center_t)
```

with the real regular solid harmonics S_lm in Racah normalization
(`M_Utils_SolidHarmonics`: S₀₀ = 1, S₁₀ = z, S₁₁ = x, S₁,₋₁ = y,
S₂₀ = (3z² − r²)/2, S₂,₋₂ = √3 xy, …; m > 0 cos-type, m < 0 sin-type). So
ρ^(n−l−1) S_lm is r^(n−1) times an angular factor — no division by ρ, regular
on the nucleus nodes. The terms are unnormalized primitives (coeff weights
them as they are); the whole orbital is normalized on the grid.

```json
{
  "orbsInit": {
    "fem": {
      "sto": {
        "orbital": [1, 2, 3, 4, 5, 5, 5],
        "center":  [1, 1, 1, 1, 1, 2, 3],
        "n":       [1, 2, 2, 2, 2, 1, 1],
        "l":       [0, 0, 1, 1, 1, 0, 0],
        "m":       [0, 0, 1, 0, -1, 0, 0],
        "zeta":    [7.66, 2.25, 2.23, 2.23, 2.23, 1.0, 1.0],
        "coeff":   [1, 1, 1, 1, 1, 0.3, 0.3]
      }
    }
  }
}
```

| Key | Type | Default | Validation |
|-----|------|---------|------------|
| `orbital` | integer(:) | required | ≥ 1; every index 1..max(orbital) has a term |
| `center` | integer(:) | required | 1 ≤ center ≤ number of `grid.fem` centers |
| `n` | integer(:) | required | n > l |
| `l` | integer(:) | required | l ≥ 0 |
| `m` | integer(:) | required | \|m\| ≤ l |
| `zeta` | real(:) | required | > 0 |
| `coeff` | real(:) | 1 each | — |

All arrays need equal sizes. The model defines `OrbsInit_Fem_nOrbitals =
max(orbital)` orbitals; requesting more (e.g. a larger
`method.mb.orbBased.nOrbs`) stops in `InitializeOrb` with a clear message.

## Initialization Sequence

```fortran
! 1. Fabricate (reads JSON, binds procedure pointers)
call OrbsInit_Fabricate

! 2. Setup (allocates buffers if needed — currently no-op for most backends)
call OrbsInit_Setup

! 3. Initialize orbitals (ends with Orbs_ProjectOnSymmetry(orbs))
call OrbsInit_Initialize(orbs)

! Alternative: initialize single orbital
call OrbsInit_InitializeOrb(orb(:), orbital_index, body_type)
```

`OrbsInit_Initialize` ends with `Orbs_ProjectOnSymmetry(orbs)`: with declared
m (`method.mb.orbBased.mOfOrbs`) the channels of other m are zeroed, and the
program stops with an error if more than 1e-8 of an orbital's norm would be
removed (orbitals that do not have the declared m); without declaration it is
a NoOp. Apart from `EmbedOrbitals` of the Be/Ne programs, which replaces
this initialization, it is the only projection on the declared m
(src/Method/A_Method.md).
`OrbsInit_InitializeOrb` alone does not project.

## Data Flow

```
┌─────────────────────────────────────────────────────────────────┐
│                      JSON Configuration                          │
│  "orbsInit": { "linear": { "harmonic": { "omega": 1.0 } } }     │
└───────────────────────────────┬─────────────────────────────────┘
                                │
                                ▼
┌───────────────────────────────────────────────────────────────────┐
│                     OrbsInit_Fabricate                            │
│  • Reads JSON tree to select backend                              │
│  • Binds OrbsInit_Initialize → Initialize (generic loop)         │
│  • Delegates to backend Fabricate (e.g., Linear → Harmonic)      │
│  • Backend binds OrbsInit_InitializeOrb and InitFunction         │
└───────────────────────────────┬───────────────────────────────────┘
                                │
                                ▼
┌───────────────────────────────────────────────────────────────────┐
│                     OrbsInit_Initialize(orbs)                     │
│  • Loops: do ibt = 1, nBodyTypes; do ind = 1, nOrbs(ibt)         │
│  • Calls OrbsInit_InitializeOrb(orbs(:,i1), ind, ibt)            │
└───────────────────────────────┬───────────────────────────────────┘
                                │
                                ▼
┌───────────────────────────────────────────────────────────────────┐
│              OrbsInit_InitializeOrb (backend-specific)            │
│  • Loops over grid points                                         │
│  • Calls InitFunction(coord, index, bt_) for each point          │
│  • Normalizes via Grid_InnerProduct                              │
└───────────────────────────────────────────────────────────────────┘
```

## Dependencies

| Module | Purpose |
|--------|---------|
| `M_Grid`, `M_Grid_*` | Grid coordinates, inner product, normalization |
| `M_Method_Mb`, `M_Method_Mb_OrbBased` | Body type count, orbital count per type |
| `M_Utils_Json` | JSON configuration parsing |
| `M_Utils_SfGslLib` | Hermite polynomials, Laguerre polynomials |
| `M_Utils_DataStorage` | Binary file I/O for Load backend |
| `M_Grid_Fem`, `M_Utils_SolidHarmonics` | Nodes and centers, real solid harmonics (Fem/Sto) |

## Adding a New Backend

1. **Create directory** under `src/OrbsInit/` (e.g., `PlaneWave/`)

2. **Create interface module** `M_OrbsInit_PlaneWave.f90`:
   ```fortran
   module M_OrbsInit_PlaneWave
     use M_Utils_Types
     implicit none
     interface
       module subroutine OrbsInit_PlaneWave_Fabricate
       end subroutine
     end interface
     ! Optional: expose procedure pointer for sub-backends
   end module
   ```

3. **Create submodule** `S_OrbsInit_PlaneWave.f90`:
   ```fortran
   submodule(M_OrbsInit_PlaneWave) S_OrbsInit_PlaneWave
   contains
     module subroutine OrbsInit_PlaneWave_Fabricate
       use M_OrbsInit
       OrbsInit_InitializeOrb => InitializeOrb
     end subroutine
     
     subroutine InitializeOrb(orb, ind, bt_)
       ! Implementation here
     end subroutine
   end submodule
   ```

4. **Register in S_OrbsInit.f90**:
   ```fortran
   use M_OrbsInit_PlaneWave
   ! ...
   else if (Json_GetExistence("orbsInit.planeWave")) then
     call OrbsInit_PlaneWave_Fabricate
   ```

5. **CMake** collects new sources automatically (`GLOB_RECURSE`).

## Common Pitfalls

1. **Forgetting normalization**: All backends must normalize via `Grid_InnerProduct` which includes proper quadrature weights and metric factors.

2. **Index conventions**: Orbital indices are 1-based. For harmonic oscillator, quantum number `n = index - 1`.

3. **Ylm grid structure**: Grid points encode composite `(r, l, m)` coordinates. The `InitFunction` must return zero for non-matching `(l, m)`.

4. **Body types**: The `bt_` parameter enables species-dependent orbitals (e.g., spin-up vs spin-down). Pass it through even if unused.

## Testing

Tests using OrbsInit are found in:
- `test/T_He1d_*.json` — 1D helium with Linear/Harmonic
- `test/T_FermiHubbard_*.json` — Hubbard model with Lattice/OnSite
- `test/T_H3d_*.json` — 3D hydrogen with Ylm/HydrogenLike
- `test/componentTests/OrbsInit/Load/` — loading with declared m: round-off
  admixtures are removed, orbitals of another m stop with an error
  (`ctest-pass-regex` test)
- `test/componentTests/OrbsInit/Fem/Sto/T_OrbsInit_Fem_Sto` — two-center
  mesh (N = 257594): six orbitals (1s, LCAO of two 1s, 2p_x, 3d_xy, 3s,
  4f (x² − y²)z) normalized to 8e-14 and equal to the analytic functions
  (normalized on the grid) to 1.4e-15; 1s at the nucleus node vs 1/√π:
  1.6e-6 (the grid norm error); p, d, f reflection symmetries to 8e-14; the
  same orbitals for every body type and without `bt_`.

Run tests:
```bash
cd build && ctest -R T_
```
