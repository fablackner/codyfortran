# OrbsInit Module — Agent Onboarding Guide

## Overview

The **OrbsInit** module provides orbital initialization for quantum many-body simulations in CodyFortranRDM. It constructs initial single-particle wavefunctions (orbitals) on various spatial discretizations, serving as the starting point for time-dependent and ground-state calculations.

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
            └── Ylm/             ─── Spherical harmonics R(r)·Y_lm
                └── HydrogenLike/─── Hydrogenic radial functions R_nl(r)
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

### Load from File
```json
{
  "orbsInit": {
    "load": { }
  }
}
```
Reads binary files `orb{bt}_{ind}.in` via DataStorage routines.

### GridPoint (δ-basis)
```json
{
  "orbsInit": {
    "gridPoint": { }
  }
}
```
Produces `φ_i(j) = δ_{ij}` normalized via grid inner product.

## Initialization Sequence

```fortran
! 1. Fabricate (reads JSON, binds procedure pointers)
call OrbsInit_Fabricate

! 2. Setup (allocates buffers if needed — currently no-op for most backends)
call OrbsInit_Setup

! 3. Initialize orbitals
call OrbsInit_Initialize(orbs)

! Alternative: initialize single orbital
call OrbsInit_InitializeOrb(orb(:), orbital_index, body_type)
```

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

5. **Update CMakeLists.txt** to include new files.

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

Run tests:
```bash
cd build && ctest -R T_
```
