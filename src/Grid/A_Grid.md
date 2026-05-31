# Grid Module – Agent Onboarding Guide

## Purpose

The **Grid** module provides the spatial discretization layer for quantum many-body simulations in CodyFortranRDM. It abstracts away coordinate-system details so that higher-level code (orbitals, Hamiltonians, propagators) can operate on generic field vectors without knowing whether the underlying representation is a 1D uniform mesh, a 3D spherical-harmonic expansion, or a discrete Hubbard lattice.

Key responsibilities:

| Capability | Description |
|------------|-------------|
| **Inner product** | Metric-aware ⟨f\|g⟩ respecting quadrature weights / volume elements |
| **Orthonormalization** | Gram–Schmidt with the correct metric |
| **Subspace projection** | Remove components along an orthonormal basis (gauge constraints) |
| **Coordinate access** | Back-end–specific arrays (xCoord, rCoord, weights, …) |

---

## Architecture Overview

The Grid module follows the **Interface Module + Implementation Submodule** pattern:

```
src/Grid/
├── M_Grid.f90                 # Public API: procedure pointers, Grid_nPoints
├── S_Grid.f90                 # Fabrication dispatcher + shared Gram–Schmidt
├── AGENTS.md                  # This file
│
├── Linear/                    # 1D Cartesian grids
│   ├── M_Grid_Linear.f90      # Interface: xmin, xmax, xCoord, weights
│   ├── S_Grid_Linear.f90      # Dispatcher for Const / Fedvr
│   ├── Const/                 # Uniform spacing
│   │   ├── M_Grid_Linear_Const.f90
│   │   └── S_Grid_Linear_Const.f90
│   └── Fedvr/                 # Finite-Element DVR (high accuracy)
│       ├── M_Grid_Linear_Fedvr.f90
│       └── S_Grid_Linear_Fedvr.f90
│
├── Square/                    # 2D Cartesian (x,y)
│   ├── M_Grid_Square.f90
│   ├── S_Grid_Square.f90
│   └── Const/
│
├── Polar/                     # 2D polar (r,φ)
│   ├── M_Grid_Polar.f90
│   ├── S_Grid_Polar.f90
│   └── Const/
│
├── Spherical/                 # 3D spherical (r,θ,φ)
│   ├── M_Grid_Spherical.f90
│   ├── S_Grid_Spherical.f90
│   └── Const/
│
├── Ylm/                       # Radial × spherical harmonics Y_lm
│   ├── M_Grid_Ylm.f90         # lmax, radial arrays, (l,m) accessors
│   ├── S_Grid_Ylm.f90         # SpatialProduct via Gaunt coefficients
│   ├── Const/
│   ├── Fedvr/
│   └── FedvrEcs/              # FEDVR + Exterior Complex Scaling
│
└── Lattice/                   # 3D discrete lattice (Hubbard models)
    ├── M_Grid_Lattice.f90     # xSize, ySize, zSize, periodicQ, code(ix,iy,iz)
    └── S_Grid_Lattice.f90
```

### Key Concepts

| Component | Role |
|-----------|------|
| **M_Grid** | Public interface: `Grid_nPoints`, procedure pointers (`Grid_Setup`, `Grid_InnerProduct`, `Grid_Orthonormalize`, `Grid_ProjectOnSubspace`) |
| **Procedure Pointers** | Runtime polymorphism – concrete back-ends bind these during fabrication |
| **Grid_Fabricate** | Reads JSON, selects back-end, binds pointers |
| **Grid_Setup** | Allocates coordinate / weight arrays (called once after fabrication) |
| **Grid_InnerProduct** | Computes ⟨f\|g⟩ = Σ_i conj(f_i) g_i w_i (metric-dependent) |

---

## Data Flow

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                        JSON Configuration                                    │
│  { "grid": { "linear": { "xmin": -10, "xmax": 10, "const": { ... } } } }    │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                          Grid_Fabricate                                      │
│  1. Bind shared routines (Orthonormalize, ProjectOnSubspace)                │
│  2. Probe JSON keys: grid.linear → grid.square → ... → grid.lattice         │
│  3. Dispatch to family Fabricate (e.g., Grid_Linear_Fabricate)              │
│     └─► Dispatch to variant Fabricate (e.g., Grid_Linear_Const_Fabricate)   │
│  4. Bind: Grid_Setup → Setup, Grid_InnerProduct → InnerProduct              │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                            Grid_Setup                                        │
│  1. Allocate weights(Grid_nPoints), xCoord(Grid_nPoints), ...               │
│  2. Fill coordinate and quadrature arrays                                   │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│               Client code (Orbs, SysKinetic, SysPotential, ...)             │
│  - Call Grid_InnerProduct(fConjg, f) to compute overlaps                    │
│  - Call Grid_Orthonormalize(orbSet) to enforce ⟨φ_i|φ_j⟩ = δ_ij            │
│  - Call Grid_ProjectOnSubspace(dOrbs, orbs) to remove gauge components      │
│  - Access coordinates/weights via M_Grid_Linear, M_Grid_Ylm, etc.           │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## JSON Configuration Reference

### 1D Linear Grid (Constant Spacing)

```json
{
  "grid": {
    "linear": {
      "xmin": -20.0,
      "xmax":  20.0,
      "const": {
        "nPoints": 256
      }
    }
  }
}
```

| Path | Type | Default | Description |
|------|------|---------|-------------|
| `grid.linear.xmin` | real | −10.0 | Left domain boundary |
| `grid.linear.xmax` | real | +10.0 | Right domain boundary |
| `grid.linear.const.nPoints` | int | 100 | Number of grid points |

### 1D Linear Grid (FEDVR)

```json
{
  "grid": {
    "linear": {
      "xmin": -20.0,
      "xmax":  20.0,
      "fedvr": {
        "nElements": 10,
        "order": 7
      }
    }
  }
}
```

| Path | Type | Default | Description |
|------|------|---------|-------------|
| `grid.linear.fedvr.nElements` | int | — | Number of finite elements |
| `grid.linear.fedvr.order` | int | — | Polynomial order per element |

### Ylm Grid (Radial + Spherical Harmonics)

```json
{
  "grid": {
    "ylm": {
      "lmax": 3,
      "rmax": 50.0,
      "fedvr": {
        "nElements": 20,
        "order": 9
      }
    }
  }
}
```

| Path | Type | Default | Description |
|------|------|---------|-------------|
| `grid.ylm.lmax` | int | 1 | Maximum angular momentum |
| `grid.ylm.rmax` | real | 20.0 | Radial box size |

### 3D Lattice Grid

```json
{
  "grid": {
    "lattice": {
      "xSize": 4,
      "ySize": 4,
      "zSize": 1,
      "xPeriodicQ": true,
      "yPeriodicQ": true,
      "zPeriodicQ": false
    }
  }
}
```

| Path | Type | Default | Description |
|------|------|---------|-------------|
| `grid.lattice.xSize` | int | 1 | Sites along x |
| `grid.lattice.ySize` | int | 1 | Sites along y |
| `grid.lattice.zSize` | int | 1 | Sites along z |
| `grid.lattice.{x,y,z}PeriodicQ` | bool | false | Periodic BCs |

---

## Mathematical Background

### Inner Product

For continuous grids, the inner product approximates the L² integral via quadrature:

```
⟨f|g⟩ = ∫ f*(x) g(x) dV  ≈  Σ_i  f*_i  g_i  w_i
```

where `w_i` encodes the quadrature weight and any Jacobian factors (e.g., r² dr for radial grids).

For lattice grids (discrete Hilbert space), the inner product is simply the dot product with unit weights.

### Gram–Schmidt Orthonormalization

Given a set {|φ₁⟩, …, |φ_N⟩}, the classical Gram–Schmidt procedure is:

```
for i = 1 to N:
    for j = 1 to i−1:
        |φ_i⟩  ←  |φ_i⟩  −  ⟨φ_j|φ_i⟩ |φ_j⟩
    end
    |φ_i⟩  ←  |φ_i⟩ / √⟨φ_i|φ_i⟩
end
```

This ensures ⟨φ_i|φ_j⟩ = δ_ij with respect to the metric defined by `Grid_InnerProduct`.

### Subspace Projection (Gauge Constraint)

In MCTDH-style methods, orbital time derivatives must remain orthogonal to the occupied subspace. The projection removes unwanted components:

```
|ψ⟩  ←  |ψ⟩  −  Σ_i  ⟨φ_i|ψ⟩ |φ_i⟩
```

**Numerical caveat:** This is a delicate operation. Small errors accumulate over thousands of time steps and interact nonlinearly with the coefficient dynamics. Periodic re-orthonormalization may be necessary.

---

## Back-end–Specific Features

### Ylm Grid: Spatial Product and Gaunt Coefficients

The Ylm back-end represents functions as f(r,Ω) = Σ_{l,m} f_{lm}(r) Y_{lm}(Ω). The product of two such functions involves Gaunt coefficients:

```
Y_{l₁,m₁}(Ω) · Y_{l₂,m₂}(Ω)  =  Σ_{l,m}  G(l₁,m₁; l₂,m₂; l,m)  Y_{l,m}(Ω)
```

`Grid_Ylm_SpatialProduct` computes this expansion efficiently, which is essential for evaluating two-body interaction potentials (Coulomb, Hartree) in atomic simulations.

### Lattice Grid: Site Indexing

For a 3D lattice with dimensions (Lx, Ly, Lz), sites are linearized in row-major order (x fastest):

```
index = ix + (iy-1)*Lx + (iz-1)*Lx*Ly
```

The `Grid_Lattice_code(ix,iy,iz)` array provides a lookup table for this mapping.

---

## Typical Usage in Client Code

```fortran
! Initialization sequence (once at startup)
call Grid_Fabricate          ! Reads JSON, wires back-end
call Grid_Setup              ! Allocates coordinates, weights

! Access grid size
print *, "Grid has", Grid_nPoints, "points"

! Compute overlap of two orbitals
overlap = Grid_InnerProduct(orb1, orb2)

! Orthonormalize a set of orbitals
call Grid_Orthonormalize(Orbs_orbs)

! Project dOrbs onto complement of occupied subspace (gauge)
call Grid_ProjectOnSubspace(dOrbs, Orbs_orbs)
```

---

## Adding a New Grid Back-end

### Example: Adding a "Cylindrical" grid family

1. **Create directory**: `src/Grid/Cylindrical/`

2. **Interface module** (`M_Grid_Cylindrical.f90`):
   ```fortran
   module M_Grid_Cylindrical
     use M_Utils_Types
     implicit none

     interface
       module subroutine Grid_Cylindrical_Fabricate
       end subroutine
     end interface

     real(R64) :: Grid_Cylindrical_rmin, Grid_Cylindrical_rmax
     real(R64), allocatable :: Grid_Cylindrical_weights(:)
     ! ... coordinates ...
   end module
   ```

3. **Implementation submodule** (`S_Grid_Cylindrical.f90`):
   ```fortran
   submodule(M_Grid_Cylindrical) S_Grid_Cylindrical
   contains
     module subroutine Grid_Cylindrical_Fabricate
       use M_Utils_Json
       use M_Grid
       ! Read JSON, set Grid_nPoints, bind Grid_InnerProduct => InnerProduct, etc.
     end subroutine

     pure function InnerProduct(fConjg, f) result(res)
       ! Implement with cylindrical metric weights
     end function

     subroutine Setup
       ! Allocate and fill coordinate/weight arrays
     end subroutine
   end submodule
   ```

4. **Register in S_Grid.f90**:
   ```fortran
   else if (Json_GetExistence("grid.cylindrical")) then
     call Grid_Cylindrical_Fabricate
   ```

5. **Update CMakeLists.txt** to include new source files.

---

## Dependencies

| Depends On | Purpose |
|------------|---------|
| `M_Utils_Types` | R64, I32 type kinds |
| `M_Utils_Json` | Configuration parsing (`Json_Get`, `Json_GetExistence`) |
| `M_Utils_Say` | Logging (`Say_Fabricate`, `Say_Setup`) |
| `M_Utils_NoOpProcedures` | Default no-op for Grid_Setup |
| `M_Utils_SphericalHarmonics` | Gaunt coefficients (Ylm only) |

---

## Testing Considerations

- **Inner product accuracy**: Verify ⟨f\|f⟩ ≈ 1 for normalized Gaussians / known functions
- **Orthonormalization**: After `Grid_Orthonormalize`, check |⟨φ_i\|φ_j⟩ − δ_ij| < ε
- **Coordinate ranges**: Ensure xCoord spans [xmin, xmax] with correct spacing
- **Weight sums**: For uniform grids, Σ w_i ≈ domain volume
- **Back-end independence**: Same physics with different grids (regression tests)

---

## Common Pitfalls

1. **Forgetting to call Grid_Setup**: Coordinate/weight arrays remain unallocated.
2. **Metric mismatch**: Using raw `dot_product` instead of `Grid_InnerProduct` gives wrong norms on non-uniform grids.
3. **Ylm indexing confusion**: The flattened order is (l=0,m=0), (l=1,m=−1), (l=1,m=0), (l=1,m=+1), … — not sorted by m first.
4. **Lattice periodicity**: Operators like `SysKinetic` must respect the `PeriodicQ` flags; mismatches cause edge artifacts.
5. **Projection drift**: `Grid_ProjectOnSubspace` is sensitive to accumulated errors; consider periodic re-orthonormalization during long propagations.

---

## Performance Notes

- **Linear grids**: Inner product is O(N) with vectorized `sum()`.
- **Lattice grids**: Uses intrinsic `dot_product` (no weights needed), typically BLAS-optimized.
- **Ylm spatial product**: O(l_max⁴ × N_radial) due to Gaunt coefficient loops; precompute where possible.
- **FEDVR grids**: Sparse structure can be exploited by kinetic operators but inner products remain dense.
