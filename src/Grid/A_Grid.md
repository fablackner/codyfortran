# Grid Module – Agent Onboarding Guide

## Purpose

The **Grid** module provides the spatial discretization layer for quantum many-body simulations in CodyFortran. It abstracts away coordinate-system details so that higher-level code (orbitals, Hamiltonians, propagators) can operate on generic field vectors without knowing whether the underlying representation is a 1D uniform mesh, a 3D spherical-harmonic expansion, or a discrete Hubbard lattice.

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
├── Prolate/                   # Prolate spheroidal (homonuclear diatomics)
│   ├── M_Grid_Prolate.f90     # R, mmax, (xi,eta) arrays, m-channel accessors
│   ├── S_Grid_Prolate.f90     # SpatialProduct via azimuthal m-convolution
│   └── Fedvr/                 # FEDVR in xi, Gauss-Legendre DVR in eta
│
├── Fem/                       # 3D Cartesian spectral elements (molecules)
│   ├── M_Grid_Fem.f90         # axes (nodes, weights, eigenvalues), centers, operator API
│   ├── S_Grid_Fem.f90         # mesh: fabrication and grading, Grid_Setup, inner product, node layout
│   ├── S_Grid_Fem_Laplacian.f90  # 1D operator tables (private), Laplacian, FDM transforms, shifted solves
│   └── S_Grid_Fem_Poisson.f90 # child of S_Grid_Fem_Laplacian: multipole tables, Poisson solves
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
| `grid.ylm.mmax` | int | `lmax` | Maximum \|m\| of orbitals (clamped to `lmax`) |
| `grid.ylm.rmax` | real | 20.0 | Radial box size |

**m truncation.** Orbitals store the channels 0 ≤ l ≤ lmax, |m| ≤ min(l, mmax),
l-major with m ascending and radial points innermost
(`Grid_nPoints = nRadial · Grid_Ylm_NChannels(lmax, mmax)`; the channel
layout and count delegate to `M_Utils_SphericalHarmonics`). Lowering `mmax`
is useful when the dynamics conserves the orbital m quantum numbers, e.g. a
linearly z-polarized laser acting on an atom whose initial orbitals have
|m| ≤ mmax. Potential-like fields carry their own caps, just like their own
lmax: `SysInteraction_Ylm_mmax` (`sysInteraction.ylm.mmax`, default
min(2·mmax, lmaxPot), exact since a product of two orbitals has |m| ≤ 2·mmax)
and `SysPotential_Ylm_mmax` (set by the model; models with mmax > 0 refuse to
run when `grid.ylm.mmax < lmax`).

`Grid_Ylm_{Get,Set,Add}LmComponent` take `mmax_`; `Grid_Ylm_SpatialProduct`
takes `lMaxOut_`, `lMaxConjg_`, `lMax_`, `mMaxOut_`, `mMaxConjg_`, `mMax_`.
All default to the orbital layout (`Grid_Ylm_lmax`, `Grid_Ylm_mmax`); pass only
the caps of non-orbital operands. The accessors do not range-check: callers
must request |m| ≤ min(l, cap), which is guaranteed by validating inputs at
fabrication (e.g. OrbsInit checks every |m| ≤ mmax).

**Node representation of fixed-m fields.** A field whose channels all have
the same m is Σ_l f_lm(r) Θ_lm(θ) e^{imφ} with real Θ_lm (Y_lm = Θ_lm e^{imφ},
Condon–Shortley phase). `Grid_Ylm_ToNodes` evaluates Σ_l f_lm(r) Θ_lm(x_q) on
the Gauss–Legendre nodes x_q in cos θ and leaves the phase implicit: the
pointwise product of node fields with m₁ and m₂ represents m₁ + m₂, the complex
conjugate −m, and a product is integrated over φ (factor 2π, included in the
node weights) only when its m add up to zero. Such products are polynomials in
cos θ (the sin^|m| θ factors pair up), so the quadrature is exact for total
degree ≤ 4 lmax + 5, as for m = 0. The table sizes are set by `Grid_Ylm_Fabricate`, the
tables in `Grid_Ylm` Setup (`SetupQuadrature`, normalized associated Legendre
recurrences) and are used by
`SpatialProductM0` and by Method's `ApplyCorrelationOpYlmNodes` and
`FillH2YlmNodes` (`S_Method_Mb_OrbBased.f90`, orbitals of declared m,
`method.mb.orbBased.mOfOrbs`); see also the Ylm spatial
product bullet under Performance Notes.

**Parity split.** The nodes are symmetrized exactly in Setup, x_{n+1−q} = −x_q
(odd n, middle node x = 0), and Θ_lm(−x) = (−1)^{l+m} Θ_lm(x).
`Grid_Ylm_ChannelsToNodes` therefore transforms the degrees with even and
odd l − |m| (every second column of fM) separately on the nodes x_q ≥ 0 and
restores the mirrored half by a butterfly (f = e ± o);
`Grid_Ylm_NodesToChannels` folds the node field into its even and odd parts
first (a local, i.e. thread-private, buffer). Both halve the GEMM work of the
unsplit tables, which are kept for `SpatialProductM0`; the result differs
from the unsplit transform by round-off only
(`test/componentTests/Grid/Ylm/T_Grid_Ylm_Nodes.f90`). The symmetrization
also changes the unsplit tables and `Grid_Ylm_nodeWeights` by round-off, so
`SpatialProductM0` (every Ylm run, with or without
`method.mb.orbBased.mOfOrbs`) is reproduced to round-off, not bitwise.

| Name | Kind | Contract |
|------|------|----------|
| `Grid_Ylm_ToNodes(fNodes, f, m)` | subroutine | fNodes(r, q) = Σ_l f_lm(r) Θ_lm(x_q), \|m\| ≤ l ≤ `Grid_Ylm_lmax`, for f in the orbital layout; other m are ignored; requires \|m\| ≤ `Grid_Ylm_mmax` |
| `Grid_Ylm_FromNodes(f, fNodes, m)` | subroutine | Projection onto Y_lm, \|m\| ≤ l ≤ `Grid_Ylm_lmax`, into the orbital layout, all other channels zero; exact for angular degree ≤ 2 nNodes − 1 − lmax |
| `Grid_Ylm_ChannelsToNodes(fNodes, fM, m)` | subroutine | as `Grid_Ylm_ToNodes` on the channel matrix fM(r, l − \|m\| + 1), \|m\| ≤ l ≤ \|m\| + size(fM, 2) − 1 (real GEMMs on the complex arrays viewed as real, split by parity, see below) |
| `Grid_Ylm_NodesToChannels(fM, fNodes, m)` | subroutine | as `Grid_Ylm_FromNodes` into the channel matrix fM, without metric weights |
| `Grid_Ylm_ProjectOnM(f, m)` | pure subroutine | zeroes all channels of f (orbital layout) with an m other than m (used by `Orbs_ProjectOnSymmetry`) |
| `Grid_Ylm_nNodes` | integer | 2 lmax + 3 Gauss–Legendre nodes in cos θ |
| `Grid_Ylm_lmaxNodes` | integer | 2 lmax + 2, largest tabulated l |
| `Grid_Ylm_mmaxNodes` | integer | min(2 mmax, 2 lmax + 2), largest tabulated \|m\| (pair densities and potentials of orbitals with \|m\| ≤ mmax) |
| `Grid_Ylm_nodeWeights` | real(:) | 2π w_q, Σ_q = 4π |

### Prolate Grid (Homonuclear Diatomic Molecules)

Two-center prolate spheroidal coordinates (ξ, η, φ) with the nuclei at the
foci z = ±a, a = R/2. The volume element a³(ξ²−η²) dξ dη dφ cancels both
nuclear Coulomb singularities analytically. Fields are stored in azimuthal
channels ψ = Σ_m f_m(ξ,η) e^{imφ}/√(2π) (m is a good quantum number of the
linear molecule); channels are laid out contiguously in the interleaved order
m = 0, −1, +1, −2, +2, … Discretization: FEDVR in ξ (ξ = 1 is a grid point,
ξmax carries a Dirichlet condition) and a single-element Gauss–Legendre DVR
in η (interior nodes only).

```json
{
  "grid": {
    "prolate": {
      "R": 2.0,
      "ximax": 12.0,
      "mmax": 0,
      "nEta": 12,
      "fedvr": {
        "nElements": 8,
        "nLocals": 10
      }
    }
  }
}
```

| Path | Type | Default | Description |
|------|------|---------|-------------|
| `grid.prolate.R` | real | 2.0 | Internuclear distance (foci at ±R/2) |
| `grid.prolate.ximax` | real | 10.0 | Box boundary in ξ (z-extent is R/2·ξmax) |
| `grid.prolate.mmax` | int | 0 | Maximum azimuthal quantum number (σ: 0, π: 1, …) |
| `grid.prolate.nEta` | int | 8 | Gauss–Legendre points in η |
| `grid.prolate.fedvr.nElements` | int | 10 | Finite elements in ξ |
| `grid.prolate.fedvr.nLocals` | int | 7 | Lobatto points per element |

Companion back-ends: `sysKinetic.prolate.laplacian.fedvr`,
`sysPotential.prolate.coulomb` (charge per nucleus, −2Zξ/(a(ξ²−η²))),
`sysInteraction.prolate.coulomb.stdImpl` (Neumann-expansion Poisson solve per
(τ,m) channel), `orbsInit.prolate.lcao` (gerade/ungerade hydrogen-like LCAO).
Validated against the exact H₂⁺ 1sσg energy at R = 2 and the H₂ Hartree–Fock
limit at R = 1.4 (`test/simulationTests/H2/`). Note: odd-|m| channels carry a
non-polynomial factor ((ξ²−1)(1−η²))^{|m|/2} and converge only algebraically;
even-m channels are fully spectral.

### FEM Grid (3D Cartesian Spectral Elements for Molecules)

A finite-element grid for molecules (H₂, H₂O, …) in the box
[xmin, xmax] × [ymin, ymax] × [zmin, zmax]. The mesh is the tensor product of
three 1D element partitions graded toward the *centers* (nuclei); every center
coordinate is an element boundary on every axis, so every nucleus sits at an
element vertex where the C⁰ basis captures its cusp. Basis: Lagrange
polynomials on `nLocals` Gauss–Lobatto–Legendre (GLL) nodes per element edge
(degree p = nLocals − 1), continuous through the shared element-boundary
nodes. GLL quadrature makes the mass matrix diagonal, W = Wx ⊗ Wy ⊗ Wz (lumped
mass, as in FEDVR). Orbitals obey homogeneous Dirichlet conditions on the box
surface, so the two box nodes of each axis are removed:
n_a = nElements_a · (nLocals − 1) − 1 interior nodes per axis.

```json
{
  "grid": {
    "fem": {
      "nLocals": 6,
      "xmin": -10.0, "xmax": 10.0,
      "ymin": -10.0, "ymax": 10.0,
      "zmin": -10.0, "zmax": 10.0,
      "centers": {
        "x": [0.0, 1.4305227, -1.4305227],
        "y": [0.0, 0.0, 0.0],
        "z": [0.0, -1.1071202, -1.1071202],
        "hMin": [0.05, 0.2, 0.2]
      },
      "growth": 1.8,
      "hMax": 2.0,
      "multipoleLmax": 6
    }
  }
}
```

| Path | Type | Default | Description |
|------|------|---------|-------------|
| `grid.fem.nLocals` | int | 6 | GLL nodes per element edge (polynomial degree nLocals − 1), ≥ 2 |
| `grid.fem.xmin` / `xmax` (`y`, `z` alike) | real | −10.0 / 10.0 | Box bounds (Dirichlet surface), min < max |
| `grid.fem.centers.x` / `.y` / `.z` | real(:) | none | Center (nucleus) coordinates, parallel arrays of equal size, each at least 0.5·hMin (its own) inside every box face (error stop otherwise) |
| `grid.fem.centers.hMin` | real(:) | 0.2 each | Smallest element size at each center, > 0 |
| `grid.fem.growth` | real | 1.5 | Geometric growth of the element size away from the centers, ≥ 1 |
| `grid.fem.hMax` | real | 2.0 | Largest element size, > 0 |
| `grid.fem.xBounds` / `yBounds` / `zBounds` | real(:) | graded | Explicit element boundaries of one axis (override the grading): strictly increasing, first and last on the box bounds (within 1e-10 of the box length, then snapped; NaN-safe checks); a center coordinate missing from them is reported (that nucleus is not at a vertex) |
| `grid.fem.multipoleLmax` | int | 6 | Largest l of the multipole boundary values of `Grid_Fem_SolvePoisson`, ≥ 0 |
| `grid.fem.multipoleCenter` | real(3) | centroid of the centers (box center without centers) | Expansion center of those boundary values, strictly inside the box; put it at the charge (see *convergence* below) |

JSON supports only 1D arrays, hence the parallel center arrays.

**Mesh grading** (per axis, in `Grid_Fem_Fabricate`, so `Grid_nPoints` is set
after fabrication). Target element size
s(x) = min(hMax, min_k(h_k + (growth − 1)|x − c_k|)) over the center
coordinates c_k of the axis. Center coordinates closer than
0.25·min(h_k, h_l) are merged into one (in ascending order, each against the
kept coordinate of its cluster): the cluster keeps the center with the
smallest hMin (on ties the lower index), and s(x) runs over the kept centers
only, so a near duplicate (e.g. two nuclei 1e-7 apart from rounded geometry
files) gives exactly the mesh of one center instead of a sliver element whose
1/δ² stiffness destroys the FDM (δ = 1e-7 used to give a negative eigenvalue
and a Poisson error of 0.33). A note names every merged center off the kept
coordinate: it is not at an element vertex, but its Galerkin point charge
stays exact through the basis values. Mandatory boundaries are the box bounds
and the kept coordinates; a gap [p, q] between them gets n = max(1, ⌈F − 1e-9⌉)
elements, F = ∫_p^q dx/s(x) (trapezoidal rule on 4000 subintervals), with the
interior boundaries where the cumulative integral reaches jF/n (linear inverse
interpolation). Guards against degenerate meshes: every center must lie at
least 0.5·hMin inside the box faces, `Grid_Setup` warns if an element is
smaller than 0.1·min(hMin) (possible with explicit bounds) and stops if the
smallest FDM eigenvalue of an axis is not positive, and both Poisson solves
check that the spectrum D of −lap is positive. The 1D nodes and weights come
from `Fedvr_CreateCtxFromBounds` (`M_Utils_Fedvr`, explicit non-uniform
element boundaries; element boundaries are nodes exactly). `Grid_Setup`
builds the 1D operators (`Grid_Fem_SetupAxes`) and the 3D weights, prints a
mesh summary (elements, nodes, smallest/largest element per axis), then
builds the tables of the multipole boundary values (`Grid_Fem_SetupPoisson`)
and ends the summary with them (expansion center, convergence radius, memory
of the face-harmonics table).

**Layout.** Fields are nodal values f(ix, iy, iz), flattened with x fastest,
i = ix + nx(iy − 1) + nx·ny(iz − 1) (`Grid_Fem_Index`), and
`Grid_InnerProduct(f, g)` = Σ_i conj(f_i) g_i w_i with the product weights
w = wx(ix) wy(iy) wz(iz) (`Grid_Fem_weights`). The generic Gram–Schmidt and
subspace projection of `S_Grid.f90` apply unchanged.

| Name | Kind | Contract |
|------|------|----------|
| `T_Grid_Fem_Axis` | type | geometry and 1D spectrum of an axis: `nElements`, `nPoints` (interior nodes), `lowerBound`, `upperBound`, `bounds(0:nElements)`, `points`, `weights` (assembled GLL weights), `eigenvalues` (ascending, of W⁻¹K = −d²/dx² with Dirichlet conditions). The 1D operator tables (FDM eigenvectors, stiffness band, boundary columns) are private to `S_Grid_Fem_Laplacian` (see *Code* below) |
| `Grid_Fem_axes(3)` | type(T_Grid_Fem_Axis) | x, y, z axes |
| `Grid_Fem_nLocals`, `Grid_Fem_nx`, `Grid_Fem_ny`, `Grid_Fem_nz` | integer | nodes per element edge; interior nodes per axis (Grid_nPoints = nx·ny·nz) |
| `Grid_Fem_nCenters`, `Grid_Fem_centers(3, nCenters)`, `Grid_Fem_centerHMin(nCenters)` | integer, real | centers (nuclei) and their smallest element sizes |
| `Grid_Fem_weights(Grid_nPoints)` | real(:) | 3D product weights |
| `Grid_Fem_multipoleLmax`, `Grid_Fem_multipoleCenter(3)` | integer, real | multipole boundary values of the Poisson solve |
| `Grid_Fem_multipoleRadius` | real | convergence radius of those boundary values: distance from the expansion center to the nearest face node (set in `Grid_Setup`) |
| `Grid_Fem_SetupAxes`, `Grid_Fem_SetupPoisson` | subroutines | steps of `Grid_Setup`, called by it only (not API): 1D operators of the axes; multipole tables and `Grid_Fem_multipoleRadius` |
| `Grid_Fem_ApplyLaplacian(lapF, f)` | subroutine | lapF = −Σ_a (W_a⁻¹K_a) f along axis a (banded stencils) |
| `Grid_Fem_ToEigenbasis(c, f)` | subroutine | c = S3ᵀ W3 f, S3 = Sx ⊗ Sy ⊗ Sz |
| `Grid_Fem_FromEigenbasis(f, c)` | subroutine | f = S3 c |
| `Grid_Fem_SolveShiftedLaplacian(u, f, shift)` | generic subroutine | u = (−lap + shift)⁻¹ f = S3 (D + shift)⁻¹ S3ᵀ W3 f, D(i,j,k) = lx(i) + ly(j) + lz(k), for a `real(R64)` shift (specific `Grid_Fem_SolveShiftedLaplacianReal`; shift > −min D) or a `complex(R64)` shift (`Grid_Fem_SolveShiftedLaplacianComplex`, e.g. 2/(i dt) for Crank–Nicolson; the mode products stay real, only (D + shift)⁻¹ mixes real and imaginary coefficients; a zero imaginary part takes the real path) |
| `Grid_Fem_SolvePoisson(pot, load)` | subroutine | Galerkin solve of −lap pot = 4πρ for load = wρ (the weighted pair density of SysInteraction): K_II pot = 4π load − K_IB V_B, V_B the multipole expansion of the load (moments of the real and imaginary parts in one pass, V_B = T Q with the precomputed face table T); warns once per run if more than 1e-4 of Σ\|load\| lies outside the convergence radius |
| `Grid_Fem_SolvePoissonPointCharges(pot, charges, positions)` | subroutine | Galerkin potential of point charges q_I at R_I = positions(:, I) (strictly inside the box): load_i = Σ q_I φ_i(R_I), exact V_B = Σ q_I/\|r_B − R_I\| |
| `Grid_Fem_Index(ix, iy, iz)` | pure function | flattened node index |
| `Grid_Fem_FillNodeCoordinate(coord, direction)` | subroutine | coord(i) = direction · r_i at every node (`real(R64)`, allocatable intent(out), allocated to Grid_nPoints; after `Grid_Setup`), evaluated as d₁x + (d₂y + d₃z): an axis unit vector gives x, y or z exactly. Used for e_pol·r of `sysGauge.fem.length`, the coordinates of the FEM dipole printer and the kicks e^(iκ r_a) of the Fem tests and `app/fem_molecules` |

All field arguments are complex(R64), contiguous, of size `Grid_nPoints`
(the real node coordinates of `Grid_Fem_FillNodeCoordinate` excepted). The
routines keep no module-level scratch (local allocations only; the tables of
`Grid_Setup` are read only) and may be called inside OpenMP parallel regions;
their own OpenMP loops (over z) carry `if (.not. omp_in_parallel())`. The
transforms, shifted solves and `Grid_Fem_SolvePoisson` evaluate whether the
input has an imaginary part and extract both parts before they write the
output, so an in-place call (same array for input and output) also works;
it violates the Fortran aliasing rules, though, and the code base does not
rely on it (`Grid_Fem_ApplyLaplacian` needs distinct arrays).

**Code.** Three submodules of `M_Grid_Fem`, layered like their call graphs
(the mesh calls the other two only through the setup steps of `Grid_Setup`,
`Grid_Fem_SetupAxes` and `Grid_Fem_SetupPoisson`, and uses none of their
tables or kernels; the Poisson solves need only the shared kernels and
boundary columns of the Laplacian layer):

| File | Submodule | Contents |
|------|-----------|----------|
| `S_Grid_Fem.f90` | `submodule(M_Grid_Fem) S_Grid_Fem` | mesh: `Grid_Fem_Fabricate` (grading, explicit bounds), `Grid_Setup` (orchestration: `Grid_Fem_SetupAxes`, 3D weights, mesh summary, `Grid_Fem_SetupPoisson`), inner product, `Grid_Fem_Index`, `Grid_Fem_FillNodeCoordinate` |
| `S_Grid_Fem_Laplacian.f90` | `submodule(M_Grid_Fem) S_Grid_Fem_Laplacian` | `Grid_Fem_SetupAxes`; the private 1D operator tables `axisOperators(3)` (type `T_AxisOperators`: FDM eigenvectors, stiffness band, boundary columns); `Grid_Fem_ApplyLaplacian`, the eigenbasis transforms, the shifted solves and their real kernels (mode products, spectral inverse, real and imaginary parts, argument checks) |
| `S_Grid_Fem_Poisson.f90` | `submodule(M_Grid_Fem:S_Grid_Fem_Laplacian) S_Grid_Fem_Poisson` | `Grid_Fem_SetupPoisson` (multipole tables, convergence radius, the multipole summary lines), `Grid_Fem_SolvePoisson`, `Grid_Fem_SolvePoissonPointCharges` and their helpers (face nodes and boundary coupling, each one loop over the normal axis a with the tangential axes b < c of `TangentialAxes`; moments; point loads); its tables are private to it |

`S_Grid_Fem_Poisson` is a *nested submodule*, a child of
`S_Grid_Fem_Laplacian` (the only one in the code base). It reaches the
operator tables and the FDM kernels of its parent by host association, so
none of them is public. This is standard Fortran 2008; gfortran compiles it,
and CMake's Makefile generator builds the child after its parent from the
`.smod` dependency, whatever the file names. `Grid_Fem_SetupAxes` and
`Grid_Fem_SetupPoisson` are declared in `M_Grid_Fem` only so that
`Grid_Setup` in the sibling `S_Grid_Fem` can call them (the pattern of the
`Method_Mb_OrbBased_SetupFieldFree*` entry points).

**Numerics.**

- *Lumped-mass stiffness.* K_a(i, j) = ∫ L_i′ L_j′ dx is exact with GLL
  quadrature (degree 2p − 2), assembled from the reference element matrix
  (`M_Utils_DerivativeFedvr`) divided by the element size. For a tensor mesh
  with GLL quadrature K3 = Kx⊗Wy⊗Wz + Wx⊗Ky⊗Wz + Wx⊗Wy⊗Kz exactly, so
  lap = −W⁻¹K3 = −Σ_a W_a⁻¹K_a, self-adjoint in the W metric. Only nodes
  sharing an element couple: a node at local position k = mod(i, p) of its
  element couples to i − k … i − k + p, an element-boundary node to i − p … i + p.
- *Fast diagonalization (FDM).* Per axis A = W^(−1/2) K_II W^(−1/2) = QΛQᵀ
  (DSYEV), S = W^(−1/2) Q. Then (K3 + sW3)⁻¹ = S3 (D + s)⁻¹ S3ᵀ and
  g(−lap) f = S3 g(D) S3ᵀ W3 f for any function g. The shifted solves
  (exposed as `SysKinetic_SolveShifted`) precondition the LOBPCG eigensolver
  and the GMRES of the Cayley time step; the transforms
  `Grid_Fem_ToEigenbasis` / `Grid_Fem_FromEigenbasis` remain for spectral
  analysis and tests. grid.fem deliberately has no exact kinetic propagation
  and no spectral cut: the dynamics propagates h0 = T + V_static by the
  Cayley exact part (A_Method.md, "Exact Field-Free Part on grid.fem: Cayley
  Step"), and `exactFieldFreeCut` stops at fabricate. The mode products split
  complex fields into real and imaginary parts (the imaginary transform is
  skipped when that part is identically zero) and use real DGEMMs: one along
  x on the (nx, ny·nz) matrix, nz along y on the (nx, ny) slices (OpenMP over
  z), one along z on the (nx·ny, nz) matrix — 2N(nx + ny + nz) flops per real
  transform.
- *Poisson solve with free-space boundary values.* −lap V = 4πρ in weak form
  on the interior test functions, K_II V_I = 4π load_I − K_IB V_B, solved
  exactly by FDM. With the lumped mass only face nodes couple to interior
  rows (edges and corners do not): on the x faces
  (K_IB V_B)(ix, iy, iz) = Kx_full(ix, box node) wy(iy) wz(iz) V_B, nonzero for
  ix in the first/last element. V_B(r) = Σ_lm Q_lm S_lm(r − c)/\|r − c\|^(2l+1)
  with real regular solid harmonics S_lm (Racah normalization,
  Helgaker–Jørgensen–Olsen recursion; S₁₀ = z, S₁₁ = x, S₁,₋₁ = y), whose
  polynomial coefficients are tabulated once in `Grid_Fem_SetupPoisson`
  (`SolidHarmonics_FillCoefficients` of `M_Utils_SolidHarmonics`); the moments
  Q_lm = Σ_i load_i S_lm(r_i − c) come from the Cartesian moments
  Σ_i load_i x^a y^b z^c (three contractions with 1D power tables, O(L·N);
  real and imaginary parts in one pass, one DGEMM along x). `Grid_Setup`
  precomputes the face-node geometry and the irregular solid harmonics
  T_lm(r_B) = S_lm(r_B − c)/\|r_B − c\|^(2l+1) of all n_B face nodes (one real
  n_B × (L+1)² table, 8·n_B·(L+1)² bytes: 12.3 MB on the 376119-node test mesh,
  25 MB on the 1.1M-node H₂O mesh for L = 6, ≈ 3.5× that for L = 12; printed in
  the setup summary), so the boundary values are one DGEMV (DGEMM for both
  parts), V_B = T Q.
- *Convergence of the multipole boundary values (limitation).* The expansion
  converges at a face node only for charge closer to c than that node, i.e.
  for a load inside the sphere of radius R = `Grid_Fem_multipoleRadius`
  (distance from c to the nearest face node). Charge outside makes the series
  diverge, and raising L makes it worse: in an elongated box (x, y ±6,
  z ±20, two Gaussians at z = ±12, typical of ionisation set-ups) the error is
  2.1 for L = 6 and 86 for L = 12, against 1.2e-3 in a ±20 cube. Off-center
  expansion centers cost accuracy as well (x ∈ [−6, 14], charge at the
  origin: 1.6e-3 about the box center, 1.2e-6 about the charge), hence the
  default center at the centroid of the nuclei. `Grid_Fem_SolvePoisson`
  therefore computes the fraction of Σ\|load\| at nodes with \|r − c\| ≥ R and,
  above 1e-4 (`MULTIPOLE_WARNING_FRACTION`), prints a warning once per run
  (OpenMP-safe) recommending a box centered on the charge or a larger box.
  The boundary-value error is of the order of the fraction times the total
  load over R, so the exponential tails of bound densities are harmless:
  1s density in a box of half-width 10: 3.8e-7 (the error stays the
  discretization error, 1.9e-6, for L = 0 … 12); H₂ RHF density (±10,
  T_Fem_03): 4.9e-7; H₂O with the diffuse STO start guess (±8, T_Fem_04):
  1.55e-5. Loads that reach the faces (elongated boxes, ionized or
  off-center densities) give 1e-3 … 1. Use boxes of comparable extents
  centered on the charge (per-block expansions for elongated boxes are not
  implemented).
- *Galerkin point charges.* The load of a charge q at R is q φ_i(R) (tensor
  Lagrange basis of the element containing R; q δ_iJ if R is the node J) and
  the box values are exact, Σ q/\|r_B − R\|. With q = −Z this is the Galerkin
  nuclear potential; its value at the nucleus node is finite, roughly
  −30 … −50 Z/hMin for nLocals 5–6, depending on the mesh (measured: −27.6
  for nLocals 5, growth 1.8; −40.9 … −49.6 for nLocals 6, growth 1.8 … 1.6;
  about −85 for nLocals 8), and its moments up to degree p per axis are exact
  for charges off the nodes.

**Validation** (`test/componentTests/Grid/Fem/`, mesh: box ±8, nLocals 6,
centers (0, 0, ±0.7), hMin 0.2, growth 1.6, hMax 2 → n = (69, 69, 79),
N = 376119; identical to the numpy prototype):

| Check | Result |
|-------|--------|
| Σ w e^(−r²) − π^1.5 | −1.11e-8 |
| −½⟨g\|lap g⟩/⟨g\|g⟩, g = e^(−r²) (exact 1.5) | 1.4999993336 |
| ⟨f\|lap h⟩ vs ⟨lap f\|h⟩ (relative) | 3e-14 |
| (−lap + s) u − f, s = 0.7 / s = 3 − 40i (relative residual) | 5.4e-12 / 3.2e-13 |
| unit Gaussian charge: max \|V − erf(r)/r\| (r < 4) / self-energy error | 5.9e-7 / 4.3e-8 |
| Gaussian at (0.5, −0.3, 0.8), multipoleLmax 6: max \|V − erf(r)/r\| | 4.8e-6 (6.7e-3 with l = 0 only) |
| point charge −1 at (0, 0, 0.7): max \|V + 1/r\| for r > 1.5 / r > 3 | 6.7e-6 / 1.3e-6 |

`T_Grid_Fem_Asymmetric` tells the axes apart (a mesh symmetric in x and y
cannot catch x/y mix-ups): nLocals 5, box [−7, 8.5] × [−6, 6.5] × [−8, 7],
centers (0.4, −0.3, 0.9) and (−0.8, 0.5, −0.6), explicit y bounds,
multipoleLmax 12 about (0.3, −0.2, 0.4) → n = (59, 47, 63), all numbers of
the operators equal to the numpy prototype on the same mesh:

| Check | Result |
|-------|--------|
| −½⟨g\|lap g⟩/⟨g\|g⟩, Gaussian at (0.5, −0.4, 0.2) | 1.500000952185 |
| unit Gaussian charge there: max \|V − erf(r)/r\| / self-energy error | 2.12e-5 / −2.15e-7 |
| node-delta load: multipole (l ≤ 12) vs exact boundary values, max \|ΔV\| | 1.8e-9 (7.4e-9 for l ≤ 11, 5.9e-6 for l ≤ 6) |
| point charge −1 on a center node: r > 1.5 / r > 3 / value at the node | 3.7e-5 / 4.2e-6 / −108.7928651665 |
| operators inside an OpenMP parallel region vs serial calls (relative); a deliberately wide load (e^(−0.15 r²), fraction 1.8e-2 beyond R) raises the once-only multipole warning from inside the region | ≤ 5.2e-16 |
| `Grid_Fem_FillNodeCoordinate` along x, y, z (exact) / along d = (0.6, −1.3, 0.45) against d · r_i | 0 / 3.6e-15 |

`T_Grid_Fem_Sliver` checks the merge (two centers 1e-7 apart give exactly the
mesh of one center and accurate Poisson solves, also for the off-node point
charge), `T_Grid_Fem_CenterNearFace` the face-distance check.

End-to-end ground states (`test/simulationTests/Fem/`, meshes identical to the
numpy prototype, multipole boundary values about the origin with l ≤ 6):

| Test | Mesh | Result (prototype on the same mesh) |
|------|------|--------------------------------------|
| `T_Fem_01_HydrogenAtomLobpcg` | ±12, nLocals 6, hMin 0.4 → N = 205379 | E = −0.500012546111 (−0.500012546113; exact −0.5) |
| `T_Fem_02_H2plusLobpcg` | ±12 × ±12 × ±13, nLocals 6, hMin 0.2 → N = 447534 | E = −1.102636133119 (−1.102636133105; exact −1.1026342145) |
| `T_Fem_03_H2ScfRhf` | ±10, nLocals 5, hMin 0.3 → N = 130331 | E_el = −1.847953860878 (−1.847953860883; HF limit −1.8479152851) |
| `T_Fem_04_H2oScfRhf` | ±8, nLocals 4, hMin 0.1 (O) / 0.3 (H) → N = 72380 | E_tot = −76.403755960854 (−76.403755960801) |

Galerkin point nuclei are not variational: the coarse meshes lie below the
exact or HF-limit energies (H₂O by 0.34 Ha).

Real-time TDHF on the mesh of T_Fem_03 (restricted MCTDHX with nOrbs =
nBodies = [1, 1], orbital loaded from T_Fem_03, `propagator.splitStep.order2`
with the Cayley exact part of `method.mb.orbBased.exactFieldFreeQ`, RK4,
A_Method.md; T_Fem_10: MCTDHF on a coarser mesh):

| Test | Run | Result |
|------|-----|--------|
| `T_Fem_05_H2FieldFree` | dt = 0.01, 0.5 a.u. | ground state stationary: \|ΔE\| 2.3e-10, \|N − 2\| 3e-12, 1 − \|⟨φ0\|φ⟩\| 1.3e-12; 5 GMRES iterations per Cayley solve |
| `T_Fem_06_H2Kick` | kick e^(iκz), κ = 1e-3, dt = 0.02, 1 a.u. | d⟨z⟩/dt(0+) = Nκ (TRK) to −1.3e-6 on the grid, propagated slope to −1.5e-5 (second order in dt) |
| `T_Fem_07_H2KickDouble` | κ = 2e-3 | dipole response linear: \|ratio − 2\| 2.1e-8 |
| `T_Fem_08_ExponentialUnsupported` | `propagator.exponential` | refused at fabricate (no phi functions for the Cayley part) |
| `T_Fem_09_H2KickMeanFieldReference` | frozen mean-field reference (`exactFieldFreeMeanFieldQ`): ground state at dt = 0.1, kick at dt = 0.02 | ground state stationary (\|ΔE\| 6e-13); dipole equal to the plain split of T_Fem_06 to 2.7e-7 of the response; second order (p = 2.00) |
| `T_Fem_10_H2MctdhfMeanFieldReference` | MCTDHF, nOrbs = [2, 2], N = 15341, κ = 0.1, dt = 0.005 | the reference split follows the plain split to 1.6e-7 of the dipole response and 3.4e-7 Ha |

Companion back-ends (`sysKinetic.fem.laplacian` with the exact shifted
solve `SysKinetic_SolveShifted`, `sysPotential.fem.coulomb`,
`sysGauge.fem.length`, `sysInteraction.fem.coulomb`, `orbsInit.fem.sto`,
`absorber.fem.cosinus`) build on this API; see their module docs.

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
| `M_Utils_Fedvr`, `M_Utils_DerivativeFedvr` | GLL nodes/weights on (explicit) element bounds and reference derivative matrices (FEDVR back-ends, Fem) |
| `M_Utils_LapackLib`, `M_Utils_BlasLib` | FDM eigendecompositions (DSYEV) and mode products (DGEMM) (Fem) |
| `M_Utils_SolidHarmonics` | Polynomial coefficients of the real solid harmonics (Fem multipole boundary values) |
| `stdlib_sorting` | Sorting the center coordinates of the mesh grading (Fem) |

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
3. **Ylm indexing confusion**: The flattened order is (l=0,m=0), (l=1,m=−1), (l=1,m=0), (l=1,m=+1), … — not sorted by m first. The channel order is `SphericalHarmonics_ChannelIndex` (`M_Utils_SphericalHarmonics`).
4. **Lattice periodicity**: Operators like `SysKinetic` must respect the `PeriodicQ` flags; mismatches cause edge artifacts.
5. **Projection drift**: `Grid_ProjectOnSubspace` is sensitive to accumulated errors; consider periodic re-orthonormalization during long propagations.

---

## Performance Notes

- **Linear grids**: Inner product is O(N) with vectorized `sum()`.
- **Lattice grids**: Uses intrinsic `dot_product` (no weights needed), typically BLAS-optimized.
- **Ylm spatial product**: O(l_max⁴ × N_radial) due to Gaunt coefficient loops; precompute where possible.
  When all three operands are axially symmetric (`mMax*_ = 0`, e.g. `grid.ylm.mmax = 0`
  with a z-polarized laser) an exact Gauss–Legendre quadrature in cos θ is used instead:
  both factors are evaluated on the nodes, multiplied pointwise and projected back
  (O(l_max² × N_radial), matrix products). Tables cover l ≤ 2 lmax + 2 and products of
  total degree ≤ 4 lmax + 5; larger layouts fall back to the Gaunt sum.
  A spherically symmetric first factor (`lMaxConjg_ = 0`, e.g. a central external
  potential) takes precedence: Y₀₀ Y_lm = Y_lm/√(4π), so the product is an exact
  radial scaling of every channel of the second factor (O(N), no transforms).
  All three paths share one metric-weight step (`withWeightsQ_`).
- **Gaunt table** (`M_Utils_SphericalHarmonics`): SysInteraction and SysPotential
  tabulate only |m| ≤ max(grid mmax, operator mmax) in a compact (l, m) index;
  other m fall back to the direct Wigner-3j evaluation. With `mmax = 0` and
  lmax = 32 this is ~1 MB built in milliseconds instead of 2.4 GB and minutes.
- **FEDVR grids**: Sparse structure can be exploited by kinetic operators but inner products remain dense.
- **FEM grid**: on the 376119-node test mesh with 2 threads (Pentium G4560, OpenBLAS):
  `Grid_Fem_ApplyLaplacian` ≈ 10 ms (complex field), `Grid_Fem_ToEigenbasis` ≈ 24 ms
  (complex, two real transforms), `Grid_Fem_SolveShiftedLaplacian` ≈ 44 ms (complex),
  `Grid_Fem_SolvePoisson` ≈ 25 ms (real load; ≈ 52 ms complex),
  `Grid_Fem_SolvePoissonPointCharges` ≈ 21 ms. The transforms dominate:
  2N(nx + ny + nz) flops each, so the cost grows like N^(4/3). The multipole
  boundary values (precomputed face table, one DGEMV/DGEMM) add ≈ 1–3 ms to the
  two transforms of a real load (before: ≈ 11 ms on this mesh, ≈ 25 ms on the
  1.1M-node H₂O mesh, where a real solve now takes ≈ 94 ms).
