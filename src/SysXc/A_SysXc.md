# SysXc Module — Exchange-Correlation Functionals (DFT, TDDFT)

## Overview

The **SysXc** module is the density-functional approximation of the
electron-electron interaction for Kohn–Sham DFT and real-time TDDFT. It
provides the local Kohn–Sham potential of the spin densities

```
v_KS,s = J[ρ] + v_xc,s[ρ↑, ρ↓],    ρ = ρ↑ + ρ↓,
```

with the Hartree potential J of SysInteraction (one Poisson solve of the total
density) and the exchange-correlation potential v_xc = δE_xc/δρ_s of a
functional of **Libxc** (LDA, GGA and their global hybrids), plus the fraction
`SysXc_exactExchangeFraction` of exact exchange of a hybrid, which the
consumers apply with the exchange of SysInteraction:

| Consumer | Operator |
|----------|----------|
| `method.mb.orbBased.tdks` (TDKS, real-time TDDFT; `src/Method/A_Method.md`) | i ∂φ_i/∂t = (h + v_KS,s − a K_s) φ_i |
| `groundSolver.scf.stdImpl` with a `sysXc` block (Kohn–Sham SCF; `src/GroundSolver/A_GroundSolver.md`) | F̂_KS = T̂ + V̂_ext + v_KS − a K̂ |

The Kohn–Sham energy (`Method_GetEnergy` of TDKS) is

```
E = Σ_i ⟨g_i|h φ_i⟩ + E_H + E_xc − a E_x,   E_H = ½ ∫ρ J,   E_x = ½ Σ_s Σ_i ⟨g_i|K_s φ_i⟩,
```

g = φ R the orbitals weighted with the 1-RDM blocks.

## Architecture

```
M_SysXc.f90                 ← contract: xc points, maps, functional, drivers
S_SysXc.f90                 ← fabricate: nodal defaults, dispatch to the family
S_SysXc_KohnSham.f90        ← grid-agnostic drivers on the contract
├── Fem/                    ← grid.fem: nodes, FEM nodal gradient (GGA)
├── Ylm/                    ← grid.ylm: angular product quadrature, spherical gradient (GGA)
├── Linear/                 ← grid.linear: nodes (1D LDA)
└── Libxc/                  ← the functional (shared by all families)
M_Utils_XcLib (src/Utils/XcLib/) ← bind(C) layer of the Libxc C API
```

Two orthogonal parts behind one contract:

- the **grid family** (`sysXc.fem`, `sysXc.ylm`, `sysXc.linear`) binds the
  *xc points* and the maps between orbitals, densities and potentials on them;
- the **functional** (`sysXc.<family>.libxc`) binds the point-wise evaluation
  of the energy density and its partial derivatives.

The family fabricate dispatches to the functional with its JSON path, so one
functional implementation serves every grid. The Kohn–Sham drivers compose
both and know neither.

## Contract

### Xc points

Densities and local potentials are `real(R64)` arrays on the `SysXc_nPoints`
xc points with quadrature weights `SysXc_weights` (Σ_p w_p f_p ≈ ∫ f d^d r):

| Family | Xc points | Weights |
|--------|-----------|---------|
| Fem | interior nodes (= `Grid_nPoints`) | lumped-mass weights `Grid_Fem_weights` |
| Linear | nodes (= `Grid_nPoints`) | `Grid_Linear_weights` |
| Ylm | nRadial × nTheta × nPhi angular product quadrature | w_r 2π w_q / nPhi |

### Pointers bound by the family

| Pointer | Contract |
|---------|----------|
| `SysXc_Setup` | xc points, weights, tables (after `Grid_Setup`) |
| `SysXc_AddDensity(density, orbs, weightedOrbs)` | density += Re Σ_l conj(φ_l) g_l: the density of γ(r, r') = Σ_l g_l(r) φ_l*(r') |
| `SysXc_MultiplyWithPotential(resOrb, pot, orb)` | resOrb = v·orb, Hermitian in the grid metric: ⟨a|v b⟩ = Σ_p w_p conj(a_p) v_p b_p |
| `SysXc_FillInteractionSrc(src, density)` | source of a density in the layout of `SysInteraction_FillInteractionSrc` |
| `SysXc_AddInteractionPotential(pot, V)` | pot += Re V on the xc points for an interaction potential V |
| `SysXc_ApplyGradient(gradF, f)` | G f in orthonormal directions (GGA; null if the family has none) |
| `SysXc_AddGradientAdjoint(f, gradF)` | f += W⁻¹ Gᵀ W gradF, the adjoint in the quadrature metric |

The nodal families (Fem, Linear) use the defaults of `S_SysXc.f90`
(point-wise products, src = ρ w, V real part); Ylm overrides all maps.

### Pointer bound by the functional

`SysXc_EvaluateFunctional(energyDensity, vRho, vSigma, rho, sigma)` on point
arrays (rho(:, nSpin), sigma(:, 1 or 3)): energy density e = ρ ε_xc and its
partial derivatives ∂e/∂ρ_s, ∂e/∂σ_c.

### Grid-agnostic drivers (`S_SysXc_KohnSham.f90`)

| Procedure | Result |
|-----------|--------|
| `SysXc_FillXcPotential(xcPot, xcEnergy, density)` | v_xc,s and E_xc of densities (nSpin 1: total density, unpolarized functional; 2: ρ↑, ρ↓) |
| `SysXc_FillHartreePotential(hartreePot, hartreeEnergy, density, time)` | J and E_H from one SysInteraction solve (zero without `sysInteraction`) |
| `SysXc_FillKohnShamPotential(ksPot, eH, eXc, density, time)` | v_KS,s = J + v_xc,s |
| `SysXc_Integrate(f)` | Σ w f |

### Data

| Name | Meaning |
|------|---------|
| `SysXc_nPoints`, `SysXc_weights` | xc points and weights |
| `SysXc_dimension` | spatial dimension of the grid (functionals of another dimension are refused) |
| `SysXc_gradientQ` | the functional needs gradients (any GGA) |
| `SysXc_exactExchangeFraction` | a of a hybrid (0 otherwise) |
| `SysXc_hartreeQ` | a `sysInteraction` block exists (J included) |
| `SysXc_PARALLEL_MIN_WORK` | smallest loop that runs OpenMP-parallel |

## Discrete functional and its exact derivative

The energy is the quadrature of the energy density with σ from the gradient
operator G of the family,

```
E_xc = Σ_p w_p e(ρ_p, σ_p),    σ = |Gρ|²   (spin: σ↑↑, σ↑↓, σ↓↓),
```

and the potential is its **exact derivative** with respect to the density
values, v = (1/w_p) ∂E_xc/∂ρ_p:

```
v_s = ∂e/∂ρ_s + W⁻¹Gᵀ W (2 ∂e/∂σ_ss Gρ_s + ∂e/∂σ_ud Gρ_s').
```

This holds for **any** linear G, so the Kohn–Sham orbital equations are the
Hamiltonian flow of the discrete Kohn–Sham energy (TDDFT conserves it, the SCF
minimizes it), and the GGA machinery is grid-agnostic: a family provides G and
its metric adjoint, nothing else. Negative densities (round-off of weighted
orbitals) count as zero; below Libxc's density threshold the energy density
and potential are zero.

## Grid families

### Fem (`sysXc.fem`)

Nodal xc points. The gradient is the nodal FEM derivative of `Grid_Fem_ApplyGradient`:
per axis W_a⁻¹ B_a with the assembled weak 1D derivative B(i, j) = ∫ L_i L_j′
dx (GLL quadrature), i.e. the derivative of the element interpolants averaged
with their weights over the elements sharing a node, the lumped-mass L²
projection of the piecewise derivative. It is exact for polynomials of degree
nLocals − 1 per axis (continuous across elements), banded (the element
couplings of the stiffness band) and OpenMP-parallel over z;
`Grid_Fem_AddGradientAdjoint` applies W_a⁻¹ B_aᵀ.

### Ylm (`sysXc.ylm`)

Every radial point carries an angular product quadrature of its own (not the
node tables of `Grid_Ylm`): `nTheta` Gauss–Legendre nodes in cos θ × `nPhi`
uniform angles. Tables Θ_lm(x_q) (Condon–Shortley, the recurrences of
Grid_Ylm), 2π w_q Θ_lm, dΘ_lm/dθ (ladder relations) and 2π w_q dΘ_lm/dθ up to
the largest l and |m| needed. One real DGEMM per m maps a channel matrix
(nRad × nL) to the nodes (nRad × nTheta) and back.

| Key | Default | Meaning |
|-----|---------|---------|
| `sysXc.ylm.nTheta` | 2 lmax + 3 | Gauss–Legendre nodes (≥ lmax + 1) |
| `sysXc.ylm.nPhi` | 1 if mmax = 0 or `method.mb.orbBased.mOfOrbs`, else 4 mmax + 1 | φ points; 1 selects the axial mode |

- **General mode** (nPhi > 1): orbitals are synthesized at every φ_k from all
  their m. Densities have |m| ≤ 2 mmax, so nPhi ≥ 4 mmax + 1 resolves them
  exactly.
- **Axial mode** (nPhi = 1): every m block of an orbital is transformed and
  multiplied separately (no m coupling); the density is the φ-average
  Σ_m |f_m|², exact for orbitals of definite m (atoms in linearly polarized
  fields). Without `mOfOrbs` and mmax > 0 it is the axial approximation (a
  note is printed).
- **Hartree potential**: the source is the channel projection of the xc
  density onto the layout of `sysInteraction.ylm` (radial metric weights
  included) and the potential channels are evaluated back on the xc points,
  one radial Poisson solve per channel and no Gaunt products. For
  nTheta ≥ 2 lmax + 1 this equals the Gaunt-product path of SysInteraction to
  round-off (`T_SysXc_Ylm`).
- **Gradient** (GGA, requires `grid.ylm.fedvr`): spherical frame,
  ∂_r f = (D(r f) − f)/r with the banded nodal FEDVR derivative D (r f vanishes
  at the excluded node r = 0), (1/r) ∂_θ f by projection onto the channels
  l ≤ 2 lmax of a density and differentiation of Θ_lm, (1/(r sin θ)) ∂_φ f by
  the discrete Fourier series in φ; the adjoints by transposition (analysis
  with dΘ/dθ, synthesis with Θ; −∂_φ). Exact for band-limited fields.

Refused: `grid.ylm.fedvrEcs` (complex-scaled radial grids have no real
density).

### Linear (`sysXc.linear`)

Nodal xc points of `grid.linear` (const or fedvr), dimension 1: only the 1D
LDA functionals of Libxc qualify (e.g. `lda_x_1d_soft` with `beta`,
`lda_c_1d_csc` with `interaction` and `beta`, for the soft-Coulomb interaction
of `sysInteraction.linear.softYukawa`). No gradients (Libxc has no 1D GGA).

## The functional: Libxc (`sysXc.<family>.libxc`)

```json
"sysXc": { "fem": { "libxc": { "functionals": ["gga_x_pbe", "gga_c_pbe"] } } }
```

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `functionals` | string(:) | required | Libxc names (case-insensitive, optional prefix `xc_`), e.g. `["lda_x", "lda_c_vwn"]`, `["gga_x_pbe", "gga_c_pbe"]`, `["hyb_gga_xc_pbeh"]` (PBE0), `["hyb_gga_xc_b3lyp"]`; may be empty |
| `coefficients` | real(:) | 1 each | weights c_i of the functionals |
| `densityThreshold` | real | Libxc's | densities below it count as zero |
| `extParams.<name>` | real(:) | Libxc's | external parameters of functional `<name>` |
| `exactExchangeFraction` | real | Σ_i c_i a_i | fraction of exact exchange (overrides that of the hybrids) |

An empty list with `exactExchangeFraction: 1` is pure exact exchange: the
Kohn–Sham energy is then the Hartree–Fock energy (T_He3d_03 reproduces the HF
limit of He).

Refused at fabricate: meta-GGAs (τ or Laplacian), range-separated hybrids
(CAM, LC: would need an erf-attenuated interaction), VV10 nonlocal
correlation, kinetic functionals, functionals without energy or potential, and
functionals of another spatial dimension than the grid.

Evaluation: chunks of `SysXc_Libxc_CHUNK_SIZE` = 2048 points, OpenMP-parallel
over the chunks with per-thread buffers in Libxc's interleaved spin layout
(Libxc only reads the functional handles, so the calls are reentrant). Each
functional is initialized unpolarized and polarized at fabricate.

### Libxc binding (`M_Utils_XcLib`)

`T_XcLib_Functional` holds the opaque `xc_func_type *` of `xc_func_alloc`, so
the binding does not depend on Libxc's struct layout and needs no Fortran
modules of Libxc (any Fortran compiler with a C build of the library):
`XcLib_FunctionalId`, `XcLib_Init`, `XcLib_End`, `XcLib_SetDensityThreshold`,
`XcLib_SetExtParams`, `XcLib_FillExtParamsDefault`,
`XcLib_ExactExchangeFraction`, `XcLib_EvaluateLda`, `XcLib_EvaluateGga`. The
CMake option `CODY_USE_LIBXC` (default ON) selects `S_Utils_XcLib_Libxc`; with
OFF the stub `S_Utils_XcLib_Stub` stops at the first use and the tests marked
`!> ctest-requires: libxc` are not registered.

## Usage

```fortran
call Grid_Fabricate
call SysKinetic_Fabricate
call SysPotential_Fabricate
call SysInteraction_Fabricate
call SysXc_Fabricate            ! after SysInteraction (Hartree layout)
call Method_Fabricate           ! method.mb.orbBased.tdks
...
call Grid_Setup
call SysKinetic_Setup
call SysPotential_Setup
call SysInteraction_Setup
call SysXc_Setup                ! after Grid_Setup, before Method_Setup
...
```

## Performance notes

- Loops around BLAS calls (the per-m GEMMs of Ylm) never run inside OpenMP
  regions of SysXc (the GEMMs thread themselves); element-wise loops are
  parallel only above `SysXc_PARALLEL_MIN_WORK` points. On a loaded 4-core
  machine nested OpenMP/BLAS regions around small GEMMs cost 50× (6.6 ms
  instead of 0.13 ms per Kohn–Sham application on Ne, lmax 1).
- Ylm: the Hartree potential costs one Poisson solve per channel from the xc
  density instead of one Gaunt product per orbital (O(lmax³ nRad) each).
- Kohn–Sham application per orbital: Fem one point-wise pass; Ylm 2(2 mmax + 1)
  GEMMs of (2 nRad × nL × nTheta) (general mode) or per m block (axial).
- The cost of a TDKS derivative is one Poisson solve (J), one functional
  evaluation and one application per orbital; a hybrid adds n_s(n_s + 1)/2
  pair potentials per spin.
- Real-time GGA: the gradient makes the kernel grow like k², and the split
  step of TDKS is stable only if every propagated h0 mode has E Δt < π. On
  Ylm, use `method.mb.orbBased.exactFieldFreeCut` with Ec Δt < π; on
  grid.fem, use a ~10× smaller time step than LDA (H2 PBE: 0.005). LDA has
  no such limit. Measurements and the mechanism: A_Method.md, "Time step
  with a GGA" in the TDKS section.

## Tests

| Test | Checks |
|------|--------|
| `test/componentTests/Utils/T_Utils_XcLib` | names, Slater exchange in closed form (unpolarized, polarized, ζ = 1), PBE at s = 0, v_σ by finite differences, PBE0/B3LYP fractions, external parameters |
| `test/componentTests/SysXc/Fem/T_SysXc_Fem` | density, Hartree energy 2√(2/π) (5e-7), Slater exchange in closed form (6e-7), gradient exact for polynomials (6e-12), adjoint identity (5e-14), PBE potential = derivative of E_xc (3e-10, polarized 5e-10), Kohn–Sham potential = J + v_xc, hybrid fractions |
| `test/componentTests/SysXc/Ylm/T_SysXc_Ylm` | general mode with mixed l, m: density integral, Hartree energy and J·φ equal the Gaunt path (8e-16), Hermiticity, angular gradient exact (1e-13), radial gradient exact for polynomials (6e-14), adjoint identity, PBE derivative (1e-10), axial = general mode for an axial density |
| `test/simulationTests/He3d/T_He3d_01` | He LDA vs NIST: E = −2.8348356 (−2.834836), ε = −0.5704248 (−0.570425) |
| `test/simulationTests/He3d/T_He3d_02` | He PBE: −2.8929349 (−2.892935) |
| `test/simulationTests/He3d/T_He3d_03` | pure exact exchange = HF limit −2.8616800, ε −0.9179556 |
| `test/simulationTests/He3d/T_He3d_04` | real-time TDKS with PBE at Δt 0.05 and the spectral cut Ec = 40 of the exact part: energy conserved to 5e-7 over 8 a.u., no weight above the cut (Ec 120 > π/Δt: 2e-2) |
| `test/simulationTests/Ne3d/T_Ne3d_18` | Ne LDA vs NIST: −128.2334813 (−128.233481), ε 1s, 2s, 2p to 2e-6 |
| `test/simulationTests/FemDft/T_FemDft_01` | He LDA and PBE on grid.fem: mesh error −1.8e-5; E_PBE − E_LDA equals the radial value to 8e-7 (independent gradient implementations) |
| `test/simulationTests/FemDft/T_FemDft_02` | H2 LDA and PBE0 (ACE and exact exchange agree to 1e-13) |
| `test/simulationTests/FemDft/T_FemDft_03` | H2 TDKS kick: TRK, Kohn–Sham energy conserved to 2e-11 with the frozen Kohn–Sham reference |
| `test/simulationTests/H3d/T_H3d_02` | H atom LSDA (one body type: polarized functional, ρ↓ = 0) by imaginary-time TDKS: −0.4786707 (NIST LSD −0.478671; an independent radial solver with the same Libxc functionals: −0.4786708) |
| `test/simulationTests/He1d/T_He1d_11` | TDKS with exact exchange only = TDHX (1D He, restricted): −2.4819265615 to 1e-10 |
| `test/simulationTests/He1d/T_He1d_12`, `T_He1d_13` | 1D LDA (`lda_x_1d_soft`, `lda_c_1d_csc`, soft-Coulomb): Kohn–Sham SCF −2.2013755814; imaginary-time TDKS reaches it (O(Δτ²) bias 2.6e-6 at Δτ = 0.01; Richardson to 1e-9) |

Apps: `app/molecules/{h2,h2o}/groundstate/dft/P_GroundKsFem` (Kohn–Sham
ground state) and `dynamics/tddft/P_KickTddftFem` (kick spectrum by
real-time TDDFT), see `app/README.md`.

## Adding a functional source or a grid family

- **Grid family**: a folder `src/SysXc/<Family>/` whose fabricate sets
  `SysXc_dimension` and `SysXc_nPoints`, binds `SysXc_Setup` (weights) and
  the maps that differ from the nodal defaults, optionally the gradient pair,
  and dispatches to `SysXc_Libxc_Fabricate("sysXc.<family>.libxc")`; register
  it in `S_SysXc.f90`.
- **Functional source**: a folder `src/SysXc/<Source>/` that binds
  `SysXc_EvaluateFunctional` and sets `SysXc_gradientQ` and
  `SysXc_exactExchangeFraction`; the families dispatch to it like to Libxc.
