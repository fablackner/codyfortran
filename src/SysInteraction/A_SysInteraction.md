# SysInteraction Module — Agent Onboarding Guide

## Overview

The **SysInteraction** module implements particle-particle (two-body) interactions
for the CodyFortran quantum many-body simulation framework. It provides the
operator **Ŵ** that appears in the many-body Hamiltonian:

```
Ĥ = T̂ + V̂_ext + Ŵ
```

where `Ŵ = ½ Σᵢⱼ w(rᵢ,rⱼ)` describes electron-electron repulsion, Hubbard-U
interactions, or other pairwise potentials.

---

## Architecture

### Module/Submodule Pattern

This module follows the framework's **Interface Module + Implementation Submodule**
paradigm:

```
M_SysInteraction.f90          ← Public interface (abstract contracts)
S_SysInteraction.f90          ← Fabrication logic (runtime dispatch)
S_SysInteraction_MeanField.f90 ← mean field of a factored 1-RDM (generic, on the contract):
                                 FillDirectPotential, FillExchangeAce, AddMeanField
├── Lattice/
│   ├── M_SysInteraction_Lattice.f90
│   ├── S_SysInteraction_Lattice.f90
│   └── OnSite/...
├── Linear/
│   ├── M_SysInteraction_Linear.f90
│   ├── S_SysInteraction_Linear.f90
│   └── SoftYukawa/...
├── Ylm/
│   ├── M_SysInteraction_Ylm.f90
│   ├── S_SysInteraction_Ylm.f90
│   └── Coulomb/...
└── Fem/
    ├── M_SysInteraction_Fem.f90      ← nodal source, multiply, conjugate
    ├── S_SysInteraction_Fem.f90
    └── Coulomb/                      ← strength * Grid_Fem_SolvePoisson(src)
```

**Key rule:** Interface modules (`M_*.f90`) only import other `M_*` modules.
Implementation submodules (`S_*.f90`) contain the actual algorithms.

### Procedure Pointer Binding

At runtime, `SysInteraction_Fabricate` reads the JSON configuration and binds
procedure pointers to concrete implementations. This enables:

- Swapping algorithms without recompilation
- JSON-driven physics configuration
- Clean separation of interface from implementation

---

## Core API

### Public Procedure Pointers

| Pointer | Signature | Purpose |
|---------|-----------|---------|
| `SysInteraction_Setup` | `()` | One-time initialization (precompute kernels, FFT plans) |
| `SysInteraction_FillInteractionSrc` | `(src, orbConjg, orb)` | Build source ρ = ψ*·ψ |
| `SysInteraction_FillInteractionPotential` | `(pot, src, t, bt1, bt2)` | Solve for V from ρ |
| `SysInteraction_MultiplyWithInteractionPotential` | `(dOrb, pot, orb)` | Apply V·ψ |

### Module Procedures: Mean Field of a Factored 1-RDM

Three grid-agnostic helpers in `S_SysInteraction_MeanField.f90`, written only
against the contract above and `Grid_InnerProduct` (any grid), for the
one-body density matrix γ(r, r') = Σ_l g_l(r) φ_l*(r') (φ = `orbs`,
g = `weightedOrbs`): its direct potential J, its exchange operator K
compressed with ACE, and the application of J − K_ACE to a vector.

**Direct potential.** `SysInteraction_FillDirectPotential(pot, orbs,
weightedOrbs, time, factor_, bt1_, bt2_)`: J = V[factor ρ],
ρ = Σ_l φ_l* g_l, from one solve of the total source (the potential is
linear in the source; the sources are accumulated in the order of the
orbitals; `factor_` defaults to 1, an exact scaling, as is 2).

**Adaptively compressed exchange.**
`SysInteraction_FillExchangeAce(zeta, positiveDefiniteQ, orbs, weightedOrbs,
time, bt1_, bt2_)` compresses the exchange operator of γ,

```
K v = Σ_l V[φ_l* v] g_l,   ξ_i = K φ_i,   M = ⟨φ|ξ⟩ = L Lᴴ,   ζ = ξ L⁻ᴴ,
K_ACE = Σ_k |ζ_k⟩⟨ζ_k|  (= K on span{φ}, Hermitian, positive semidefinite),
```

with n(n+1)/2 potential solves with `conjSymmetricQ` (the swapped pair by
`ConjugateInteractionPotential`), else n². If M is numerically not positive
definite it uses the pseudo-inverse on the eigenvalues above 1e-12 λ_max and
returns `positiveDefiniteQ = .false.` (the caller warns).

**Application.** `SysInteraction_AddMeanField(w, scratch, z, directPotential,
zeta)`: w ← w + J z − Σ_k ζ_k ⟨ζ_k|z⟩ with J z by
`MultiplyWithInteractionPotential` in `scratch`: the n inner products in
parallel over k, then one fused pass, per element
((w + J z) − o_1 ζ_1) − o_2 ζ_2 − …; reentrant, parallel when called outside
a parallel region.

Users:
- the SCF (`groundSolver.scf.stdImpl`): J of the closed-shell density
  (g = φ, `factor_ = 2`) in every iteration; with `aceQ` K_ACE of the
  occupied orbitals of one spin (g = φ), and every Fock action applies J and
  K_ACE in one `AddMeanField` call;
- the frozen mean-field reference of the Cayley exact part of Method
  (`method.mb.orbBased.exactFieldFreeMeanFieldQ`): J of the spatial 1-RDM
  ρS (g = φ ρS), K_ACE per body type (g = φ R with the 1-RDM block R), built
  once and applied per vector with the ζ columns of the body type.

### Public Data

| Variable | Type | Meaning |
|----------|------|---------|
| `SysInteraction_timeIndependentQ` | `logical` | True if Ŵ is time-independent (enables caching) |
| `SysInteraction_bodyTypeIndependentQ` | `logical` | True if Ŵ treats all species identically |

---

## Grid Back-Ends

### 1. Lattice (Discrete)

**Use case:** Hubbard models, tight-binding systems

**Physics:** On-site density-density interaction
```
Ŵ = U Σᵢ n↑ᵢ n↓ᵢ
```

**Implementations:**
- `OnSite/StdImpl`: Direct V(i) = U·ρ(i)

**JSON path:** `sysInteraction.lattice.onSite`

### 2. Linear (1D Real-Space)

**Use case:** 1D model systems, trapped gases

**Physics:** Convolution with interaction kernel
```
V(x) = ∫ w(|x-x'|) ρ(x') dx'
```

**Implementations:**
- `SoftYukawa/StdImpl`: O(N²) direct convolution
- `SoftYukawa/Fftw3`: O(N log N) FFT convolution

**Kernel form:**
```
w(r) = strength × exp(-dampening×r) / (√(r² + softening1²) + softening2)
```

**JSON path:** `sysInteraction.linear.softYukawa`

### 3. Ylm (Spherical Harmonics)

**Use case:** Atoms, molecules with spherical symmetry

**Physics:** Coulomb 1/r expanded in multipoles
```
1/|r₁-r₂| = Σₗₘ (4π/(2l+1)) × (r<ˡ/r>ˡ⁺¹) × Yₗₘ*(Ω₁)Yₗₘ(Ω₂)
```

**Implementations:**
| Variant | Complexity | Notes |
|---------|------------|-------|
| `StdImpl` | O(N²) | Reference, Green's function |
| `TwoScan` | O(N) | Forward/backward prefix sums |
| `FullEq` | O(N³) | FEDVR + LU per call |
| `FullEqEcs` | O(N³) | FEDVR-ECS contour + complex LU per call (requires `grid.ylm.fedvrEcs`) |
| `BlockEq` | O(nE² + nE×nLoc²) | FEDVR + element inverses and Schur LU precomputed per l |

**JSON path:** `sysInteraction.ylm.coulomb`

**Layout module data** (set in `SysInteraction_Ylm_Fabricate`):
| Name | Meaning |
|------|---------|
| `SysInteraction_Ylm_lmax`, `SysInteraction_Ylm_mmax` | Caps of the source/potential expansion |
| `SysInteraction_Ylm_nChannels` | Number of (l,m) channels, `Grid_Ylm_NChannels(lmax, mmax)` |
| `SysInteraction_Ylm_nPoints` | Length of a source/potential array, `nChannels * Grid_Ylm_nRadial` |

**Potential of one m:** `SysInteraction_Ylm_FillInteractionPotentialM(potM,
srcM, m, time)` maps the source channels of one m (channel matrix
`srcM(:, iL)`, l = |m| + iL − 1, radial metric weights included) of a
body-type-independent interaction to the potential channels by one radial
Poisson solve per l
(`SysInteraction_Ylm_FillInteractionPotentialRadial`); channels with a
negligible source (all |ρ_lm| < 1e-14) are zero. It is parallel over l unless
called inside a parallel region. The definite-m node path of Method calls it
per pair (on the channel matrix of the pair's M).
`SysInteraction_FillInteractionPotential` runs one parallel loop over all
(l, m) channels. Both go through the same per-channel routine
(`FillInteractionPotentialLm`: Poisson solve or zero), so they share the solve
and the threshold.

### 4. Fem (3D Cartesian Spectral Elements)

**Use case:** molecules on `grid.fem` (Coulomb e-e interaction of the
orbital pairs, Hartree and exchange potentials).

Sources and potentials are nodal arrays of `Grid_nPoints` entries. The family
binds the parts shared by nodal kernels:

| Pointer | Fem implementation |
|---------|--------------------|
| `FillInteractionSrc(src, a, b)` | src = conjg(a) · b · w — exactly the Galerkin load of the pair density ρ = conjg(a) b |
| `MultiplyWithInteractionPotential` | point-wise V · orb (reentrant: OpenMP-parallel at top level, e.g. the mean-field reference of the Cayley step and ACE, serial inside Method's OpenMP loops) |
| `ConjugateInteractionPotential` | conjg(V) |

**Coulomb** (`sysInteraction.fem.coulomb`) binds
`FillInteractionPotential`: V = strength × `Grid_Fem_SolvePoisson(src)`, the
Galerkin solution of −lap V = 4πρ (K_II V = 4π load − K_IB V_B) with
free-space box values V_B from the multipole expansion of the source up to
`grid.fem.multipoleLmax` (see `src/Grid/A_Grid.md`), solved exactly by fast
diagonalization. A real source costs one FDM solve, a complex one two (real
and imaginary parts separately). The kernel and the weights are real, so
`conjSymmetricQ = .true.` (the swapped pair is the conjugate, bitwise), and
`timeIndependentQ = bodyTypeIndependentQ = .true.`; `time`, `bt1_` and `bt2_`
are ignored (Method passes `bt1_ = bt2_ = −1`). The solver is linear in the
source, so the SCF Hartree term may use the total density in one solve.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `sysInteraction.fem.coulomb.strength` | real | 1.0 | Prefactor of 1/\|r1 − r2\| |

**Cost** (2 threads, N = 493039): 45 ms per real source, 91 ms per complex
source. Over orbital pairs Method needs n(n+1)/2 solves per evaluation
(conjSymmetricQ).

**Test:** `test/componentTests/SysInteraction/Fem/Coulomb/T_SysInteraction_Fem_Coulomb`
(box ±10, nLocals 6, one center at the origin, hMin 0.2, growth 1.6,
N = 493039): Gaussian density e^(−r²) against π^1.5 erf(r)/r, max error
3.4e-6 (r < 4), 4.8e-8 relative at the origin node (limit 2π); Hartree
integral 24.7394343 (exact √2 π^2.5 = 24.7394295, 2e-7 relative);
J(1s, 1s) = 0.625000041 (exact 5/8); swapped pair conjg(V_ab) − V_ba = 0;
strength 0.5 (refabricated): the potential is exactly half of the strength-1
one and J(1s, 1s) = 0.312500021. The density tails of this box lie far below
the multipole warning threshold of `Grid_Fem_SolvePoisson` (1e-4; 1s
density: 3.8e-7), so no warning is printed.

---

## Workflow

### Initialization Sequence

```fortran
call SysInteraction_Fabricate    ! Parse JSON, bind pointers
call SysInteraction_Setup        ! Precompute kernels/plans
```

### Typical Usage (in Method/Hamiltonian)

```fortran
! Build source from orbital product
call SysInteraction_FillInteractionSrc(src, conjg(orb_i), orb_j)

! Solve for interaction potential
call SysInteraction_FillInteractionPotential(pot, src, time)

! Apply to orbital
call SysInteraction_MultiplyWithInteractionPotential(dOrb, pot, orb_k)
```

---

## JSON Configuration Examples

### Hubbard Model (Lattice)
```json
{
  "sysInteraction": {
    "lattice": {
      "onSite": {
        "strength": 4.0,
        "stdImpl": {}
      }
    }
  }
}
```

### 1D Soft-Coulomb (Linear)
```json
{
  "sysInteraction": {
    "linear": {
      "softYukawa": {
        "strength": 1.0,
        "softening1": 0.5,
        "fftw": {}
      }
    }
  }
}
```

### Helium Atom (Ylm)
```json
{
  "sysInteraction": {
    "ylm": {
      "lmax": 6,
      "coulomb": {
        "strength": 1.0,
        "twoScan": {}
      }
    }
  }
}
```

### Molecule (Fem, requires `grid.fem`)
```json
{
  "sysInteraction": {
    "fem": {
      "coulomb": {
        "strength": 1.0
      }
    }
  }
}
```

---

## Adding a New Interaction

1. **Create directory structure:**
   ```
   src/SysInteraction/<Grid>/<NewModel>/
   ├── M_SysInteraction_<Grid>_<NewModel>.f90
   └── S_SysInteraction_<Grid>_<NewModel>.f90
   ```

2. **In the interface module (`M_*.f90`):**
   - Declare any model-specific parameters
   - Declare the `_Fabricate` subroutine interface

3. **In the implementation submodule (`S_*.f90`):**
   - Implement `_Fabricate` to read JSON and bind pointers
   - Implement the actual physics routines

4. **Wire into parent:**
   - Add `use M_SysInteraction_<Grid>_<NewModel>` in parent fabricate
   - Add JSON existence check and dispatch call

5. **Add tests:**
   - Create `test/T_SysInteraction_<NewModel>.f90`
   - Create corresponding JSON fixture

---

## Performance Considerations

| Grid | Recommendation |
|------|----------------|
| Lattice | Always fast (local interaction) |
| Linear | Use `fftw` for N > 500 |
| Ylm | Use `twoScan` for large grids; `blockEq` if many repeated solves |
| Fem | One FDM Poisson solve per pair (2 for complex pairs); parallelism comes from inside the solve (DGEMMs, OpenMP over z) |

### Caching Opportunities

When `timeIndependentQ = .true.` and `bodyTypeIndependentQ = .true.`:
- The interaction potential can be computed once and reused
- Higher-level code (Method) should exploit this

---

## Dependencies

| Module | Usage |
|--------|-------|
| `M_Grid*` | Grid points, weights, nPoints |
| `M_Utils_Json` | Configuration parsing |
| `M_Utils_Fftw3Lib` | FFT for convolution |
| `M_Utils_LapackLib` | LU factorization for FEDVR solvers |
| `M_Utils_Fedvr` | Finite element DVR infrastructure |
| `M_Utils_ConvolutionIntegral` | Direct convolution utility |
| `M_Utils_ConvolutionFftw` | FFT convolution utility |
| `M_Grid_Fem` | Weights and the Galerkin Poisson solve (Fem) |

---

## Common Pitfalls

1. **Forgetting `SysInteraction_Setup`:** Some implementations (FFT, BlockEq)
   require setup to precompute plans/factorizations.

2. **Wrong `lmax` for Ylm:** The potential expansion needs `lmax_pot ≥ 2×lmax_orb`
   to capture the full density product.

3. **Ylm `mmax`:** `sysInteraction.ylm.mmax` caps |m| of sources/potentials
   (default min(2·`grid.ylm.mmax`, lmax_pot), which is exact). Smaller values
   truncate the m-transfer of the interaction.

4. **Missing quadrature weights:** Linear/Ylm source terms include weights;
   Lattice does not (discrete sum).

5. **Body type arguments:** Pass `bt1_`, `bt2_` if simulating multiple species
   with different interactions.

---

## File Summary

| File | Lines | Purpose |
|------|-------|---------|
| `M_SysInteraction.f90` | ~250 | Public interface, procedure pointers, mean-field helper interfaces |
| `S_SysInteraction.f90` | ~60 | Top-level fabrication dispatch |
| `S_SysInteraction_MeanField.f90` | ~230 | `SysInteraction_FillDirectPotential`, `SysInteraction_FillExchangeAce` (ACE), `SysInteraction_AddMeanField` (on the contract) |
| `*/M_*.f90` | ~30-50 | Grid/model interface modules |
| `*/S_*.f90` | ~50-450 | Implementation submodules |

Total: ~25 files, ~2000 lines of physics code.
