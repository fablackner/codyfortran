# Method Module — Architecture & Onboarding Guide

## Overview

The **Method** module defines the unified dynamics interface for quantum many-body
simulations in CodyFortran. It provides the time-derivative contract
`Method_TimeDerivative(dState, state, t)` that all propagators consume, enabling
the framework to evolve arbitrary quantum states through time.

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                            Propagator / Integrator                          │
│                      calls Method_TimeDerivative(dState, state, t)          │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                              M_Method (interface)                            │
│         Method_state(:)  |  Method_Setup  |  Method_TimeDerivative          │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
              ┌───────────────────────┴───────────────────────┐
              ▼                                               ▼
       ┌─────────────┐                               ┌───────────────┐
       │  Sb (single │                               │ Mb (many-body)│
       │   particle) │                               │               │
       └─────────────┘                               └───────────────┘
                                                              │
                          ┌───────────────────────────────────┼───────────────┐
                          ▼                                   ▼               ▼
                   ┌─────────────┐                    ┌─────────────┐  ┌─────────────┐
                   │  OrbBased   │                    │  GridBased  │  │  GemBased   │
                   │ (orbital    │                    │ (direct     │  │ (Gaussian   │
                   │  expansion) │                    │  product)   │  │  expansion) │
                   └─────────────┘                    └─────────────┘  └─────────────┘
                          │
       ┌──────────┬───────┴───────┬────────────┬─────────────────┐
       ▼          ▼               ▼            ▼                 ▼
    ┌──────┐  ┌──────┐       ┌────────┐   ┌────────┐    ┌─────────────────┐
    │ TDCI │  │ TDHX │       │ MCTDHX │   │TD-2RDM │    │TD-2RDM-Static   │
    │      │  │      │       │        │   │        │    │     Orbs        │
    └──────┘  └──────┘       └────────┘   └────────┘    └─────────────────┘
```

---

## Key Concepts

### 1. Packed State Vector

All methods operate on a **single contiguous 1D complex array** `Method_state(:)`.
The internal layout is method-specific:

| Method              | State Layout                                      |
|---------------------|---------------------------------------------------|
| Sb                  | Single wavefunction ψ(r) on the grid              |
| TDCI                | CI coefficients c_I only (static orbitals)        |
| TDHX                | Orbitals φ_i(r) only (single determinant)         |
| MCTDHX              | Orbitals φ_i(r) ⊕ CI coefficients c_I             |
| TD-2RDM             | Orbitals φ_i(r) ⊕ 2-RDM D²_{ijkl}                 |
| TD-2RDM-StaticOrbs  | 2-RDM D²_{ijkl} only (fixed orbitals)             |
| GridBased/Full      | Full many-body wavefunction Ψ(r₁,r₂,...) on grid  |

### 2. Procedure Pointer Binding (Fabrication Pattern)

At runtime, the JSON configuration determines which implementation is bound:

```fortran
! In Method_Fabricate:
if (Json_GetExistence("method.sb")) then
  call Method_Sb_Fabricate     ! binds Sb procedures
else if (Json_GetExistence("method.mb")) then
  call Method_Mb_Fabricate     ! descends into many-body branch
end if
```

Each `*_Fabricate` routine:
1. Reads method-specific JSON parameters
2. Binds the module's procedure pointers (`Method_Setup`, `Method_TimeDerivative`, etc.)
3. Calls the next level's fabricate if further branching is needed

### 3. Time Derivative Contract

The central interface every method must implement:

```fortran
subroutine Method_TimeDerivative(dState, state, time)
  complex(R64), intent(out) :: dState(:)   ! output: -i·H·|Ψ⟩
  complex(R64), intent(in)  :: state(:)    ! input:  |Ψ⟩
  real(R64), intent(in)     :: time        ! current time
end subroutine
```

**Semantics:**
- Computes dΨ/dt = −i·Ĥ(t)·Ψ (ħ=1 units)
- Overwrites `dState`; no aliasing with `state`
- May use cached operators or module-level data

---

## Module Hierarchy

### Top Level: `M_Method` / `S_Method`

| Export                    | Purpose                                          |
|---------------------------|--------------------------------------------------|
| `Method_state(:)`         | The packed quantum state (allocated by Setup)    |
| `Method_Setup`            | Procedure pointer → initialize state & caches    |
| `Method_TimeDerivative`   | Procedure pointer → compute −i·H·|Ψ⟩             |
| `Method_GetEnergy`        | Procedure pointer → compute ⟨Ψ|H|Ψ⟩              |
| `Method_Fabricate`        | Entry point; branches to Sb or Mb                |

### Single-Body Branch: `M_Method_Sb`

Solves the single-particle TDSE under an external potential V(r,t):
```
i ∂ψ/∂t = (T̂ + V̂) ψ
```
No mean-field or exchange terms. Useful for testing grids and basic physics.

### Many-Body Branch: `M_Method_Mb`

Shared metadata for all many-body methods:

| Data                           | Purpose                                    |
|--------------------------------|--------------------------------------------|
| `Method_Mb_nBodyTypes`         | Number of particle species/spin types      |
| `Method_Mb_nBodies(:)`         | Particle count per body type               |
| `Method_Mb_bodyStatistics(:)`  | 'f' (fermion) or 'b' (boson) per type      |
| `Method_Mb_nBodiesStart/End(:)`| Index ranges for packed body-type storage  |

Branches to: **OrbBased**, **GridBased**, or **GemBased**.

---

## OrbBased Methods (Most Feature-Rich)

### Shared Infrastructure: `M_Method_Mb_OrbBased`

Provides orbital indexing, RDM storage, and operator-application interfaces:

| Procedure Pointer                              | Signature / Purpose                        |
|------------------------------------------------|--------------------------------------------|
| `Method_Mb_OrbBased_FillH1`                    | Build 1-body Hamiltonian h¹_{ij}           |
| `Method_Mb_OrbBased_FillH2`                    | Build 2-body Hamiltonian h²_{ijkl}         |
| `Method_Mb_OrbBased_FillRdm1`                  | Extract 1-RDM from state                   |
| `Method_Mb_OrbBased_FillRdm2`                  | Extract 2-RDM from state                   |
| `Method_Mb_OrbBased_ApplyKineticOp`            | T̂ · orbitals (accumulate)                  |
| `Method_Mb_OrbBased_ApplyPotentialOp`          | V̂_ext · orbitals (accumulate)              |
| `Method_Mb_OrbBased_ApplyInteractionOp`        | Ŵ · orbitals for index pair (i₂,j₂)        |
| `Method_Mb_OrbBased_ApplyHartreeFockOp`        | Mean-field + exchange (Hartree–Fock)       |
| `Method_Mb_OrbBased_ApplyCorrelationOp`        | Correlation term `(dOrbs, orbs, invRdm1, rdm2, time, h2_)` |
| `Method_Mb_OrbBased_FillSpatialRdm1` (module procedure) | Spatial 1-RDM (spin-summed for restricted) |
| `Method_Mb_OrbBased_FillInverseRdm1` (module procedure) | Regularized spatial inverse U n/(n²+ε) U† of the 1-RDM |
| `Method_Mb_OrbBased_SaveRdm2Blocks` (module procedure) | Write body-type blocks `rdm2BtXX_YY.in` of a 2-RDM (read by `twoRdmInit.load`) |
| `Method_Mb_OrbBased_ApplySingleBodyOp`         | Combined T̂ + V̂_ext                         |
| `Method_Mb_OrbBased_TimeDerivativeOrbsLin`     | Linear (static-orbital) time derivative    |
| `...TimeDerivativeCoeffsPlusOrbsNonLin`        | Non-linear part for MCTDHX-like schemes    |

### Concrete Methods

| Method           | Orbitals   | Coefficients | Description                                |
|------------------|------------|--------------|-------------------------------------------|
| **TDCI**         | Static     | Dynamic      | Full CI; exact within orbital basis        |
| **TDHX**         | Dynamic    | Implicit     | Single determinant; mean-field + exchange  |
| **MCTDHX**       | Dynamic    | Dynamic      | Multiconfiguration; Dirac–Frenkel VP       |
| **TD-2RDM**      | Dynamic    | 2-RDM        | Propagates 2-RDM with adaptive orbitals    |
| **TD-2RDM-Static**| Static    | 2-RDM        | Propagates 2-RDM only; fixed basis         |

### Split-Step Propagation (MCTDHX, TD-2RDM)

`Method_Mb_OrbBased_TimeDerivativeOrbsLin` and
`Method_Mb_OrbBased_TimeDerivativeCoeffsPlusOrbsNonLin` split the time
derivative for `propagator.splitStep` (stiff linear part by SIL, the rest by
an explicit integrator):

- **linear**: dφ = −i h φ for every orbital, *unprojected*; the many-body part
  (CI coefficients or 2-RDM) is constant. Propagating all orbitals with
  exp(−i h t) at fixed many-body part is the exact one-body evolution.
- **non-linear**: only the interaction — CI coefficients with h², or the
  2-RDM with [h², D²] plus the collision operator; orbitals with the projected
  correlation term.

The sum is a gauge transform of the full derivative (in-space orbital rotation
by h¹, compensated in the many-body part). Adding h¹ to the many-body part of
the non-linear step would count the one-body Hamiltonian twice.

### Exact Field-Free One-Body Propagation

With `"exactFieldFreeQ": true` in `method.mb.orbBased` (MCTDHX, TD-2RDM), the
linear step is split once more: the field-free h0 = T + V_ext is propagated
exactly and the linear integrator only sees the gauge coupling
(`Method_Mb_OrbBased_ApplyLinOp` binds the gauge operator instead of the full
one-body operator). The methods bind `Method_PropagateExactPart`, and the
split-step propagators compose every linear step from exact and integrator-1
substeps (`propagator.splitStep.exactPartScheme`; the default `strang` is
exp(−i h0 Δt/2) · integrator 1(Δt) · exp(−i h0 Δt/2), local error O(Δt³), as
the outer Strang splitting). Only MCTDHX and TD-2RDM (without weak orbitals)
bind `Method_PropagateExactPart`; for TDCI, TDHX and TD-2RDM-StaticOrbs the
flag has no effect.

`S_Method_Mb_OrbBased_FieldFree.f90` builds h0 on the first call: it applies
the kinetic and potential operators to radial unit vectors placed in all
(l,m) channels at once, symmetrizes each block with the radial metric
weights (w r² on FEDVR) and diagonalizes it per l (real symmetric, m
independent). Propagation is then two real GEMMs per l on the (orbital, m)
columns of that block, with phases exp(−i E Δt) (or exp(−E Δτ) in imaginary
time). With declared m (`method.mb.orbBased.mOfOrbs`, see "Declared m of the
orbitals" below) only the column of the declared m of every orbital is
transformed (one of the 2 min(l, mmax) + 1 columns; built once with the
blocks), also by the phi functions and the spectral cut; without it all
columns are.

Why: on FEDVR grids h0 is very stiff (centrifugal term l(l+1)/2r² at the
first radial point, ~10⁵–10⁶ a.u. at lmax = 32), so SIL needs hundreds of
Krylov vectors per step; the gauge coupling alone needs about ten. For the Be
HHG runs this is 5–8× faster per time step with the default Strang
composition, but its splitting error ([h0, A·p] = iA ∂zV, largest near the
nucleus) is 7–35× larger than the SIL time-step error at the peak of a
4×10¹⁴ W/cm² pulse. Use the fourth-order composition
`propagator.splitStep.exactPartScheme = "blanesMoan6"` (A_Propagator.md): its
errors coincide with those of SIL at the same Δt (the outer Strang error
dominates) at 2.5–5× lower cost; see `app/be/README.md`.

`"exactFieldFreeGpuQ": true` runs the two GEMMs per l on the GPU with cuBLAS
(`M_Utils_CublasLib`, build with `cmake -DCODY_USE_CUDA=ON`): the eigenvectors
stay on the device, each call gathers all blocks on the host, sends them once,
applies the phases on the host between the two transformations and returns
them. Results agree with the CPU path to round-off.

**Spectral cut.** `"exactFieldFreeCut": Ec` (MCTDHX only, default 0 = off)
restricts the orbitals to the h0 eigenvectors with E < Ec, P = Σ_{E<Ec} |u⟩⟨u|.
The highest FEDVR eigenvalues (~10⁵ a.u., from the dense Lobatto points at
the element boundaries) are unphysical, but they make the Strang composition
of h0 and the gauge coupling inaccurate. With the cut the exact part
transforms with the kept eigenvectors only (it projects onto P), and
`Method_Mb_OrbBased_ApplySpectralCut` applies P to the linear and nonlinear
orbital derivatives of MCTDHX, so P is invariant under the whole step. The pointer is
bound in `Method_Mb_OrbBased_Fabricate` to
`Method_Mb_OrbBased_ProjectBelowFieldFreeCut` (two GEMMs with the few dropped
eigenvectors, P = 1 − U_H U_Hᵀ) or to a no-op. The cut changes the model:
relax the ground state with the same cut (the initial projection of an uncut
state is a jump). For Be HHG, Ec = 1000 with `strang` is ~2× faster per step
than `blanesMoan6Swap` and shifts harmonics 21–60 by 0.02–0.04 decades
(`app/be/README.md`). T_Ne3d_11 checks that no weight appears above the cut.

**Phi functions.** MCTDHX with `exactFieldFreeQ` also binds
`Method_ApplyExactPartPhi(res, vecs, dt)`: res = Σ_m φ_m(−i H0 dt) vecs(:, m+1)
with φ_0(z) = e^z and φ_m(z) = (φ_{m−1}(z) − 1/(m−1)!)/z, H0 = h0 on the
orbitals and 0 on the CI coefficients (φ_m(0) = 1/m!). The orbital part,
`Method_Mb_OrbBased_ApplyFieldFreePhi`, transforms every vector into the h0
eigenbasis (kept eigenvectors only), weights it with φ_m(−i E dt) (Taylor
series for |z| < 1, recurrence otherwise) and transforms the sum back once
(CPU only). It is consumed by the exponential integrators
(`propagator.exponential`, src/Propagator/A_Propagator.md); real time only.

Requirements (checked on the first exact step): Hermitian (non-ECS) Ylm grid, time- and
body-type-independent kinetic and potential, central potential
(`SysPotential_Ylm_lmax = 0`), no TD-2RDM weak orbitals. Memory: one real
nRadial² matrix per l. T_Ne3d_09 (Strang) and T_Ne3d_10 (Blanes–Moan) check
the agreement with the SIL scheme of T_Ne3d_08.

### TD-2RDM Active/Weak Natural-Orbital Split

The TD-2RDM orbital equations contain the inverse 1-RDM, so weakly occupied
natural orbitals amplify every error of the reconstructed 2-RDM by 1/n_min:
with more orbitals the equations get stiff and noisy. `method.mb.orbBased.td2rdm.nWeakOrbs` (int, default 0 =
off, standard TD-2RDM) makes the `nWeakOrbs` lowest natural orbitals of each
body type (spin-summed 1-RDM if restricted) *weak*: they get no Q-space drive.
Typical: `nWeakOrbs = 2` for 6–7 orbitals.

With G_u = Σ_j hφ_j D1(j,u) + Γ_u (Γ_u the two-body part of the correlation
operator) the orbitals evolve with

```
dφ_i/dt = −i Q Σ_u G_u M(u,i),     M = Σ_{u active} |u⟩ n_u/(n_u²+ε) ⟨u|
```

For no weak orbitals this is the usual −iQ(hφ_i + Σ_u Γ_u M(u,i)); with weak
orbitals D1 M projects on the active natural orbitals, so weak natural
orbitals only rotate in span. **Energy conservation** holds for every
Hermitian M (field free): dE/dt = 2 Re Σ_i ⟨φ̇_i|G_i⟩ = 2 Re(i Tr(M A)) = 0
with A(u,i) = ⟨G_u|Q|G_i⟩ Hermitian (the trace of a product of two Hermitian
matrices is real); the 2-RDM part conserves E by itself.

**Split step.** The linear part rotates the orbitals with
H_L = (1−P_w) h (1−P_w) + P_w h P_w (P_w the grid projector on the weak
natural orbitals, frozen during the step); the non-linear part adds the
missing active-weak blocks of h1 to the 2-RDM, −i[Π h1 (1−Π) + (1−Π) h1 Π, D2]
(Π the weak projector in coefficient space). The sum is a gauge transform of
the equation above. Γ M is computed by `Method_Mb_OrbBased_ApplyCorrelationOp`
with the truncated inverse passed as `invRdm1` (instead of the full
regularized inverse from `Method_Mb_OrbBased_FillInverseRdm1`).

**`Method_Mb_OrbBased_Td2rdm_RefreshWeakOrbs`** (`M_Method_Mb_OrbBased_Td2rdm`,
no-op for `nWeakOrbs = 0`) recomputes P_w from `Method_state`. `Method_Setup`
and every non-linear evaluation refresh it (the last RK stage approximates the
end of the step; validated for `splitStep.order2` with an RK non-linear
integrator only). A program that changes the state outside the propagator
(absorber, purification, kick) must call it before the next
`Propagator_Propagate`, otherwise the split step drops to first order; for
other propagators call it after every `Propagator_Propagate`.
`nWeakOrbs` must be smaller than every `nOrbs`.
Tests: `T_He1d_06_TimePropagationTd2rdmWeakOrbs` (restricted orbitals, RK4:
exact energy conservation and natural occupations against the unrestricted
reference) and `...WeakOrbsSplit` (unrestricted, split step against that
reference).

### Declared m of the Orbitals and the Definite-m Ylm Fast Path

`"mOfOrbs": [...]` in `method.mb.orbBased` (optional, Ylm grids only) declares
the m of every spin orbital, in the flat layout of `nOrbs` (all orbitals of
body type 1, then body type 2; restricted: equal halves), e.g. Ne 1s 2s 2p 3s
3p with `nOrbs = [9, 9]`: `[0,0,-1,0,1,0,-1,0,1, 0,0,-1,0,1,0,-1,0,1]`. It
states a symmetry of the problem (orbitals of an atom in a linearly polarized
field, or m = 0 orbitals on a `grid.ylm.mmax = 0` grid): every operator
conserves m, so each orbital stays in the Ylm channels of its m. It is
validated in `Method_Mb_OrbBased_Fabricate` (|m| ≤ `grid.ylm.mmax`, one entry
per spin orbital; `Orbs_Fabricate`: equal spin halves if restricted) and
stored in `Method_Mb_OrbBased_mOfOrb` (`Method_Mb_OrbBased_definiteMQ`). It is
independent of `orbsInit.ylm.hydrogenLike.m`, which only sets initial
orbitals. Without it everything runs on the generic path (orbitals of
mixed m, e.g. a field that is not z-polarized); on a Ylm grid the fabricate
step then prints "mOfOrbs not declared: generic correlation operator and h2".

With the declaration:
- **Initialization**: the orbitals are projected on their m once, at the end
  of `OrbsInit_Initialize` (`Orbs_ProjectOnSymmetry`, `Grid_Ylm_ProjectOnM`);
  this removes round-off admixtures of loaded orbitals and stops with an
  error if more than 1e-8 of an orbital's norm would be removed (orbitals
  that do not have the declared m, e.g. mOfOrbs in another order than the
  files). The only other call is the end of `EmbedOrbitals` of the Be/Ne
  programs, which replaces this initialization when orbitals are embedded
  from a grid with smaller lmax (`program.lmaxIn`).
- **Exact zeros**: inner products of orbitals of different m are exactly 0
  and every operator is Δm = 0, so the other channels stay exactly zero. The
  only operations that could seed admixtures are eigen-decompositions in
  orbital space with degenerate blocks of different m; they act per block of
  `Orbs_blockOfOrb` (body type and m; `LapackLib_DiagonalizeHermitianBlocks`):
  the natural orbitals of the inverse 1-RDM
  (`Method_Mb_OrbBased_FillInverseRdm1`), the Loewdin steps of the Be/Ne
  programs and the alignment of the Ylm-optimized ground solvers
  (`Orbs_AlignOnReference(..., blockOf_)`; the standard ground solvers, which
  diagonalize on the full grid, are rejected with mOfOrbs), and the weak
  natural orbitals of TD-2RDM. T_Ne3d_11 (unrestricted) and T_Ne3d_13
  (restricted) check the exact zeros after propagation. Without the
  declaration the blocks are the body types: restricted orbitals form one
  block and the decompositions are the plain ones bitwise; for unrestricted
  orbitals the inverse 1-RDM is decomposed per spin instead of as a whole,
  which changes results at round-off only.
- **Binding at fabricate**: for Hermitian Ylm grids with a
  body-type-independent Ylm interaction (`YlmNodesFastPathQ`),
  `Method_Mb_OrbBased_ApplyCorrelationOp` is bound to
  `ApplyCorrelationOpYlmNodes` and `Method_Mb_OrbBased_FillH2` to
  `FillH2YlmNodes` (also used by `Method_GetEnergy`); the exact field-free
  part transforms only the declared columns (above).

The node path evaluates orbitals (in the channels of their declared m), pair
densities and pair potentials on the Gauss–Legendre nodes in cos θ
(`Grid_Ylm_ToNodes`, `Grid_Ylm_FromNodes` and their channel-matrix forms).
The pair (l2, k2) carries M = m_k2 − m_l2; only the terms m_k1 + M = m_j1 in
the correlation sum and m_j1 + M = m_i1 in h2 are kept. The dropped h2 terms
vanish by the φ integral, the dropped correlation terms for an M-conserving
many-body state (their 2-RDM weights vanish). Pairs with |M| above the |m|
cap of the potentials are dropped, as in the generic operator. Every product
is pointwise, and the node count makes all angular integrals exact, so
results equal the generic operator to round-off (T_Ne3d_14 runs the
restricted T_Ne3d_13 without mOfOrbs: dipoles agree to round-off) while the cost per pair drops from n
spatial products to O(1) transforms.

The pair list (ordered pairs of equal body type within the |m| cap) and the
target groups (`SetupYlmNodes`) are built once from the declaration in
`Method_Mb_OrbBased_Setup`, which every orbital-based `Method_Setup` calls
first; every call only refills the 2-RDM weights of the groups
(`FillTargetGroupWeights`). The pair potentials are computed in parallel over
the unordered pairs on the channel matrix of their single M
(`FillPairPotentialsYlmNodes`: `Grid_Ylm_NodesToChannels`, radial metric
weights, `SysInteraction_Ylm_FillInteractionPotentialM` (one radial Poisson
solve per l, the same per-channel solve and source threshold as the generic
`SysInteraction_FillInteractionPotential`), `Grid_Ylm_ChannelsToNodes`) and
stored on the nodes; the swapped pair uses the conjugate. The pointwise sums
(`AccumulateTargetSumsYlmNodes`) are grouped by the target m of the selection
rule: for every group the products X = V_iP · φ_k of its columns (iP, k) are
formed blockwise over cache-sized blocks of nodes and enter two ZGEMMs,
dφ_j += X · weights and h2 = (φ_i w)^H · X, which touch only the 2-RDM
weights allowed by the rule (Ne with 1s 2s 2p 3s 3p: 1921 of 6561);
`FillH2YlmNodes` forms only the h2 columns.

One MCTDHX non-linear derivative of Ne with 9 orbitals (1s 2s 2p 3s 3p,
lmax 40, mmax 1, 761 radial points, four threads on 2 cores): 0.38 s (0.28 s
correlation operator), against 3.7 s with the generic operator. One energy
evaluation (`Method_GetEnergy`, h2 by `FillH2YlmNodes`): 0.25 s instead of
2.2 s.

---

## GridBased Methods

Represents the full many-body wavefunction Ψ(r₁,r₂,...,rₙ) directly on a tensor
product grid. Exponentially expensive but numerically exact for small N.

| Procedure Pointer                              | Purpose                            |
|------------------------------------------------|------------------------------------|
| `Method_Mb_GridBased_ApplyKineticOp`           | Laplacian on 1D grid slice         |
| `Method_Mb_GridBased_ApplyPotentialOp`         | V(r) on 1D grid slice              |
| `Method_Mb_GridBased_ApplyInteractionOp`       | W(r₁,r₂) on 2D grid slice          |

Currently supports: **Full** (complete tensor-product representation).

---

## GemBased Methods (Placeholder)

Intended for Gaussian-Expansion Methods (GEM) where the wavefunction is expanded
in explicitly correlated Gaussian basis functions. Infrastructure is in place
but concrete implementations are not yet finalized.

---

## JSON Configuration Examples

### Single-Body (Sb)
```json
{
  "grid": { "linear": { "const": { "nPoints": 128, "xMin": -10, "xMax": 10 } } },
  "sysPotential": { "linear": { "harmonic": { "omega": 1.0 } } },
  "method": { "sb": {} }
}
```

### Many-Body MCTDHX (2 fermions, 4 orbitals)
```json
{
  "method": {
    "mb": {
      "nBodyTypes": 1,
      "nBodies": [2],
      "bodyStatistics": ["f"],
      "orbBased": {
        "nOrbs": [4],
        "mctdhx": {}
      }
    }
  }
}
```

### Many-Body TD-2RDM with Static Orbitals
```json
{
  "method": {
    "mb": {
      "nBodyTypes": 2,
      "nBodies": [2, 2],
      "bodyStatistics": ["f", "f"],
      "orbBased": {
        "nOrbs": [4, 4],
        "td2rdmStaticOrbs": {}
      }
    }
  }
}
```

---

## Adding a New Method

1. **Create directory structure:**
   ```
   src/Method/Mb/OrbBased/MyMethod/
       M_Method_Mb_OrbBased_MyMethod.f90   (interface module)
       S_Method_Mb_OrbBased_MyMethod.f90   (submodule implementation)
   ```

2. **Define the interface module** with a `*_Fabricate` subroutine.

3. **Implement the submodule:**
   - Bind `Method_Setup` → allocate `Method_state`, initialize
   - Bind `Method_TimeDerivative` → compute −i·H·|Ψ⟩
   - Optionally bind `Method_GetEnergy`

4. **Register in parent fabricate:**
   Add a branch in `S_Method_Mb_OrbBased.f90`:
   ```fortran
   else if (Json_GetExistence("method.mb.orbBased.myMethod")) then
     call Method_Mb_OrbBased_MyMethod_Fabricate
   ```

5. **Update CMakeLists.txt** to include the new source files.

---

## Testing

Tests live in `test/` and follow the naming convention `T_Method_*.f90` with
corresponding `T_Method_*.json` configuration files.

Run tests:
```bash
cd build && ctest -R T_Method
```

---

## Common Pitfalls

1. **Forgetting to initialize `dState` to zero** — `TimeDerivative` overwrites,
   but operator-application routines *accumulate*. Ensure proper initialization.

2. **Body-type index mismatch** — When looping over orbitals, always check
   `bodyTypeOfOrb(i) == bodyTypeOfOrb(j)` before computing matrix elements
   between body-type-restricted quantities.

3. **Restricted vs. unrestricted orbitals** — In restricted calculations,
   `Orbs_nOrbsInState < Method_Mb_OrbBased_nOrbsSum`. The code duplicates
   blocks in h1/h2 matrices; be careful with index arithmetic.

4. **Regularization** — `Method_Mb_OrbBased_regularizationParameter` prevents
   division by zero in natural-orbital transformations. Keep it small (1e-10)
   but not zero.

---

## Key Equations

### TDSE (all methods)
```
i ∂|Ψ⟩/∂t = Ĥ(t) |Ψ⟩
```

### Energy Expectation
```
E(t) = ⟨Ψ(t)| Ĥ(t) |Ψ(t)⟩ = Tr[ρ¹ h¹] + ½ Tr[ρ² h²]
```

### MCTDHX Equations of Motion
```
i ∂c_I/∂t = Σ_J ⟨Φ_I| Ĥ |Φ_J⟩ c_J

i ∂|φ_j⟩/∂t = (1 - P̂) Σ_k (ρ¹)⁻¹_{jk} ĥ^{mf}_k |φ_k⟩
```
where P̂ is the projector onto the occupied orbital space.

### TD-2RDM (BBGKY Hierarchy)
```
i ∂D²/∂t = [h¹, D²] + [Ŵ, D³]   (D³ approximated via reconstruction)
```

---

## File Summary

| File                                  | Role                                      |
|---------------------------------------|-------------------------------------------|
| `M_Method.f90`                        | Top-level interface; state & pointers     |
| `S_Method.f90`                        | Fabricate dispatcher (Sb vs Mb)           |
| `Sb/M_Method_Sb.f90`                  | Single-body interface                     |
| `Sb/S_Method_Sb.f90`                  | Single-body implementation                |
| `Mb/M_Method_Mb.f90`                  | Many-body shared data (body types, etc.)  |
| `Mb/S_Method_Mb.f90`                  | Many-body fabricate dispatcher            |
| `Mb/OrbBased/M_Method_Mb_OrbBased.f90`| Orbital-based shared data & interfaces    |
| `Mb/OrbBased/S_Method_Mb_OrbBased.f90`| Orbital-based operator implementations    |
| `Mb/OrbBased/S_Method_Mb_OrbBased_FieldFree.f90` | Exact field-free one-body propagation (Ylm) |
| `Mb/OrbBased/Tdci/*`                  | Full CI with static orbitals              |
| `Mb/OrbBased/Tdhx/*`                  | Time-dependent Hartree–Fock/Exchange      |
| `Mb/OrbBased/Mctdhx/*`                | Multiconfiguration TDHX                   |
| `Mb/OrbBased/Td2rdm/*`                | TD-2RDM with dynamic orbitals             |
| `Mb/OrbBased/Td2rdmStaticOrbs/*`      | TD-2RDM with static orbitals              |
| `Mb/GridBased/M_Method_Mb_GridBased.f90` | Grid-based interfaces                  |
| `Mb/GridBased/S_Method_Mb_GridBased.f90` | Grid-based operator implementations    |
| `Mb/GridBased/Full/*`                 | Full tensor-product grid method           |
| `Mb/GemBased/*`                       | Gaussian expansion (placeholder)          |

---

## Contact & Further Reading

- Main documentation: `README.md` in repository root
- Build instructions: `CMakeLists.txt`, `setBuildVars.sh`
- Related modules: `M_Propagator`, `M_Orbs`, `M_Coeffs`, `M_TwoRdm`
