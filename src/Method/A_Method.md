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

Solves the single-particle TDSE under an external potential V(r,t) and, if a
`sysGauge` block exists, the laser coupling Ĥ_A(t) of SysGauge (velocity gauge
on Ylm, length gauge `sysGauge.fem.length` on grid.fem; bound at fabricate like
`Method_Mb_OrbBased_ApplyGaugeOp`):
```
i ∂ψ/∂t = (T̂ + V̂ + Ĥ_A(t)) ψ
```
No mean-field or exchange terms. Useful for testing grids and basic physics.
T_Method_Sb_LengthGaugeFem checks dψ/dt = −i (T + V + E(t) e_pol·r) ψ at the
envelope peak (e_pol = z there).

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
| `Method_Mb_OrbBased_H1TimeIndependentQ` (module procedure) | Whether h¹ = T̂ + V̂_ext (+ Ĥ_A if `sysGauge` is configured) is time-independent (the `*_timeIndependentQ` flags at call time); TDCI and TD-2RDM-StaticOrbs refill h¹ at every derivative otherwise |
| `Method_Mb_OrbBased_FillOneBodyMatrix` (module procedure) | h¹(i,j) = ⟨φ_i\|hφ_j⟩ of a one-body operator from its action on the stored orbitals (equal body types, nOrbsSum layout, restricted: both spin blocks); the matrix of `FillH1` and of the mean-field reference |
| `Method_Mb_OrbBased_SetupFieldFreeMeanField`, `_BuildFieldFreeMeanField`, `_ApplyFieldFreeMeanField`, `_ApplyFieldFreeMeanFieldOnOrb`, `_FieldFreeMeanFieldBuiltQ` (module procedures, `S_Method_Mb_OrbBased_FieldFreeMeanField`) | Frozen mean-field reference U = J0 − K0 of the Cayley exact part (`exactFieldFreeMeanFieldQ`): setup (releases U; called at the end of the Cayley setup), build from orbitals and 1-RDM, dOrbs += U orbs and ⟨φ\|Uφ⟩ (MCTDHX nonlinear part), dOrb += U orb per vector once built and active, a no-op otherwise (Cayley step), query; the setup and the per-vector form are internal plumbing of the Cayley exact part (called only by `S_Method_Mb_OrbBased_FieldFreeCayley`, not for programs) |
| `Method_Mb_OrbBased_SetupFieldFree`, `_PropagateFieldFree`, `_ApplyFieldFreePhi`, `_ProjectBelowFieldFreeCut` | Exact field-free part (pointers bound at fabricate with `exactFieldFreeQ`: `*Eigen` on Ylm grids, `*Cayley` setup/propagation on grids with `SysKinetic_SolveShifted`, phi and cut null there; all four null without `exactFieldFreeQ`); the setup is called only by `Method_Mb_OrbBased_Setup` when `Method_PropagateExactPart` is bound (error stop if the setup pointer is null then) |
| `Method_Mb_OrbBased_ApplySingleBodyOp`         | Combined T̂ + V̂_ext                         |
| `Method_Mb_OrbBased_TimeDerivativeOrbsLin`     | Linear (static-orbital) time derivative    |
| `...TimeDerivativeCoeffsPlusOrbsNonLin`        | Non-linear part for MCTDHX-like schemes    |
| `Method_Mb_OrbBased_OrthonormalizeKeepingState` | Loewdin orthonormalization of the orbitals with the compensating CI transformation, state unchanged (MCTDHX; error stop otherwise); after absorber masks (`absorber.keepOrthonormalQ`) and in imaginary time |

### Concrete Methods

| Method           | Orbitals   | Coefficients | Description                                |
|------------------|------------|--------------|-------------------------------------------|
| **TDCI**         | Static     | Dynamic      | Full CI; exact within orbital basis (h1 refilled every derivative unless `Method_Mb_OrbBased_H1TimeIndependentQ`: T, V or a configured gauge time-dependent) |
| **TDHX**         | Dynamic    | Implicit     | Single determinant; mean-field + exchange  |
| **MCTDHX**       | Dynamic    | Dynamic      | Multiconfiguration; Dirac–Frenkel VP       |
| **TD-2RDM**      | Dynamic    | 2-RDM        | Propagates 2-RDM with adaptive orbitals    |
| **TD-2RDM-Static**| Static    | 2-RDM        | Propagates 2-RDM only; fixed basis (h1 refilled as for TDCI) |

T_Method_Tdci_LengthGaugeFem and T_Method_Td2rdmStaticOrbs_LengthGaugeFem
check the h1 refill with the length-gauge laser of grid.fem: after a time
derivative at the envelope peak h1(t) − h1(0) = E(t) D (dipole matrix of the
static orbitals; for TD-2RDM-Static in both spin blocks, h2 unchanged).

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

`S_Method_Mb_OrbBased_FieldFreeEigen.f90` builds h0 once in
`Method_Mb_OrbBased_Setup` (`Method_Mb_OrbBased_SetupFieldFree`): it applies
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

Requirements (checked at setup): Hermitian (non-ECS) Ylm grid, time- and
body-type-independent kinetic and potential, central potential
(`SysPotential_Ylm_lmax = 0`), no TD-2RDM weak orbitals. Memory: one real
nRadial² matrix per l. T_Ne3d_09 (Strang) and T_Ne3d_10 (Blanes–Moan) check
the agreement with the SIL scheme of T_Ne3d_08.

**Binding.** The four entry points `Method_Mb_OrbBased_SetupFieldFree`,
`_PropagateFieldFree`, `_ApplyFieldFreePhi` and `_ProjectBelowFieldFreeCut`
are procedure pointers owned by `Method_Mb_OrbBased_Fabricate`: it nullifies
all four, and with `exactFieldFreeQ` its `BindExactFieldFree` binds those of
the variant the grid supports, the others staying null (SysKinetic must be
fabricated before Method):

| Grid | Implementation | Phi functions | Spectral cut | Imaginary time |
|---|---|---|---|---|
| `grid.ylm` | `*Eigen` (`S_Method_Mb_OrbBased_FieldFreeEigen.f90`, dense h0 blocks per l, exp(−i h0 dt)) | yes | yes | yes |
| other grids whose kinetic backend binds `SysKinetic_SolveShifted` (grid.fem) | `*Cayley` (`S_Method_Mb_OrbBased_FieldFreeCayley.f90`) | no (pointer null) | no | no |
| otherwise | error stop at fabricate | | | |

The variant is decided only there; the later checks test the bound pointers:
`exactFieldFreeMeanFieldQ` requires `_PropagateFieldFree` to be the Cayley
step, `exactFieldFreeCut` requires `_ProjectBelowFieldFreeCut`, and MCTDHX
binds `Method_ApplyExactPartPhi` only if `Method_Mb_OrbBased_ApplyFieldFreePhi`
is associated (eigen mode). `Method_PropagateExactPart` and
`Method_ApplyExactPartPhi` are owned by `Method_Fabricate`, which nullifies
them before the method fabricates (MCTDHX, TD-2RDM) bind them; programs and
tests reach those fabricates only through `Method_Fabricate`, so no binding
of an earlier fabricate or of a test survives. `Method_Mb_OrbBased_Setup`
stops with a message if `Method_PropagateExactPart` is bound but
`Method_Mb_OrbBased_SetupFieldFree` is not (e.g. a test binding combined
with an orbital-based method).

### Exact Field-Free Part on grid.fem: Cayley Step

On grids without a dense h0 eigen-decomposition `exactFieldFreeQ` keeps its
meaning, h0 = T + V_static, and every stored orbital (operators of its body
type; restricted storage: body type 1) is propagated by the Crank–Nicolson
(Cayley) step

```
φ ← (1 + A)⁻¹ (1 − A) φ = 2 y − φ,   (1 + A) y = φ,   A = i Δt h0 / 2,
```

unitary in the grid metric for every Δt and time-reversible
(Cayley(−Δt) = Cayley(Δt)⁻¹); the phase per h0 eigenvalue E is 2 atan(E Δt/2)
instead of E Δt, i.e. second order like the Strang composition it is used in.
(1 + A) y = φ is solved by the flexible GMRES(m) of `M_Utils_Gmres`
(src/Utils/A_Utils.md) in the grid metric (`Grid_InnerProduct`),
right-preconditioned with the exact inverse of the kinetic factor,
P⁻¹ = (1 + i Δt T/2)⁻¹ = σ (T + σ)⁻¹, σ = 2/(i Δt) (`SysKinetic_SolveShifted`
with a complex shift, fast diagonalization on grid.fem), so the Krylov
iteration only sees (1 + A) P⁻¹ = 1 + (i Δt/2) V P⁻¹: one preconditioner solve
and one application of V per iteration. V_static is always applied through the
SysPotential contract (`SysPotential_MultiplyWithExternalPotential` on the
cached `FillExternalPotential(t = 0)`; the backend owns its representation, a
nodal product on grid.fem). Initial guess y0 = P⁻¹ φ; the preconditioned
Arnoldi vectors are kept (no extra solve for the update); modified
Gram–Schmidt, complex Givens rotations. Every convergence decision, the first
one included, is made on the true residual φ − (1 + A) y with h0 applied
explicitly (`SysKinetic_MultiplyWithKineticOp` + V_static): for y0 it is
computed before the first cycle (it equals −(i Δt/2) V y0 up to the round-off
of P P⁻¹ = 1, so an inexact preconditioner could not be accepted unchecked),
and every restart cycle ends with it. The solver is generic; the Cayley
submodule supplies the physics as its three callbacks (module procedures that
read the step data Δt and body type): the preconditioner P⁻¹, the Arnoldi
operator v + (i Δt/2)(V + U) z for z = P⁻¹ v, and (1 + A) y with h0 applied
explicitly for the true residuals. It keeps the controls, the GMRES context,
the non-convergence message with its error stop, and the solve counters.

| Key (`method.mb.orbBased.`) | Default | Meaning |
|---|---|---|
| `exactFieldFreeTol` | 1e-12 | relative residual ‖φ − (1 + A) y‖ / ‖φ‖ (grid metric) |
| `exactFieldFreeKrylovDim` | 30 | GMRES restart length m (2m + 4 work vectors of Grid_nPoints: 2m + 2 in the GMRES context, plus the solution and a scratch vector; module data allocated at the first propagation after the setup and kept; pages touched only for the iterations done) |
| `exactFieldFreeMaxRestarts` | 20 | restarts; not converged after 1 + maxRestarts cycles → error stop |
| `exactFieldFreeMeanFieldQ` | false | frozen mean-field reference: the step propagates h0 + U, U = J0 − K0 of the state at the first propagation (MCTDHX; next section) |

Setup (`Method_Mb_OrbBased_SetupFieldFreeCayley`) requires a time-independent
kinetic operator and a time- and body-type-independent potential (static h0;
V_static = `SysPotential_FillExternalPotential` at t = 0 is cached, zero
without `sysPotential`); a laser belongs to SysGauge (`sysGauge.fem.length`,
E(t) e_pol·r), which the linear integrator propagates (`ApplyLinOp` = gauge).
Diagnostics: `Method_Mb_OrbBased_cayleyNSolves`, `_cayleyNIterations`,
`_cayleyMaxIterations` (since setup). Because the GMRES context, the work
vectors and the step data are shared,
`Method_Mb_OrbBased_PropagateFieldFreeCayley` is not reentrant: it stops when
called inside an active OpenMP parallel region (the methods call it serially;
its vector passes and the kinetic solves are parallel inside). The orbital
norms drift by about the GMRES tolerance per solve (T_Fem_06: N − 2 = 4.9e-11
after 200 solves);
long runs can re-orthonormalize with
`Method_Mb_OrbBased_OrthonormalizeKeepingState` (MCTDHX, grid-agnostic).

Not available (clear error stops): imaginary time (Cayley does not damp the
high modes; use the SCF ground state), phi functions (`propagator.exponential`
stops at fabricate), `exactFieldFreeCut` (stops at Method fabricate on non-Ylm
grids), `exactFieldFreeGpuQ` (stops at setup). Use `propagator.splitStep` with `exactPartScheme: "strang"` (second
order; the fourth-order compositions gain nothing over the second-order Cayley
step). Each order-2 step makes 4 Cayley solves per orbital (exact(Δt/4)
twice per linear half step).

Tests (H2 on the mesh of T_Fem_03, N = 130331): T_Fem_05_H2FieldFree: the
RHF state stays stationary (|ΔE| 2e-10, N − 2 4e-12, 1 − |⟨φ0|φ⟩| 1e-12
over 0.5 a.u.) and the GMRES iteration counts stay within bounds (5 per
solve at Δt/4 = 0.0025 in the propagation, bound 7; 7 at Δt = 0.05, bound 9:
a wrong Krylov operator or initial residual still converges through the
true-residual restarts, only with more iterations); T_Fem_06/07: dipole kick
κ = 1e-3, 2e-3 (TRK initial slope, linearity; T_Fem_07 checks the run header
of T_Fem_06's `dipoleKick.dat`: κ/2, steps, mesh size, E_el(0),
completeness); T_Fem_08: `propagator.exponential` refused; T_Fem_09 and
T_Fem_10: the frozen mean-field reference in TDHF and in MCTDHF (next
section). Component tests (H atom, N = 12167, under 1 s each):
T_Method_FieldFreeCayleyFem runs the step with Krylov dimension 30 and with
Krylov dimension 2 and 50 restarts (restart path) at Δt = 0.0025 and 0.05
and checks the Crank–Nicolson residual, unitarity of the overlap matrix,
Cayley(−Δt) Cayley(Δt) = 1, the agreement of both (all ~1e-12, tolerance
1e-11) and the iteration counts; T_Method_FieldFreeCayleyImagTimeFem and
T_Method_FieldFreeCutFem check the refusals of imaginary time and of
`exactFieldFreeCut`. Timings (per FDM solve, Cayley solve, GMRES iteration
and time step, H2 and H2O) are kept in one place, app/fem_molecules/README.md.

### Frozen Mean-Field Reference of the Cayley Step (`exactFieldFreeMeanFieldQ`)

**Why.** The split above propagates h0 = T + V_nuc exactly and the mean field W
(Hartree, exchange and correlation terms of the orbital equations) with the
nonlinear integrator. In a heavy core both are large and cancel only in their
sum: for the O 1s of H2O ⟨φ|h0|φ⟩ = −33.1 Ha and ⟨φ|J − K|φ⟩ = +12.5 Ha
against ε = −20.6 Ha (valence orbitals: about −7 and +6.4 Ha against ε of
−0.5 to −1.3 Ha), so at the ground state the nonlinear part carries
−i Q W φ = i Q h0 φ of the size of the core energies and the commutators of
the Strang split are huge. Under this split the ground state of H2O is far
from stationary, and its TDHF energy drifts by several Ha per a.u. at
Δt = 0.05 (the measured drifts: app/fem_molecules/README.md, "H₂O
dynamics"); H2 (no core) was fine.

**What.** `"exactFieldFreeQ": true, "exactFieldFreeMeanFieldQ": true`
(MCTDHX with the Cayley part) splits instead

```
exact part      A = h0 + U       (Cayley/GMRES, orbitals only, unprojected)
nonlinear part  B = Ŵ − Û        (orbitals −i Q (W − U) φ, CI −i (H2 − Û) c)
U = J0 − K0     static, Hermitian in the grid metric
```

with U built once (`Method_Mb_OrbBased_BuildFieldFreeMeanField`): by default
from the state at the first exact-part call, i.e. after a program's kick
(called by the MCTDHX wrapper with `Method_Mb_OrbBased_FillRdm1` of the full
state; also at the first nonlinear call if a propagator started with that
part), unless the program built it before from another state ("Which
state" below). The factors are those of `ApplyHartreeFockOp`:
- J0: Hartree potential of the density ρ = Σ_kl ρS(k,l) φ_l* φ_k of the
  spatial (spin-summed) 1-RDM (closed shell: 2 Σ|φ|²), one potential solve
  of the total source (`SysInteraction_FillDirectPotential` with the
  weighted orbitals φ ρS), applied through `SysInteraction_MultiplyWithInteractionPotential`;
- K0: exchange operator of the 1-RDM within every body type (equal spins;
  restricted storage: the spin-up block for all shared orbitals),
  K0 v = Σ_jl rdm1(j,l) φ_j V[φ_l* v], compressed with ACE by
  `SysInteraction_FillExchangeAce` (shared with groundSolver.scf.stdImpl
  `aceQ`): K0 ≈ Σ_k |ζ_k⟩⟨ζ_k| with ζ = ξ L⁻ᴴ, ξ = K0 φ, ⟨φ|ξ⟩ = L Lᴴ,
  exact on span{φ}, Hermitian and positive semidefinite.

For a single determinant (TDHF = MCTDHX with nOrbs = nBodies) h0 + U is its
Fock operator, so the canonical SCF orbitals are eigenvectors of A and B
vanishes at the ground state; after a kick B = O(κ).

**Code.** U has its own submodule, `S_Method_Mb_OrbBased_FieldFreeMeanField.f90`,
which owns its data (J0, the ACE vectors and their ranges per body type, the
built flag) and its procedures. The mean-field algebra is SysInteraction's
(`SysInteraction_FillDirectPotential`, `SysInteraction_FillExchangeAce`, and
`SysInteraction_AddMeanField` for U v with the ACE vectors of the body type
of v; the SCF Fock action with `aceQ` uses the same helpers). The Cayley
step applies U per vector in its Arnoldi operator and true residual
(`Method_Mb_OrbBased_ApplyFieldFreeMeanFieldOnOrb`, a no-op until U is built
or without `exactFieldFreeMeanFieldQ`), MCTDHX per orbital together with its
matrix (`Method_Mb_OrbBased_ApplyFieldFreeMeanField`; the matrix by
`Method_Mb_OrbBased_FillOneBodyMatrix`, as h¹), and the Cayley setup ends with
`Method_Mb_OrbBased_SetupFieldFreeMeanField`, which releases U, checks the
interaction and prints the mean-field line.

**Which state.** For linear response build U from the unperturbed ground
state, before the kick, as app/fem_molecules does (`P_FemKickTdhf`:
`Method_Mb_OrbBased_BuildFieldFreeMeanField` with the SCF orbitals and their
1-RDM; the first exact-part call then keeps it). B(t) = W(t) − W0 is then
the pure response δW(t), O(κ) and oscillating. U from the kicked state (the
default) leaves the static part −δW(0+), the exchange change of the kicked
orbitals, in B at all times, and U from a later state of the run a static
O(κ) part of another shape. These static parts cost far more than the
oscillating response (H2O, presumably through the core): after a kick the
energy drifts about two orders of magnitude faster with U from the kicked
state than with U from the ground state, and with U rebuilt from a later
state of a kicked run it grows quadratically (app/fem_molecules/README.md,
"Which U").

**Equations of motion unchanged.** A at fixed CI coefficients is the exact
one-body flow of Ĥ0 + Û (unprojected orbitals); B is the flow of V̂_ee − Û in
the gauge without in-space orbital rotation: the orbital derivative subtracts
U φ before the projector (`Method_Mb_OrbBased_ApplyFieldFreeMeanField`), and
the CI derivative gets −Û c with Û = Σ_kj ⟨φ_k|U|φ_j⟩ a⁺_k a_j (one
`Coeffs_ApplyH1FillRdm1`). A + B is the full Hamiltonian, so the many-body
state follows the unchanged equations up to the splitting error. For
nOrbs = nBodies Û c is parallel to c and `Coeffs_ProjectOnSubspace` removes
it: the TDHF CI equation is untouched. For nOrbs > nBodies Û c is not
parallel to c (Û weighs configurations of different occupations
differently), and without the CI term the coefficients would not follow the
in-space rotation of the orbitals (T_Fem_10). The alternative with U outside
the projector also in B (orbitals −i Q W φ + i U φ, CI untouched) is the
same flow in the gauge with in-space rotation by −⟨φ|U|φ⟩: at the ground
state it leaves the rotation i P W φ = i φ w in the RK4 part (⟨1s|U|1s⟩ =
12.5 Ha for the O 1s of H2O, w Δt = 0.6 at Δt = 0.05), which RK4 does not
integrate unitarily, so B does not vanish there. Measured on H2O at
Δt = 0.05 (scratch build): an energy drift of 0.55 Ha/a.u., six orders of
magnitude above that of the same run with U inside the projector, and N
falling by 1.3e-3 per step instead of being conserved to 1e-11.

**Cost.** Build: 1 + n(n+1)/2 potential solves for n stored orbitals (H2: 2;
H2O, n = 5: 16). Per GMRES iteration and per true residual: one
multiplication by J0, n inner products ⟨ζ_k|z⟩ (parallel over k) and one
fused pass over the n ACE vectors (the preconditioner is unchanged); per
nonlinear derivative: U φ for the n orbitals and the n × n matrix ⟨φ|Uφ⟩,
negligible against the pair potentials of W. Memory: n ACE vectors and J0.
Timings: app/fem_molecules/README.md.

**Tests.** T_Fem_09_H2KickMeanFieldReference (H2, N = 130331): at the
ground state |B φ| / |−i Q W φ| = 2.5e-11 (the 1-RDM regularization
n/(n² + 1e-10) at n = 2); the ground state at Δt = 0.1 stays stationary
(|ΔE| 6e-13, 1 − |⟨φ0|φ⟩| 2e-13, 6 GMRES iterations per solve); after the
kick κ = 1e-3 at Δt = 0.02 the dipole agrees with the plain split (T_Fem_06)
to 2.7e-7 of the response over 0.6 a.u., and both splittings converge at
second order (p = 2.00 from Δt, Δt/2, Δt/4) to the same trajectory, with
errors −5.7e-6 (reference) and −6.4e-6 (plain) at t = 0.1; a reference
built before the kick from the ground state is kept at the first call (the
nonlinear derivative at the kicked state is then the response, 3.0e-3
relative, against 2.5e-11 with the default reference from the kicked
state).
T_Fem_10_H2MctdhfMeanFieldReference (MCTDHF: H2 with nOrbs = [2, 2], N =
15341, κ = 0.1 on |σg σg⟩, Δt = 0.005, 0.25 a.u.): the split with the
reference follows the plain split to 1.6e-7 of the dipole response and
3.4e-7 Ha in energy, while the part of the CI term that the projector keeps,
|(1 − |c⟩⟨c|) Û c|, reaches 4.5e-2 Ha; a scratch build without the CI term
leaves the trajectory by 4.5e-3 (dipole) and 1.7e-3 Ha (energy). Scratch
runs of the same system over 0.5 a.u. against the plain split at Δt =
0.00125: dipole within 2.0e-7 with and 3.4e-7 without the reference at
Δt = 0.005, 5.8e-6 and 9.5e-6 at Δt = 0.025.

The runs of app/fem_molecules are documented in one place,
app/fem_molecules/README.md ("Results on this machine"): H2 over 50 a.u. with
and without the reference, and H2O on the dynamics mesh (energy drift per
time step and choice of U, the ground state without kick, GMRES iterations
per solve, cost per GMRES iteration and per step, the 10 a.u. demo and its
continuation).

**Restrictions** (clear error stops): requires `exactFieldFreeQ`, a non-Ylm
grid (the Cayley part; Ylm grids propagate h0 by its eigen-decomposition),
`method.mb.orbBased.mctdhx`, `sysInteraction`, fermions (fabricate) and a
body-type-independent interaction (setup); imaginary time stops in the
Cayley step as without the reference. T_Method_FieldFreeMeanFieldTd2rdmFem
checks the refusal for TD-2RDM, whose nonlinear part does not subtract U.

**Lifetime, restarts, absorber.** Every `Method_Mb_OrbBased_SetupFieldFree`
(Method_Setup; the Cayley setup ends with
`Method_Mb_OrbBased_SetupFieldFreeMeanField`) releases the reference; the next
propagation builds it again unless the program does. Method neither saves nor
loads U, so a run continued from a checkpoint must rebuild U from the state
it was first built from, not from the checkpoint state ("Which state"); the
program owns this plumbing. `P_FemKickTdhf` does it with two of its own
procedures: `BuildAndSaveReference` builds U from the ground state before
the kick and saves that state (orbitals and 1-RDM, `meanFieldReference.dat`),
and in a continued run `LoadAndBuildReference` rebuilds U from it (it stops
if the file is missing, unreadable or from another mesh or orbital space).
The continued run then reproduces the uninterrupted one bit for bit (H2,
Δt = 0.05, restart at t = 0.2: `dipole.dat`, `observables.dat` and the
orbitals at t = 0.4 identical). A second program with restarts and
`exactFieldFreeMeanFieldQ` would be the reason to move this save and rebuild
into Method. An absorber (`Absorber_ApplyAbsorber`, with
`keepOrthonormalQ` `Method_Mb_OrbBased_OrthonormalizeKeepingState`) acts on
the state only; U stays as built (any static Hermitian U keeps the equations
exact), and both parts keep the orbitals orthonormal (A unitary, B
orthogonal to the orbital space). H2 with `absorber.fem.cosinus` and
`keepOrthonormalQ` at Δt = 0.05: N and the dipole agree with the plain split
to 1e-11 and 4e-6 (relative) over 1 a.u.

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
- **CI sector** (optional `coeffs.generic.totalM`): the CI Hamiltonian acts
  only on the configurations of the declared total M, with full-length
  storage (`src/Coeffs/A_Coeffs.md`); the initial coefficients are projected
  on it once (`Coeffs_ProjectOnSector` in `CoeffsInit_Initialize`). Ne 9
  orbitals: 4116 of 15876 configurations, −0.24 s per step.

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
`FillH2YlmNodes` forms only the h2 columns. The columns of the stored pairs
come first in every group (`T_TargetGroup % nColStored`) and are the only
ones entering h2: the swapped pairs follow from the Hermitian symmetry
h2(k, k2, i, l2) = conj(h2(i, l2, k, k2)) (`FillH2FromTargetColumns`; row k
lies in the group of the swapped column). This needs a Hermitian metric, the
precondition of the node path; a c-product (ECS) variant must compute the
swapped columns instead.

One MCTDHX non-linear derivative of Ne with 9 orbitals (1s 2s 2p 3s 3p,
lmax 40, mmax 1, 761 radial points, four threads on 2 cores, measured before
the parity split and the stored-pair h2 above): 0.38 s (0.28 s correlation
operator), against 3.7 s with the generic operator. One energy
evaluation (`Method_GetEnergy`, h2 by `FillH2YlmNodes`): 0.25 s instead of
2.2 s. The parity-split node transforms (A_Grid.md) and h2 from the stored pairs
only take 16 time steps of the Ne HHG setup (incl. setup and 3 energy
evaluations) from 48.7 s to 45.4 s.

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
| `Mb/OrbBased/S_Method_Mb_OrbBased_FieldFreeEigen.f90` | Exact field-free one-body propagation (Ylm, `*Eigen`) |
| `Mb/OrbBased/S_Method_Mb_OrbBased_FieldFreeCayley.f90` | Exact field-free part by the Cayley step (grid.fem, `*Cayley`): controls, the GMRES context and the three callbacks of `Gmres_Solve` (`M_Utils_Gmres`) |
| `Mb/OrbBased/S_Method_Mb_OrbBased_FieldFreeMeanField.f90` | Frozen mean-field reference U = J0 − K0 of the Cayley step (`exactFieldFreeMeanFieldQ`): setup, build, per-vector and per-orbital application |
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
