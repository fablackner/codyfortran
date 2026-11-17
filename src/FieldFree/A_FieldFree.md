# FieldFree Module

## Purpose

FieldFree propagates orbitals exactly with the field-free one-body
Hamiltonian h0 = T + V_ext. The orbital-based methods use it as the exactly
solved part of their split-step and exponential propagation: with a
`fieldFree` block in the input the linear integrator of the method only
propagates the gauge coupling, and the propagator composes every linear step
with `FieldFree_Propagate` (`propagator.splitStep.exactPartScheme`). Without
the block the linear integrator propagates the whole one-body part.

| | `fieldFree.eigen` | `fieldFree.cayley` |
|---|---|---|
| Grid | `grid.ylm` | grids whose kinetic backend binds `SysKinetic_SolveShifted` (grid.fem) |
| Step | exp(−i h0 Δt) from the eigen-decomposition of the h0 blocks per l | (1 + i Δt h0/2)⁻¹ (1 − i Δt h0/2) by preconditioned GMRES |
| Error in Δt | none | second order |
| Imaginary time | yes | no |
| Phi functions (`propagator.exponential`) | yes | no |
| Spectral cut | `cut` | no |
| Frozen mean-field reference | no | `meanFieldQ` |
| GPU | `gpuQ` (cuBLAS) | no |

## Architecture Overview

```
src/FieldFree/
├── M_FieldFree.f90                       # Contract: FieldFree_Fabricate, pointers, flags
├── S_FieldFree.f90                       # Dispatch: JSON → variant
├── Eigen/
│   ├── M_FieldFree_Eigen.f90             # fieldFree.eigen: cut, gpuQ
│   └── S_FieldFree_Eigen.f90             # h0 blocks per l, propagation, phi functions, cut
└── Cayley/
    ├── M_FieldFree_Cayley.f90            # fieldFree.cayley: GMRES counters
    ├── S_FieldFree_Cayley.f90            # Cayley step and its GMRES callbacks
    └── S_FieldFree_Cayley_MeanField.f90  # frozen mean-field reference U (nested submodule)
```

Programs fabricate and set up FieldFree next to the other modules:

```fortran
call SysInteraction_Fabricate
call FieldFree_Fabricate      ! before Method_Fabricate
call Method_Fabricate
...
call CoeffsInit_Setup
call FieldFree_Setup          ! after the setups of Grid, SysKinetic and SysPotential
call Method_Setup
```

`FieldFree_Fabricate` resets the pointers and flags of the module and binds
the variant of the block. For an input without a `fieldFree` block it binds
nothing and prints `de: FieldFree => noFieldFree`, as Absorber does. The
methods read a bound `FieldFree_Propagate` at their fabricate, so
`FieldFree_Fabricate` comes first. `Method_Mb_OrbBased_Fabricate` stops if
the input has a `fieldFree` block that `FieldFree_Fabricate` did not read, or
if FieldFree still holds the exact part of another input. A program without
the two calls therefore runs every input that has no `fieldFree` block.

## Key Interfaces

| Name | Purpose |
|---|---|
| `FieldFree_Fabricate` | Resets the module, reads `fieldFree`, dispatches to the variant |
| `FieldFree_Setup` | Builds the operators of the variant (eigen: h0 blocks per l; cayley: GMRES context, cached V_static), resets the GMRES counters and releases the mean-field reference; a no-op without the block |
| `FieldFree_Propagate(orbs, dt, imagTimeQ)` | orbs ← exp(−i h0 dt) orbs, or exp(−h0 dt) orbs with `imagTimeQ` (eigen); the Cayley step (cayley, real time only). Null without the block |
| `FieldFree_ApplyPhi(orbs, orbsList, dt)` | orbs = Σ_m φ_m(−i h0 dt) orbsList(:, :, m) for exponential integrators; eigen only, null otherwise |
| `FieldFree_ProjectBelowCut(orbs)` | Projects orbitals or orbital derivatives onto the h0 eigenvectors below `fieldFree.eigen.cut`; a no-op without a cut |
| `FieldFree_cutQ` | Whether a cut is set |
| `FieldFree_meanFieldQ` | `fieldFree.cayley.meanFieldQ`: the exact part propagates h0 + U, the method subtracts U in its nonlinear part |
| `FieldFree_meanFieldBuiltQ` | Whether U is built since the last `FieldFree_Setup` |
| `FieldFree_BuildMeanField(orbs, rdm1, time)` | Builds U from the stored orbitals and the 1-RDM of a state; cayley only, null otherwise |
| `FieldFree_ApplyMeanField(dOrbs, orbs, h1_)` | dOrbs(:, i) += U orbs(:, i), optionally the matrix ⟨φ_i\|U φ_j⟩; cayley only, null otherwise |
| `FieldFree_Eigen_cut`, `FieldFree_Eigen_gpuQ` | Values of the eigen keys |
| `FieldFree_Cayley_nSolves`, `FieldFree_Cayley_nIterations`, `FieldFree_Cayley_maxIterations` | GMRES solves, their total and largest iteration counts since the setup |

## JSON Configuration

```json
"fieldFree": {
    "eigen": {
        "cut": 1000.0,
        "gpuQ": false
    }
}
```

```json
"fieldFree": {
    "cayley": {
        "tolerance": 1e-12,
        "krylovDim": 30,
        "maxRestarts": 20,
        "meanFieldQ": true
    }
}
```

Every key has a default, so `"eigen": {}` and `"cayley": {}` are complete.
All keys are read at fabricate.

The block replaces seven keys of `method.mb.orbBased`. An input that still
has one of them stops at `Method_Fabricate` with a message that names the
block:

| Former key (`method.mb.orbBased.`) | Now |
|---|---|
| `exactFieldFreeQ: true` | a `fieldFree` block with `eigen` (grid.ylm) or `cayley` (grid.fem) |
| `exactFieldFreeCut` | `fieldFree.eigen.cut` |
| `exactFieldFreeGpuQ` | `fieldFree.eigen.gpuQ` |
| `exactFieldFreeTol` | `fieldFree.cayley.tolerance` |
| `exactFieldFreeKrylovDim` | `fieldFree.cayley.krylovDim` |
| `exactFieldFreeMaxRestarts` | `fieldFree.cayley.maxRestarts` |
| `exactFieldFreeMeanFieldQ` | `fieldFree.cayley.meanFieldQ` |

## Use by the Methods

With a `fieldFree` block the linear step of the split-step propagation
(src/Method/A_Method.md, "Split-Step Propagation") is split once more: the
field-free h0 = T + V_ext is propagated exactly and the linear integrator
only sees the gauge coupling (`Method_Mb_OrbBased_ApplyLinOp` binds the gauge
operator instead of the full one-body operator). The methods bind
`Method_PropagateExactPart` to a wrapper of `FieldFree_Propagate`, and the
split-step propagators compose every linear step from exact and integrator-1
substeps (`propagator.splitStep.exactPartScheme`; the default `strang` is
exp(−i h0 Δt/2) · integrator 1(Δt) · exp(−i h0 Δt/2), local error O(Δt³), as
the outer Strang splitting).

A method declares at its fabricate what it implements, and
`Method_Mb_OrbBased_Fabricate` stops for a requested option that the selected
method does not declare:

| Method | Exact part (binds `Method_PropagateExactPart`) | Spectral cut (`Method_Mb_OrbBased_implementsSpectralCutQ`) | Mean-field reference (`Method_Mb_OrbBased_implementsMeanFieldReferenceQ`) |
|---|---|---|---|
| MCTDHX | yes | yes | yes |
| TDKS | yes | yes | yes |
| TD-2RDM | yes, without weak orbitals | no | no |
| TDCI, TDHX, TD-2RDM-StaticOrbs | no | no | no |

`Method_PropagateExactPart` and `Method_ApplyExactPartPhi` are owned by
`Method_Fabricate`, which nullifies them before the method fabricates bind
them. MCTDHX and TDKS bind `Method_ApplyExactPartPhi` only if
`FieldFree_ApplyPhi` is associated (eigen), so `propagator.exponential` stops
at fabricate with the Cayley variant.

## Eigen Variant (`fieldFree.eigen`)

`S_FieldFree_Eigen.f90` builds h0 once in `FieldFree_Setup`: it applies
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
dominates) at 2.5–5× lower cost (app/atoms/be; measurements in the git history of app/be/README.md).

`"gpuQ": true` runs the two GEMMs per l on the GPU with cuBLAS
(`M_Utils_CublasLib`, build with `cmake -DCODY_USE_CUDA=ON`): the eigenvectors
stay on the device, each call gathers all blocks on the host, sends them once,
applies the phases on the host between the two transformations and returns
them. Results agree with the CPU path to round-off.

**Spectral cut.** `"cut": Ec` (MCTDHX and TDKS, default 0 = off)
restricts the orbitals to the h0 eigenvectors with E < Ec, P = Σ_{E<Ec} |u⟩⟨u|.
The highest FEDVR eigenvalues (~10⁵ a.u., from the dense Lobatto points at
the element boundaries) are unphysical, but they make the Strang composition
of h0 and the gauge coupling inaccurate. With the cut the exact part
transforms with the kept eigenvectors only (it projects onto P), and
`FieldFree_ProjectBelowCut` applies P to the linear and nonlinear
orbital derivatives of MCTDHX and TDKS (and to the full TDKS derivative), so
P is invariant under the whole step. `FieldFree_Eigen_Fabricate` binds the
pointer to the projection (two GEMMs with the few dropped eigenvectors,
P = 1 − U_H U_Hᵀ); without a cut it stays a no-op. The cut changes the model:
relax the ground state with the same cut (the initial projection of an uncut
state is a jump). For Be HHG, Ec = 1000 with `strang` is ~2× faster per step
than `blanesMoan6Swap` and shifts harmonics 21–60 by 0.02–0.04 decades
(git history of app/be/README.md). T_Ne3d_11 checks that no weight appears above the cut.

**Phi functions.** MCTDHX and TDKS with `fieldFree.eigen` also bind
`Method_ApplyExactPartPhi(res, vecs, dt)`: res = Σ_m φ_m(−i H0 dt) vecs(:, m+1)
with φ_0(z) = e^z and φ_m(z) = (φ_{m−1}(z) − 1/(m−1)!)/z, H0 = h0 on the
orbitals and 0 on the CI coefficients (φ_m(0) = 1/m!). The orbital part,
`FieldFree_ApplyPhi`, transforms every vector into the h0
eigenbasis (kept eigenvectors only), weights it with φ_m(−i E dt) (Taylor
series for |z| < 1, recurrence otherwise) and transforms the sum back once
(CPU only). It is consumed by the exponential integrators
(`propagator.exponential`, src/Propagator/A_Propagator.md); real time only.

Requirements: `grid.ylm` (checked at fabricate); a Hermitian (non-ECS) grid,
time- and body-type-independent kinetic and potential and a central
potential (`SysPotential_Ylm_lmax = 0`; checked at setup); no TD-2RDM weak
orbitals (Method fabricate). Memory: one real
nRadial² matrix per l. T_Ne3d_09 (Strang) and T_Ne3d_10 (Blanes–Moan) check
the agreement with the SIL scheme of T_Ne3d_08.

## Cayley Variant (`fieldFree.cayley`)

`fieldFree.cayley` is the variant for grids without a dense h0
eigen-decomposition: h0 = T + V_static, and every stored orbital (operators
of its body type; restricted storage: body type 1) is propagated by the
Crank–Nicolson (Cayley) step

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

| Key (`fieldFree.cayley.`) | Default | Meaning |
|---|---|---|
| `tolerance` | 1e-12 | relative residual ‖φ − (1 + A) y‖ / ‖φ‖ (grid metric) |
| `krylovDim` | 30 | GMRES restart length m (2m + 4 work vectors of Grid_nPoints: 2m + 2 in the GMRES context, plus the solution and a scratch vector; module data allocated at the first propagation after the setup and kept; pages touched only for the iterations done) |
| `maxRestarts` | 20 | restarts; not converged after 1 + maxRestarts cycles → error stop |
| `meanFieldQ` | false | frozen mean-field reference: the step propagates h0 + U, U = J0 − K0 of the state at the first propagation (MCTDHX; U = v_KS[ρ0] − a K0 for TDKS; next section) |

The keys are read at fabricate. `FieldFree_Setup` requires a time-independent
kinetic operator and a time- and body-type-independent potential (static h0;
V_static = `SysPotential_FillExternalPotential` at t = 0 is cached, zero
without `sysPotential`); a laser belongs to SysGauge (`sysGauge.fem.length`,
E(t) e_pol·r), which the linear integrator propagates (`ApplyLinOp` = gauge).
Diagnostics: `FieldFree_Cayley_nSolves`, `FieldFree_Cayley_nIterations`,
`FieldFree_Cayley_maxIterations` (since setup). Because the GMRES context,
the work vectors and the step data are shared, `FieldFree_Propagate` of this
variant is not reentrant: it stops when
called inside an active OpenMP parallel region (the methods call it serially;
its vector passes and the kinetic solves are parallel inside). The orbital
norms drift by about the GMRES tolerance per solve (T_Fem_06: N − 2 = 4.9e-11
after 200 solves);
long runs can re-orthonormalize with
`Method_Mb_OrbBased_OrthonormalizeKeepingState` (MCTDHX, grid-agnostic).

Not available: imaginary time (error stop; Cayley does not damp the high
modes; use the SCF ground state), phi functions (`propagator.exponential`
stops at fabricate), the spectral cut and the GPU path (keys of
`fieldFree.eigen` only). Use `propagator.splitStep` with `exactPartScheme: "strang"` (second
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
T_FieldFree_CayleyFem runs the step with Krylov dimension 30 and with
Krylov dimension 2 and 50 restarts (restart path) at Δt = 0.0025 and 0.05
and checks the Crank–Nicolson residual, unitarity of the overlap matrix,
Cayley(−Δt) Cayley(Δt) = 1, the agreement of both (all ~1e-12, tolerance
1e-11) and the iteration counts; T_FieldFree_CayleyImagTimeFem and
T_FieldFree_EigenFem check the refusals of imaginary time and of
`fieldFree.eigen` on grid.fem. Timings (per FDM solve, Cayley solve, GMRES iteration
and time step, H2 and H2O) of the app/molecules runs are in the
git history of app/fem_molecules/README.md.

## Frozen Mean-Field Reference (`fieldFree.cayley.meanFieldQ`)

**Why.** The split above propagates h0 = T + V_nuc exactly and the mean field W
(Hartree, exchange and correlation terms of the orbital equations) with the
nonlinear integrator. In a heavy core both are large and cancel only in their
sum: for the O 1s of H2O ⟨φ|h0|φ⟩ = −33.1 Ha and ⟨φ|J − K|φ⟩ = +12.5 Ha
against ε = −20.6 Ha (valence orbitals: about −7 and +6.4 Ha against ε of
−0.5 to −1.3 Ha), so at the ground state the nonlinear part carries
−i Q W φ = i Q h0 φ of the size of the core energies and the commutators of
the Strang split are huge. Under this split the ground state of H2O is far
from stationary, and its TDHF energy drifts by several Ha per a.u. at
Δt = 0.05 (the measured drifts: git history of app/fem_molecules/README.md, "H₂O
dynamics"); H2 (no core) was fine.

**What.** `"fieldFree": {"cayley": {"meanFieldQ": true}}` (MCTDHX and TDKS)
splits instead

```
exact part      A = h0 + U       (Cayley/GMRES, orbitals only, unprojected)
nonlinear part  B = Ŵ − Û        (orbitals −i Q (W − U) φ, CI −i (H2 − Û) c)
U = J0 − K0     static, Hermitian in the grid metric
```

with U built once (`FieldFree_BuildMeanField`): by default
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

**Code.** U has its own submodule, `S_FieldFree_Cayley_MeanField.f90`
(a nested submodule of `S_FieldFree_Cayley`: an add-on of the
Cayley step), which owns its data (J0, the ACE vectors and their ranges per
body type, the built flag) and its procedures; the mean field as a value
(`Method_Mb_OrbBased_FillMeanField`, shared with the adiabatic switching) is
`S_Method_Mb_OrbBased_MeanField.f90`. The mean-field algebra is SysInteraction's
(`SysInteraction_FillDirectPotential`, `SysInteraction_FillExchangeAce`, and
`SysInteraction_AddMeanField` for U v with the ACE vectors of the body type
of v; the SCF Fock action with `aceQ` uses the same helpers). The Cayley
step applies U per vector in its Arnoldi operator and true residual
(`ApplyMeanFieldOnOrb`, a no-op until U is built
or without `meanFieldQ`), MCTDHX per orbital together with its
matrix (`FieldFree_ApplyMeanField`; the matrix by
`Method_Mb_OrbBased_FillOneBodyMatrix`, as h¹), and the Cayley setup ends with
`SetupMeanField`, which releases U, checks the bodies and the
interaction and prints the mean-field line.

**Which state.** For linear response build U from the unperturbed ground
state, before the kick, as app/molecules does (`P_KickTdhfFem`:
`FieldFree_BuildMeanField` with the SCF orbitals and their
1-RDM; the first exact-part call then keeps it). B(t) = W(t) − W0 is then
the pure response δW(t), O(κ) and oscillating. U from the kicked state (the
default) leaves the static part −δW(0+), the exchange change of the kicked
orbitals, in B at all times, and U from a later state of the run a static
O(κ) part of another shape. These static parts cost far more than the
oscillating response (H2O, presumably through the core): after a kick the
energy drifts about two orders of magnitude faster with U from the kicked
state than with U from the ground state, and with U rebuilt from a later
state of a kicked run it grows quadratically (git history of app/fem_molecules/README.md,
"Which U").

**Equations of motion unchanged.** A at fixed CI coefficients is the exact
one-body flow of Ĥ0 + Û (unprojected orbitals); B is the flow of V̂_ee − Û in
the gauge without in-space orbital rotation: the orbital derivative subtracts
U φ before the projector (`FieldFree_ApplyMeanField`), and
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
Timings: git history of app/fem_molecules/README.md.

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

The runs of app/molecules are documented in the
git history of app/fem_molecules/README.md ("Results on this machine"): H2 over 50 a.u. with
and without the reference, and H2O on the dynamics mesh (energy drift per
time step and choice of U, the ground state without kick, GMRES iterations
per solve, cost per GMRES iteration and per step, the 10 a.u. demo and its
continuation).

**Restrictions** (clear error stops): `sysInteraction` or `sysXc`
(fabricate), `method.mb.orbBased.mctdhx` or `method.mb.orbBased.tdks` (Method
fabricate), fermions and a body-type-independent interaction (setup, the
interaction if `sysInteraction` exists); imaginary time stops in the Cayley
step as without the reference. The key belongs to `fieldFree.cayley` alone
(Ylm grids propagate h0 by its eigen-decomposition).
test/componentTests/Method/T_Method_FieldFreeMeanFieldTd2rdmFem checks the
refusal for TD-2RDM (in Method_Mb_OrbBased_Fabricate), whose nonlinear part
does not subtract U.

**Kohn–Sham reference (TDKS).** With `method.mb.orbBased.tdks`
(`Method_kohnShamQ`, set by its fabricate) U is the Kohn–Sham
mean field of the reference state, U = v_KS,s[ρ0] − a K0: the local
potential per spin channel on the xc points of SysXc
(`Method_Mb_OrbBased_Tdks_FillKohnShamPotential`, applied by
`SysXc_MultiplyWithPotential`) and, for a hybrid, the ACE vectors of K0 scaled
by √a (none for a = 0: the build costs one Poisson solve and one functional
evaluation). The TDKS nonlinear part is −i (v_KS − a K − U) φ, which vanishes
at the Kohn–Sham ground state. A phase kick leaves the density unchanged, so
for semilocal functionals the default reference of the kicked state is the
ground-state one. H2 LDA (T_FemDft_03, Δt = 0.02): the ground state is
stationary to 3e-11 over 1 a.u. (5e-6 with the plain split h0 | v_KS, whose
secular energy error is not monotonic in Δt: 4.9e-6, 4e-7 and 3.5e-6 per a.u.
at Δt = 0.04, 0.02, 0.01), and after a kick the Kohn–Sham energy is conserved
to 2e-11 over 0.6 a.u. with the dipole of the plain split at Δt = 0.005 to
5e-6 of the response.

**Lifetime, restarts, absorber.** Every `FieldFree_Setup` (the Cayley
setup ends with `SetupMeanField`) releases the reference; the next
propagation builds it again unless the program does. FieldFree neither saves
nor loads U, so a run continued from a checkpoint must rebuild U from the state
it was first built from, not from the checkpoint state ("Which state"); the
program owns this plumbing. `P_KickTdhfFem` does it with two of its own
procedures: `BuildAndSaveReference` builds U from the ground state before
the kick and saves that state (orbitals and 1-RDM, `meanFieldReference.dat`),
and in a continued run `LoadAndBuildReference` rebuilds U from it (it stops
if the file is missing, unreadable or from another mesh or orbital space).
The continued run then reproduces the uninterrupted one bit for bit (H2,
Δt = 0.05, restart at t = 0.2: `dipole.dat`, `observables.dat` and the
orbitals at t = 0.4 identical). A second program with restarts and
`meanFieldQ` would be the reason to move this save and rebuild
into FieldFree. An absorber (`Absorber_ApplyAbsorber`, with
`keepOrthonormalQ` `Method_Mb_OrbBased_OrthonormalizeKeepingState`) acts on
the state only; U stays as built (any static Hermitian U keeps the equations
exact), and both parts keep the orbitals orthonormal (A unitary, B
orthogonal to the orbital space). H2 with `absorber.fem.cosinus` and
`keepOrthonormalQ` at Δt = 0.05: N and the dipole agree with the plain split
to 1e-11 and 4e-6 (relative) over 1 a.u.

## Dependencies

FieldFree applies T and V through the contracts of SysKinetic and
SysPotential and takes the inner product from Grid. From the method layer it
reads the layout of the stored orbitals: the body type of every orbital
(`Method_Mb_OrbBased_bodyTypeOfOrb`, cayley) and the declared m
(`Method_Mb_OrbBased_mOfOrb`, eigen). The frozen mean-field reference is
built with the mean field of the method layer
(`Method_Mb_OrbBased_FillMeanField`, `Method_Mb_OrbBased_AddMeanField`; for a
Kohn–Sham method `Method_Mb_OrbBased_Tdks_FillKohnShamPotential` and SysXc).
The interface modules `M_FieldFree*` import only `M_Utils_*`.

## Testing

Component tests in `test/componentTests/FieldFree/` (H atom, under 1 s
each; grid.fem with N = 12167 unless noted):

| Test | Checks |
|---|---|
| `T_FieldFree_EigenYlm` | The eigen variant on a Ylm FEDVR grid (lmax = mmax = 2): unitarity, exp(−i h0 dt)² = exp(−2i h0 dt), the backward step, φ_0 = the propagation; the lowest eigenvector of every l block from imaginary time is an eigenvector of T + V applied through SysKinetic and SysPotential and gets the phase exp(−i E dt) and φ_0 + φ_1 = e^z + (e^z − 1)/z; the cut is a projector, and the propagation stays in its range |
| `T_FieldFree_CayleyFem` | The Cayley step on its own with one GMRES cycle and with restarts: Crank–Nicolson residual, unitarity, Cayley(−Δt) Cayley(Δt) = 1, iteration counts |
| `T_FieldFree_CayleyImagTimeFem` | Expected error: imaginary time with the Cayley step |
| `T_FieldFree_EigenFem` | Expected error: `fieldFree.eigen` on grid.fem |

Simulation tests: T_Ne3d_09 to T_Ne3d_16 (eigen with MCTDHX: Strang,
Blanes–Moan, spectral cut, exponential integrators, declared m), T_He3d_04
(eigen with TDKS, a GGA and a cut), T_Be3d_03 (eigen with TD-2RDM), T_Fem_05
to T_Fem_10 and T_FemDft_03 (cayley with and without the mean-field
reference).

## Common Pitfalls

- A program that does not call `FieldFree_Fabricate` stops at
  `Method_Fabricate` for an input with a `fieldFree` block.
- The keys are read at fabricate. After changing a GMRES control in the JSON
  call `FieldFree_Fabricate` and `FieldFree_Setup` again
  (T_FieldFree_CayleyFem).
- `FieldFree_Setup` releases the frozen mean-field reference and resets the
  GMRES counters. A program that built U before a new setup must build U
  again.
- A spectral cut changes the model: relax the ground state with the same cut.
- The Cayley step is second order in Δt, so the fourth-order compositions of
  `exactPartScheme` gain nothing over `strang`.
- The Cayley step is not reentrant and real-time only.
