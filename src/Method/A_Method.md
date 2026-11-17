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
       ┌──────────┬───────┴───────┬────────────┬─────────────────┬──────────┐
       ▼          ▼               ▼            ▼                 ▼          ▼
    ┌──────┐  ┌──────┐       ┌────────┐   ┌────────┐    ┌─────────────────┐ ┌──────┐
    │ TDCI │  │ TDHX │       │ MCTDHX │   │TD-2RDM │    │TD-2RDM-Static   │ │ TDKS │
    │      │  │      │       │        │   │        │    │     Orbs        │ │(DFT) │
    └──────┘  └──────┘       └────────┘   └────────┘    └─────────────────┘ └──────┘
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
| TDKS                | Kohn–Sham orbitals φ_i(r) only (single determinant) |
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
| `Method_ApplyHamiltonian` | H·|Ψ⟩ = i × `Method_TimeDerivative`: the `ApplyMatOnVec` of diagonalizers |
| `Method_ImagTimeDerivative` | −H·|Ψ⟩ = −i × `Method_TimeDerivative`: the `TimeDerivative` of imaginary-time integrators |
| `Method_kohnShamQ`        | Set by a Kohn–Sham method; `Method_Fabricate` refuses a `sysXc` block otherwise |
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
| `Method_Mb_OrbBased_FillStateWeights` (module procedure) | Diagonal metric of `Method_state` for the `metricWeights` of integrators and diagonalizers: 1 for the many-body part, `Grid_weights` for every orbital of the state (moving orbitals) |
| `Method_Mb_OrbBased_FillInverseRdm1` (module procedure) | Regularized spatial inverse U n/(n²+ε) U† of the 1-RDM |
| `Method_Mb_OrbBased_ApplyCiHamiltonian` (module procedure) | dCoeffs = H coeffs with the stored `Method_Mb_OrbBased_h1`/`h2` (no refill); the CI matvec of TDCI, of every `groundSolver.mcscf` backend (`GroundSolver_Mcscf_HamiltonianAction`) and of CI start vectors in programs |
| `Method_Mb_OrbBased_SaveRdm2Blocks` (module procedure) | Write body-type blocks `rdm2BtXX_YY.in` of a 2-RDM (read by `twoRdmInit.load`) |
| `Method_Mb_OrbBased_H1TimeIndependentQ` (module procedure) | Whether h¹ = T̂ + V̂_ext (+ Ĥ_A if `sysGauge` is configured) is time-independent (the `*_timeIndependentQ` flags at call time) and there is no adiabatic switching; TDCI and TD-2RDM-StaticOrbs refill h¹ at every derivative otherwise |
| `Method_Mb_OrbBased_FillH0`, `Method_Mb_OrbBased_RefillStaticH1` (module procedures) | h⁰ = ⟨T̂ + V̂_ext⟩ of the orbitals (the static h¹ of the static-orbital setups); the h¹ refill of the static-orbital methods: nothing if h¹ is time-independent, the static part (h¹ at the first call after `Method_Mb_OrbBased_Setup`: h⁰, or what a program stored) plus the gauge coupling if T and V are static, `Method_Mb_OrbBased_FillH1` otherwise, then the switching reference (N2 TD-2RDM in a pulse, one thread: the refill drops from about 0.05 to 0.02 s per derivative) |
| `Method_Mb_OrbBased_FillOneBodyMatrix` (module procedure) | h¹(i,j) = ⟨φ_i\|hφ_j⟩ of a one-body operator from its action on the stored orbitals (equal body types, nOrbsSum layout, restricted: both spin blocks); the matrix of `FillH1` and of the mean-field reference |
| `T_Method_Mb_OrbBased_MeanField`, `Method_Mb_OrbBased_FillMeanField`, `_AddMeanField` (type and module procedures, `S_Method_Mb_OrbBased_MeanField`) | Frozen Hartree–Fock mean field U = J0 − K0 of a 1-RDM on the stored orbitals as a value (J0, ACE vectors of K0 per body type) and w += U z; the reference of the Cayley step (FieldFree) and of the adiabatic switching |
| `Method_Mb_OrbBased_SwitchFactor`, `_SetSwitchingReference`, `_AddSwitchingReference`, `_AddSwitchingReferenceMatrix`, `Method_Mb_OrbBased_switchTime` (`S_Method_Mb_OrbBased_Switching`) | Adiabatic switching (`adiabaticSwitching`, below): λ(t), the reference U0 of the current state (grid for moving, matrix of h2 for static orbitals; at the end of the method setups, again by programs that change the state), dOrbs += (1 − λ) U0 orbs, h1 += (1 − λ) U0 |
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
| **TDKS**         | Dynamic    | Implicit     | Time-dependent Kohn–Sham (real-time TDDFT) with a functional of `sysXc`; see below |

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

### Exact Field-Free Part (`fieldFree`)

The exact propagation of the field-free h0 = T + V_ext is the module
FieldFree with its own JSON block `fieldFree` (src/FieldFree/A_FieldFree.md:
the eigen variant of Ylm grids with the spectral cut and the phi functions,
the Cayley step of grid.fem, the frozen mean-field reference). Method keeps
the method side:

- `Method_Mb_OrbBased_Fabricate` binds `Method_Mb_OrbBased_ApplyLinOp` to the
  gauge coupling (plus the switching reference) if `FieldFree_Propagate` is
  bound, and to the full one-body operator otherwise. It stops if the input
  has a `fieldFree` block that `FieldFree_Fabricate` did not read, and for the
  seven former keys `method.mb.orbBased.exactFieldFree*`.
- MCTDHX, TDKS and TD-2RDM (without weak orbitals) bind
  `Method_PropagateExactPart` to a wrapper that calls `FieldFree_Propagate`
  on the orbitals of the state; MCTDHX and TDKS also bind
  `Method_ApplyExactPartPhi` if `FieldFree_ApplyPhi` is associated. For
  TDCI, TDHX and TD-2RDM-StaticOrbs a `fieldFree` block stops at fabricate.
- Spectral cut (`fieldFree.eigen.cut`): MCTDHX and TDKS project their linear
  and nonlinear orbital derivatives with `FieldFree_ProjectBelowCut`.
- Frozen mean-field reference (`fieldFree.cayley.meanFieldQ`): MCTDHX and
  TDKS build U at their first exact-part or nonlinear call
  (`FieldFree_BuildMeanField` with the orbitals and the 1-RDM of the state)
  and subtract U in the nonlinear part (`FieldFree_ApplyMeanField`; MCTDHX
  also −Û c in the CI derivative).

`Method_PropagateExactPart` and `Method_ApplyExactPartPhi` are owned by
`Method_Fabricate`, which nullifies them before the method fabricates bind
them; programs and tests reach those fabricates only through
`Method_Fabricate`, so no binding of an earlier fabricate or of a test
survives.

### Adiabatic Switching (`adiabaticSwitching`)

**What.** `"adiabaticSwitching": {"switchTime": T}` in `method.mb.orbBased`
propagates for −T < t < 0

```
H(t) = h + (1 − λ) U0 + λ W,   λ = sin²(π (t + T) / (2T))
```

and from t = 0 on the physical H = h + W (λ = 0 before −T,
`Method_Mb_OrbBased_SwitchFactor`). U0 is the frozen Hartree–Fock mean field
J0 − K0 of the initial state (`Method_Mb_OrbBased_SetSwitchingReference`, at
the end of the method setup). For a Hartree–Fock determinant h + U0 is its
Fock operator: at t = −T the determinant is an eigenstate, the closures of
TD-2RDM are exact for it and the collision operator (linear in W) vanishes,
and a slow switch carries it into the correlated ground state of the method.
For TD-2RDM this is the ground state of the TD-2RDM equations themselves,
without any CI (`twoRdmInit.excited` in the Hartree–Fock orbitals plus
correlating orbitals): app/atoms/ne/groundstate/td2rdm (moving orbitals) and
app/molecules/c2h4/dynamics/td2rdm (static orbitals, h2 from a file).

**Where.** Each term enters in one place, so every method built from the
shared operators follows without code of its own (TDHX, MCTDHX, TD-2RDM,
TD-2RDM-StaticOrbs; not TDCI, which reads the stored h2 in Coeffs, and not
TDKS, which has no W):
- λ W scales the source of every pair potential (`ApplyInteraction`,
  `FillPairPotentialsYlmNodes`): h2, the correlation and Hartree–Fock
  operators. The setups at t = 0 store the physical h2; TD-2RDM-StaticOrbs
  multiplies its h2 terms (commutator and collision operator, both linear in
  h2) by λ.
- (1 − λ) U0 is part of the one-body operator: for moving orbitals the grid
  operator (`T_Method_Mb_OrbBased_MeanField`, J0 and the ACE vectors of K0,
  the same frozen mean field as the Cayley reference) in `ApplySingleBodyOp`
  and, with a `fieldFree` block, in the linear split-step part (`ApplyDrivingOp`
  = gauge coupling + reference); for static orbitals the matrix
  U0(i,j) = 2 Σ_kl (h2(i,k,j,l) − h2(i,k,l,j)) rdm1(l,k) of the stored h2,
  added by `Method_Mb_OrbBased_RefillStaticH1` (correct also for an h2 that a
  program replaced; it then calls `SetSwitchingReference` again).

Since everything at t ≥ 0 is physical, a ground-state run propagates from −T
to 0 and a dynamics run continues from there (laser pulses start at t = 0).
Programs that change the state after `Method_Setup` (orthonormalization of
loaded correlating orbitals) call `Method_Mb_OrbBased_SetSwitchingReference`
again.

**Regularization.** Use `regularizationParameter` 1e-14 (default 1e-10).
The correlating orbitals start empty, and n/(n² + ε) throttles their
equations while n < √ε; a slow switch stays longer in that range. He1d
MCTDHF (T_He1d_14) from the determinant of its natural orbitals: E(0) minus
the imaginary-time ground state is 2.8e-5, 8.8e-5 and 1.5e-3 at T = 10, 20
and 40 with ε = 1e-10, and 1.2e-5, 2.7e-6 and 2.9e-6 at T = 10, 20 and 80
with ε = 1e-16 (the rest: the start is not the Hartree–Fock determinant of
these orbitals). TD-2RDM (T_He1d_15, NY without ccQ) ends 3.8e-7 from that
state.

**Results.** Ne, 9 orbitals, moving orbitals (app/atoms/ne/groundstate/td2rdm,
NY with ccQ, ε = 1e-14): E(0) = −128.67801 at T = 20 and −128.67765 at
T = 40, 97% of the MCSCF correlation energy of these orbitals (−128.68183),
natural occupations of the correlating orbitals within 2% of MCSCF; with ε = 1e-10 it ends at
−128.67437 (T = 20) with 10% occupation oscillations afterwards. Ethylene, 10
active orbitals, static orbitals (P_KickTd2rdmFem, T = 40): −32.8278729
against FCI −32.8278135; the earlier switching of that program from the
diagonal of the Fock matrix gave −32.8278734.

### Time-Dependent Kohn–Sham (TDKS, real-time TDDFT)

`method.mb.orbBased.tdks` (`M_Method_Mb_OrbBased_Tdks`) propagates Kohn–Sham
orbitals (fermions; occupations from the static CI vector of
`coeffsInit.unary`, as TDHX) with a functional of `sysXc`
(`src/SysXc/A_SysXc.md`):

```
i ∂φ_i/∂t = (h + v_KS,s[ρ↑, ρ↓] − a K_s) φ_i,   v_KS,s = J + v_xc,s,
```

h = T + V_ext (+ the SysGauge laser), the adiabatic local Kohn–Sham potential
of the instantaneous spin densities and the fraction a =
`SysXc_exactExchangeFraction` of exact exchange of a hybrid within spin s
(K_s φ_i = Σ_l V[φ_l* φ_i] g_l, g = φ R_s, over unordered pairs). The
equations are unprojected: v_KS − a K is Hermitian, so the orbitals stay
orthonormal.

**Spin channels.** Restricted orbitals (`orbs.restrictedQ`) must be closed
shell (equal spin blocks of the 1-RDM, checked at setup): one density (R↑ +
R↓), the unpolarized functional, the exchange of R↑ for both spins.
Unrestricted orbitals: the polarized functional with ρ↑, ρ↓ of body types 1
and 2 (one body type: ρ↓ = 0); `Method_Mb_OrbBased_Tdks_nSpin`.

**Energy.** `Method_GetEnergy` is the Kohn–Sham energy
E = Σ_i ⟨g_i|hφ_i⟩ + E_H + E_xc − a E_x (E_x = ½ Σ_s Σ_i ⟨g_i|K_s φ_i⟩), its
parts in `Method_Mb_OrbBased_Tdks_oneBodyEnergy`, `_hartreeEnergy`,
`_xcEnergy`, `_exactExchangeEnergy`. The orbital equations are its
Hamiltonian flow (v_xc is the exact derivative of the discrete E_xc), so it
is conserved without a laser; the Kohn–Sham SCF (`groundSolver.scf.stdImpl`
with `sysXc`) minimizes it.

**Propagation.** `Method_TimeDerivative` (all terms; imaginary time with
`Method_ImagTimeDerivative` and re-orthonormalization by the program), and
the split-step pair of MCTDHX: `Method_Mb_OrbBased_TimeDerivativeOrbsLin`
(−i h φ, or only the gauge coupling with a `fieldFree` block) and
`Method_Mb_OrbBased_TimeDerivativeCoeffsPlusOrbsNonLin` (−i (v_KS − a K) φ).
With a `fieldFree` block TDKS binds `Method_PropagateExactPart`
(`fieldFree.eigen` on Ylm grids, also `Method_ApplyExactPartPhi` for the
exponential integrators; `fieldFree.cayley` on grid.fem), and with
`fieldFree.cayley.meanFieldQ` the frozen Kohn–Sham reference
(src/FieldFree/A_FieldFree.md): strongly recommended, it removes the
splitting error of the ground state.

**Time step with a GGA.** A GGA potential depends on ∇ρ, so its kernel
δv_xc/δρ acts like a second-order differential operator (∼ 4|v_σ|ρ k²; up to
0.05 k² for H2 PBE), while the LDA, Hartree and exchange kernels stay bounded
at large k. In the split step the exact part rotates a mode of h0 energy E by
E Δt and the nonlinear part shears it in proportion to Δt times the kernel, so
the composition resonates near E Δt = nπ. Because the GGA kernel grows with E,
those instability bands widen, and a GGA run is stable only if every
propagated h0 mode has E Δt < π. The FEDVR and FEM spectra reach 10³–10⁵ a.u.
at the dense nodes (nuclei, element boundaries). This is a property of the
splitting: an adaptive integrator for the nonlinear part (`gslOdeiv2` rkf45)
diverges the same way, and the frozen reference does not help, since it
removes the potential, not the kernel. LDA runs at the TDHF steps: H2 at
Δt 0.02 (T_FemDft_03, energy conserved to 2e-11) and H2O at Δt 0.1, whose
O core has hMin 0.1 (E − E(0+) 7e-8 after the kick, drifting by 2e-9 per
0.5 a.u. over 3 a.u., no growth).

- Ylm: set `fieldFree.eigen.cut` = Ec with Ec Δt < π. He PBE at Δt 0.05:
  Ec 30 and 60 stable (|E − E0| < 2e-6 over 8 a.u.), Ec 120 not (2e-2),
  no cut ~0.1. At Ec 40 the cut raises the PBE ground-state energy by
  1.6e-5 Ha (T_He3d_04).
- grid.fem: the Cayley step has no eigenvectors (no cut), and its rotation
  angles 2 atan(E Δt / 2) crowd towards π for the stiff nodes, so reduce
  Δt instead. H2 PBE (hMin 0.3, nLocals 5): Δt 0.005 is stable (|E − E0|
  1.5e-9 over 3 a.u., from the GMRES tolerance), Δt 0.02 and 0.05 grow
  exponentially (1e-5 at t = 3 for 0.05).

**Refusals.** `sysXc` with any other method (a density functional would count
correlation twice with a correlated wavefunction; `method.sb` and the grid-
and geminal-based branches would ignore it), TDKS without `sysXc`, bosons,
more than two body types, open-shell restricted orbitals.

**Tests.** `test/simulationTests/FemDft/T_FemDft_03_H2KickTdks` (kick: N κ²/2,
TRK slope, energy and particle number),
`test/simulationTests/He3d/T_He3d_04_TimePropagationTdksPbeCut` (PBE at
Δt 0.05 with the spectral cut: energy conservation, no weight above the
cut), and the ground states of
`test/simulationTests/He3d`, `Ne3d/T_Ne3d_18` and `FemDft` (SCF with the TDKS
energy).

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
  files). This also covers orbitals embedded from a grid with smaller lmax
  (`orbsInit.load.lmaxIn`).
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
  `FillH2YlmNodes` (also used by `Method_GetEnergy`; both in the nested
  submodule `S_Method_Mb_OrbBased_YlmNodes`); the exact field-free
  part transforms only the declared columns (`fieldFree.eigen`,
  src/FieldFree/A_FieldFree.md).
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
  "grid": { "linear": { "xmin": -10.0, "xmax": 10.0, "const": { "nPoints": 128 } } },
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
   - Declare the options that the method implements.
     `Method_Mb_OrbBased_Fabricate` refuses a requested option that the
     selected method did not declare, so a new method needs no entry in the
     parent: the exact field-free part (a `fieldFree` block) by binding
     `Method_PropagateExactPart`, `fieldFree.cayley.meanFieldQ`,
     `fieldFree.eigen.cut` and `adiabaticSwitching` by setting
     `Method_Mb_OrbBased_implementsMeanFieldReferenceQ`,
     `_implementsSpectralCutQ` and `_implementsSwitchingQ`; a Kohn–Sham
     method sets `Method_kohnShamQ`.

4. **Register in parent fabricate:**
   Add a branch in `S_Method_Mb_OrbBased.f90`:
   ```fortran
   else if (Json_GetExistence("method.mb.orbBased.myMethod")) then
     call Method_Mb_OrbBased_MyMethod_Fabricate
   ```

5. **Reconfigure** (`cmake -B build`): CMake collects `src/*.f90` by itself,
   no list of source files has to be edited.

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

TD-2RDM passes its current orbitals to `Reconstruction_SetOrbitals` before
every `Reconstruction_AddCollisionOp` (full derivative and non-linear split
part); a reconstruction that does not depend on the orbitals ignores them.
TD-2RDM-StaticOrbs does not call it.

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
| `Mb/OrbBased/S_Method_Mb_OrbBased_YlmNodes.f90` | Definite-m node path of the correlation operator and h2 (Ylm, `mOfOrbs`); nested submodule of `S_Method_Mb_OrbBased` |
| `Mb/OrbBased/S_Method_Mb_OrbBased_MeanField.f90` | Frozen Hartree–Fock mean field of a 1-RDM as a value (`Method_Mb_OrbBased_FillMeanField`, `_AddMeanField`), used by the Cayley reference of FieldFree and the adiabatic switching |
| `Mb/OrbBased/Tdci/*`                  | Full CI with static orbitals              |
| `Mb/OrbBased/Tdhx/*`                  | Time-dependent Hartree–Fock/Exchange      |
| `Mb/OrbBased/Mctdhx/*`                | Multiconfiguration TDHX                   |
| `Mb/OrbBased/Td2rdm/*`                | TD-2RDM with dynamic orbitals             |
| `Mb/OrbBased/Td2rdmStaticOrbs/*`      | TD-2RDM with static orbitals              |
| `Mb/OrbBased/Tdks/*`                  | Time-dependent Kohn–Sham (TDDFT, sysXc)   |
| `Mb/GridBased/M_Method_Mb_GridBased.f90` | Grid-based interfaces                  |
| `Mb/GridBased/S_Method_Mb_GridBased.f90` | Grid-based operator implementations    |
| `Mb/GridBased/Full/*`                 | Full tensor-product grid method           |
| `Mb/GemBased/*`                       | Gaussian expansion (placeholder)          |

---

## Contact & Further Reading

- Main documentation: `README.md` in repository root
- Build instructions: `CMakeLists.txt`, `setCodyVars.sh`
- Related modules: `M_Propagator`, `M_Orbs`, `M_Coeffs`, `M_TwoRdm`
