# Propagator Module — Agent Onboarding Guide

## Purpose

The **Propagator** module provides a unified interface for **time evolution** of quantum states. It abstracts away the differences between various propagation algorithms, allowing the rest of the codebase to advance states through a single facade.

**Core problem solved:** Given a quantum state |Ψ(t₀)⟩ and a Hamiltonian Ĥ, compute:

```
|Ψ(t₁)⟩ = Û(t₁, t₀) |Ψ(t₀)⟩ = exp(−iĤ(t₁−t₀)) |Ψ(t₀)⟩
```

The propagator determines *how* the time-evolution operator Û is approximated or computed exactly.

---

## Architecture Overview

```
src/Propagator/
├── M_Propagator.f90              # Public facade (procedure pointers + abstract interfaces)
├── S_Propagator.f90              # Factory logic (JSON dispatch to backends)
├── AGENTS.md                     # This file
├── Single/
│   ├── M_Propagator_Single.f90       # Single-integrator backend interface
│   └── S_Propagator_Single.f90       # Delegates to IntegratorList(1)
├── SplitStep/
│   ├── M_Propagator_SplitStep.f90    # Split-operator framework interface
│   ├── S_Propagator_SplitStep.f90    # Order dispatch, exact-part composition (IntegrateLin)
│   ├── Order2/
│   │   ├── M_Propagator_SplitStep_Order2.f90  # 2nd-order Strang splitting
│   │   └── S_Propagator_SplitStep_Order2.f90
│   └── Order4/
│       ├── M_Propagator_SplitStep_Order4.f90  # 4th-order Yoshida composition
│       └── S_Propagator_SplitStep_Order4.f90
├── EigenExpansion/
│   ├── M_Propagator_EigenExpansion.f90   # Diagonalization-based propagation
│   └── S_Propagator_EigenExpansion.f90   # Exact exp(−iEᵢΔt) phase evolution
└── Exponential/
    ├── M_Propagator_Exponential.f90      # Exponential integrators: factory + shared building blocks
    ├── S_Propagator_Exponential.f90      # Scheme dispatch, Krogstad step, history helpers
    ├── Adams/
    │   ├── M_Propagator_Exponential_Adams.f90        # Exponential Adams PECE
    │   └── S_Propagator_Exponential_Adams.f90
    └── AdamsLawson/
        ├── M_Propagator_Exponential_AdamsLawson.f90  # Adams PECE in the interaction picture
        └── S_Propagator_Exponential_AdamsLawson.f90
```

### Design Pattern: Strategy + Facade

| Component | Role |
|-----------|------|
| **M_Propagator** | Facade exposing `Propagator_Setup` and `Propagator_Propagate` pointers |
| **S_Propagator** | Factory reading JSON and wiring procedure pointers to backends |
| **Single** | Monolithic integration (no splitting) via IntegratorList |
| **SplitStep** | Suzuki–Trotter operator decomposition |
| **EigenExpansion** | Exact diagonalization + phase factors |
| **Exponential** | Exponential integrators around the method's exactly propagated part |

---

## Data Flow

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                           JSON Configuration                                 │
│  "propagator": { "single": {} }  OR  "splitStep": { "order2": {} }  etc.   │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                         Propagator_Fabricate()                               │
│  1. Read JSON "propagator" block                                            │
│  2. Dispatch to Single/SplitStep/EigenExpansion/Exponential based on key   │
│  3. Backend binds Propagator_Setup, Propagator_Propagate pointers          │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                         Propagator_Setup()                                   │
│  Backend-specific initialization (e.g., diagonalize Ĥ for EigenExpansion)  │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                    Propagator_Propagate(state, t0, t1)                       │
│  Evolve state in-place from t₀ to t₁ using the selected algorithm          │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## JSON Configuration

### Single Integrator

Direct delegation to `IntegratorList(1)`. Use when the integrator handles the full Hamiltonian.

```json
{
  "propagator": {
    "single": { }
  },
  "integratorList": {
    "rk4": { "nTimeSteps": 100 }
  }
}
```

### Split-Step (2nd Order)

Symmetric Strang splitting: A(½Δt) · B(Δt) · A(½Δt)

```json
{
  "propagator": {
    "splitStep": {
      "order2": { }
    }
  },
  "integratorList": {
    "kineticIntegrator": { ... },
    "potentialIntegrator": { ... }
  }
}
```

**Requires 2 integrators:** IntegratorList(1) for Â, IntegratorList(2) for B̂.

### Split-Step (4th Order)

Yoshida composition with 7 stages and O(Δt⁴) global error:

```json
{
  "propagator": {
    "splitStep": {
      "order4": { }
    }
  }
}
```

### Split-Step: exactly propagated linear part

Every Â step goes through `Propagator_SplitStep_IntegrateLin` (both orders):
if the method propagates part of Â exactly (`Method_PropagateExactPart`
associated, e.g. `method.mb.orbBased.exactFieldFreeQ`), the step is a
composition of exact and IntegratorList(1) substeps selected by
`exactPartScheme`; with the default `strang` scheme it is
exact(Δt/2) · IntegratorList(1)(Δt) · exact(Δt/2). Otherwise it is
IntegratorList(1) alone.

`"exactPartScheme"` in `propagator.splitStep` selects the composition of the
exact part and IntegratorList(1) inside that step:

| scheme | order | exact / integrator substeps |
|---|---|---|
| `"strang"` (default) | 2 | 2 / 1 |
| `"yoshida4"` | 4 | 4 / 3 (negative middle steps) |
| `"blanesMoan6"` | 4 | 7 / 6 (Blanes & Moan 2002, PRK S6; ~250× smaller error constant than Yoshida) |
| `"blanesMoan6Swap"` | 4 | 6 / 7 (same composition with the roles exchanged; starts and ends with an integrator substep) |

`blanesMoan6Swap` saves one exact substep per step, worthwhile when the exact
part is the expensive one (dense transformations of the exact field-free
propagator); for the Be HHG runs its error equals that of `blanesMoan6`
(git history of app/be/README.md).

The exact substeps carry the time; every integrator substep is centred on the
time reached, so its midpoint evaluation of a time-dependent part (e.g. the
vector potential) stays consistent with the symmetric composition.

**Exact part on grid.fem (Cayley step).** With `exactFieldFreeQ` on grids
without a dense h0 eigen-decomposition (grid.fem) the exact part is a
GMRES-solved Crank–Nicolson (Cayley) step of h0 = T + V_static. It is unitary
and time-reversible, Cayley(−Δt) = Cayley(Δt)⁻¹ (negative substeps are fine
in real time), but only second order, so use `exactPartScheme: "strang"`: the
fourth-order schemes spend 2–3.5× more exact substeps without raising the
order. Real time only (an integrator 1 with `imagTimeQ` stops in the exact
part); a laser is the length gauge `sysGauge.fem.length` (E(t) e_pol·r),
which integrator 1 (e.g. RK4) propagates. For molecules with heavier atoms
(a core) also set `method.mb.orbBased.exactFieldFreeMeanFieldQ` (MCTDHX): it
moves a frozen Fock mean field U = J0 − K0 into the exact part and subtracts
it in the nonlinear part, a pure change of the splitting, without which the
split h0 | mean field needs tiny time steps there; the composition and the
integrators stay as they are. Details: A_Method.md, "Exact Field-Free Part
on grid.fem: Cayley Step" and "Frozen Mean-Field Reference of the Cayley
Step" (which state to build U from, restarts); measured drifts and costs:
git history of app/fem_molecules/README.md.

```json
{
  "propagator": {
    "splitStep": {
      "order2": { },
      "exactPartScheme": "blanesMoan6"
    }
  }
}
```

### Eigen-Expansion

Diagonalize once, apply phases forever. Best for time-independent Ĥ:

```json
{
  "propagator": {
    "eigenExpansion": { }
  },
  "diagonalizerList": {
    "lapack": { "nEvals": -1 }
  }
}
```

### Exponential Integrators

For u' = L u + N(t, u) with L the exactly propagated part of the method
(`Method_PropagateExactPart` and its phi functions `Method_ApplyExactPartPhi`,
e.g. h0 with `method.mb.orbBased.exactFieldFreeQ`) and N the sum of the
TimeDerivative callbacks of IntegratorList(1) and IntegratorList(2) (their
Integrate procedures are not used). There is no operator splitting between L
and N; real time only.

| scheme | order | N per step | notes |
|---|---|---|---|
| `"krogstad4"` | 4 | 4 | exponential Runge–Kutta (Krogstad 2005), one step |
| `"adams"` | k + 1 | 2 | exponential Adams PECE, `nSteps` = k (1..6, default 4) |
| `"adamsLawson"` | k + 1 | 2 | Adams PECE in the interaction picture of L (Auzinger et al. 2021), `nSteps` = k |

Not available for the Cayley exact part of grid.fem (`exactFieldFreeQ` on
non-Ylm grids): it has no phi functions, MCTDHX binds no
`Method_ApplyExactPartPhi`, and `Propagator_Exponential_Fabricate` stops with a
message pointing to `propagator.splitStep` (T_Fem_08_ExponentialUnsupported).

The multistep schemes need a fixed step size: a step that does not continue
the previous one (other t0 or step size) restarts the history with Krogstad
steps. A state changed between steps (absorber) is used as it is.

They interpolate N (or, for `adamsLawson`, N in the interaction picture) in
time, so components of N that oscillate with the phases of L limit the step:
h times the largest exactly propagated frequency (for the spectral cut,
h · `exactFieldFreeCut`) should stay well below 1. On the Be HHG benchmark
(velocity gauge, cut 1000, dt 0.0025) Strang splitting was the most accurate;
`adams` was stable only with a cut of 300 (1.6× faster, harmonics 21–60 about
0.2 decades further off), `adamsLawson` stayed stable at cut 1000 (1.75×
faster, 0.4 decades), and `krogstad4` cost as much as Strang at 0.25 decades
(git history of app/be/README.md). They are options for observables less sensitive than
high harmonics, e.g. ionization yields.

```json
{
  "propagator": {
    "exponential": {
      "adams": { "nSteps": 4 }
    }
  }
}
```

---

## Public Interface

```fortran
! 1. Fabricate (select backend from JSON)
call Propagator_Fabricate()

! 2. Setup (backend-specific initialization)
call Propagator_Setup()

! 3. Propagate (time evolution)
call Propagator_Propagate(state, t0, t1)
```

### Propagate Signature

```fortran
subroutine Propagator_Propagate(state, t0, t1)
  complex(R64), intent(inout), contiguous :: state(:)  ! Modified in-place
  real(R64), intent(in) :: t0   ! Start time
  real(R64), intent(in) :: t1   ! End time (can be < t0 for backward)
end subroutine
```

---

## Backend Comparison

| Backend | Method | Error (global) | Cost/step | Best For |
|---------|--------|----------------|-----------|----------|
| **Single** | Delegate to integrator | Integrator-dependent | 1 integrator call | General purpose, adaptive |
| **SplitStep/Order2** | Strang ABA | O(Δt²) | 3 integrator calls | Separable Ĥ = T̂ + V̂ |
| **SplitStep/Order4** | Yoshida ABABABA | O(Δt⁴) | 7 integrator calls | High accuracy, separable Ĥ |
| **EigenExpansion** | Exact diagonalization | Machine precision | O(n²) matmul | Small dim, long times |
| **Exponential** | Exponential RK / Adams around the exact part | O(Δt⁴) or O(Δt^(k+1)) | 2–4 N + phi combinations | Exact part dominant, no splitting |

### When to Use Each

- **Single:** Default choice. Let the integrator (RK4, SIL, Expokit) handle everything.
- **SplitStep:** When Ĥ = T̂ + V̂ where each part is cheap to apply (e.g., FFT for T̂, diagonal for V̂).
- **EigenExpansion:** Small systems (dim ≲ 10⁴), time-independent Ĥ, need many time steps.
- **Exponential:** The method propagates a stiff part exactly and provides its phi functions (MCTDHX with `exactFieldFreeQ`); no splitting error, but the step must resolve the exactly propagated frequencies coupled by N.

---

## Split-Step Mathematics

### Order-2 (Strang Splitting)

```
exp(−i(Â+B̂)Δt) ≈ exp(−iÂΔt/2) · exp(−iB̂Δt) · exp(−iÂΔt/2) + O(Δt³)
```

- **Local error:** O(Δt³)
- **Global error:** O(Δt²)
- **Time-reversible:** Yes (symmetric)

### Order-4 (Yoshida)

Composes three order-2 steps with carefully chosen weights:

```
γ₁ = 1/(2 − 2^(1/3)) ≈ 1.3512
γ₂ = 1 − 2γ₁ ≈ −1.7024

S₄(Δt) = S₂(γ₁Δt) · S₂(γ₂Δt) · S₂(γ₁Δt)
```

- **Local error:** O(Δt⁵)
- **Global error:** O(Δt⁴)
- **Note:** γ₂ < 0 means a "backward" substep

---

## Eigen-Expansion Details

### Algorithm

1. **Setup:** Diagonalize Ĥ → eigenvalues Eᵢ, eigenvectors |φᵢ⟩
2. **Propagate:**
   - Project: cᵢ = ⟨φᵢ|Ψ(t₀)⟩
   - Evolve: |Ψ(t₁)⟩ = Σᵢ cᵢ exp(−iEᵢΔt) |φᵢ⟩

### Complexity

- **Setup:** O(dim³) for dense diagonalization (once)
- **Propagate:** O(nFound × dim) per step

### Limitations

- Stores O(dim × nFound) eigenvector matrix
- If `nFound < dim`, components orthogonal to eigenbasis are lost
- Not suitable for time-dependent Ĥ

---

## Dependencies

| Depends On | Purpose |
|------------|---------|
| `M_Utils_Types` | R64, I32 type kinds |
| `M_Utils_Json` | JSON configuration parsing |
| `M_Utils_Say` | Logging/diagnostics |
| `M_Utils_Constants` | Imaginary unit IU |
| `M_Utils_NoOpProcedures` | Default no-op for Setup |
| `M_IntegratorList` | Time integrators (RK, SIL, etc.) |
| `M_DiagonalizerList` | Eigensolvers (for EigenExpansion) |

---

## Typical Usage in CodyFortran

```fortran
program time_evolution
  use M_Propagator
  use M_Method
  use M_IntegratorList

  ! ... Grid, Hamiltonian, Method setup ...

  ! 1. Setup integrator(s) with time derivative
  type(T_IntegratorList_FabricateInput) :: intInput(1)
  intInput(1)%TimeDerivative => Method_TimeDerivative
  call IntegratorList_Fabricate(intInput)
  call IntegratorList_Setup

  ! 2. Fabricate and setup propagator
  call Propagator_Fabricate()
  call Propagator_Setup()

  ! 3. Time stepping loop
  do iStep = 1, nSteps
    call Propagator_Propagate(Method_state, t, t + dt)
    t = t + dt
    ! ... observables, output ...
  end do
end program
```

---

## Adding a New Backend

1. **Create directory:** `src/Propagator/MyBackend/`

2. **Interface module** (`M_Propagator_MyBackend.f90`):
   ```fortran
   module M_Propagator_MyBackend
     use M_Utils_Types
     implicit none
     interface
       module subroutine Propagator_MyBackend_Fabricate
       end subroutine
     end interface
   end module
   ```

3. **Implementation submodule** (`S_Propagator_MyBackend.f90`):
   ```fortran
   submodule(M_Propagator_MyBackend) S_Propagator_MyBackend
     implicit none
   contains
     module subroutine Propagator_MyBackend_Fabricate
       use M_Propagator
       Propagator_Propagate => MyPropagate
       Propagator_Setup => MySetup  ! optional
     end subroutine

     subroutine MyPropagate(state, t0, t1)
       ! Your algorithm here
     end subroutine
   end submodule
   ```

4. **Register in factory** (`S_Propagator.f90`):
   ```fortran
   use M_Propagator_MyBackend
   ...
   else if (Json_GetExistence("propagator.myBackend")) then
     call Propagator_MyBackend_Fabricate
   ```

5. **Update CMakeLists.txt** to include new source files

---

## Common Pitfalls

1. **Missing IntegratorList setup:** Propagator backends (Single, SplitStep) require IntegratorList to be fabricated and set up first with a valid TimeDerivative callback.

2. **Wrong integrator count for SplitStep:** Order2 and Order4 expect exactly 2 integrators in IntegratorList — one for Â, one for B̂.

3. **EigenExpansion with time-dependent Ĥ:** The eigenbasis is computed once at t=0. Time-dependent Hamiltonians will give incorrect results.

4. **Forgetting Propagator_Setup:** EigenExpansion requires setup (to diagonalize). Other backends may have no-op setup but calling it is harmless.

5. **Incomplete eigenbasis:** If DiagonalizerList returns `nFound < dim`, EigenExpansion will project out missing components.

6. **Exact part on grid.fem:** `exactFieldFreeQ` there is a second-order Cayley step: use `splitStep` with `exactPartScheme: "strang"`; `exponential` and imaginary time stop with an error. For molecules with heavier atoms (a core) add `exactFieldFreeMeanFieldQ` (MCTDHX), otherwise the split h0 | mean field needs tiny steps.

---

## Testing

Relevant test patterns:
- `T_*_Tdci*` — Tests using TDCI with various propagators
- `T_*_Mctdhx*` — Tests using MCTDHX dynamics

- `T_Propagator_Exponential` — convergence order and norm conservation of
  the exponential integrators on a small model with exact test bindings
- `T_Fem_05`–`T_Fem_10` — split-step TDHF with the Cayley exact part on
  grid.fem (stationary ground state, dipole kick, exponential refused, frozen
  mean-field reference: ground state stationary at Δt = 0.1, second-order
  convergence of the kick to the trajectory of the plain split; in MCTDHF
  with nOrbs > nBodies the same trajectory as the plain split)

Run propagator-related tests:
```bash
cd build && ctest -R Tdci
```

---

## Performance Tuning

### Single Backend
- Performance depends entirely on the chosen integrator
- Use adaptive integrators (GslOdeiv2) for stiff problems
- Use Expokit/SIL for Krylov-based propagation

### SplitStep
- **Order2 vs Order4:** Order4 costs ~2.3× more per step but allows ~3× larger Δt for same accuracy
- Optimize individual integrators (e.g., FFT-based for kinetic)

### EigenExpansion
- One-time O(dim³) cost amortized over many time steps
- Use `nEvals=-1` (all eigenvalues) for exact propagation
- Consider ARPACK with partial spectrum if only low-energy dynamics matter

---

## Related Modules

- **Method** — Defines `Method_TimeDerivative` callback used by integrators
- **IntegratorList** — ODE integrators (RK, SIL, Expokit, GSL)
- **DiagonalizerList** — Eigensolvers for EigenExpansion backend
- **Absorber** — Boundary absorption applied after propagation
