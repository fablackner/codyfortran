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
│   ├── S_Propagator_SplitStep.f90    # Dispatch to order variants
│   ├── Order2/
│   │   ├── M_Propagator_SplitStep_Order2.f90  # 2nd-order Strang splitting
│   │   └── S_Propagator_SplitStep_Order2.f90
│   └── Order4/
│       ├── M_Propagator_SplitStep_Order4.f90  # 4th-order Yoshida composition
│       └── S_Propagator_SplitStep_Order4.f90
└── EigenExpansion/
    ├── M_Propagator_EigenExpansion.f90   # Diagonalization-based propagation
    └── S_Propagator_EigenExpansion.f90   # Exact exp(−iEᵢΔt) phase evolution
```

### Design Pattern: Strategy + Facade

| Component | Role |
|-----------|------|
| **M_Propagator** | Facade exposing `Propagator_Setup` and `Propagator_Propagate` pointers |
| **S_Propagator** | Factory reading JSON and wiring procedure pointers to backends |
| **Single** | Monolithic integration (no splitting) via IntegratorList |
| **SplitStep** | Suzuki–Trotter operator decomposition |
| **EigenExpansion** | Exact diagonalization + phase factors |

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
│  2. Dispatch to Single/SplitStep/EigenExpansion based on key               │
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

### When to Use Each

- **Single:** Default choice. Let the integrator (RK4, SIL, Expokit) handle everything.
- **SplitStep:** When Ĥ = T̂ + V̂ where each part is cheap to apply (e.g., FFT for T̂, diagonal for V̂).
- **EigenExpansion:** Small systems (dim ≲ 10⁴), time-independent Ĥ, need many time steps.

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

## Typical Usage in CodyFortranRDM

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

---

## Testing

Relevant test patterns:
- `T_*_Tdci*` — Tests using TDCI with various propagators
- `T_*_Mctdhx*` — Tests using MCTDHX dynamics

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
