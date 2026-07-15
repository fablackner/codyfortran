# IntegratorList Module — Agent Onboarding Guide

## Overview

The **IntegratorList** module provides a pluggable ODE integrator framework for
time-propagating quantum states in CodyFortranRDM. It follows the project's
_Interface Module + Implementation Submodule_ paradigm: a single abstract API
(`M_IntegratorList`) defines the contract, while concrete backends live in
separate submodules—each implementing a different numerical method.

At runtime the user's JSON configuration drives **fabrication**: the appropriate
integrator(s) are instantiated and stored in a global registry, ready to be
invoked by the `Propagator` layer.

---

## Architecture

```
M_IntegratorList (interface)
│
├── T_IntegratorList_E           (abstract base type)
│   ├─ Fabricate()               read config, bind pointers
│   ├─ Setup()                   allocate work arrays
│   └─ Integrate(state,t0,t1)    advance state in-place
│
├── T_IntegratorList_Container   (polymorphic wrapper)
│
└── integratorList(:)            module-global registry
```

### Concrete Backends

| Directory       | Key Type                          | Method                                  |
|-----------------|-----------------------------------|-----------------------------------------|
| `Rk/O1Expl`     | `T_IntegratorList_E_Rk_O1Expl`    | Forward Euler (1st-order explicit)      |
| `Rk/O2Expl`     | `T_IntegratorList_E_Rk_O2Expl`    | Midpoint (2nd-order explicit)           |
| `Rk/O2Impl`     | `T_IntegratorList_E_Rk_O2Impl`    | Implicit midpoint (2nd-order implicit)  |
| `Rk/O4Expl`     | `T_IntegratorList_E_Rk_O4Expl`    | Classical RK4 (4th-order explicit)      |
| `CrankNicolson` | `T_IntegratorList_E_Cn`           | Crank–Nicolson (A-stable, symplectic)   |
| `Expokit`       | `T_IntegratorList_E_Expokit`      | Krylov subspace exp(Δt·A)v              |
| `GslOdeiv2`     | `T_IntegratorList_E_GslOdeiv2`    | GSL adaptive steppers (rkf45, rk8pd, …) |
| `Sil`           | `T_IntegratorList_E_Sil`          | Short Iterative Lanczos                 |

All backends share the **three-stage lifecycle**:

1. **Fabricate** — parse JSON, store config path, bind `TimeDerivative` pointer.
2. **Setup** — allocate internal buffers once dimensions are known.
3. **Integrate** — perform the actual time step.

---

## Data Flow

```
┌────────────────────────────────────────────────────────────────────────────┐
│                           JSON Configuration                               │
│  "integratorList": [ { "rk": { "o4Expl": {} } }, { "expokit": {...} } ]   │
└────────────────────────────────────────────────────────────────────────────┘
                                     │
                                     ▼
             ┌───────────────────────────────────────────────┐
             │        IntegratorList_Fabricate(input)        │
             │  • Reads child names from JSON                │
             │  • Dispatches to backend *_Allocate routines  │
             │  • Attaches TimeDerivative callback           │
             │  • Calls element%Fabricate()                  │
             └───────────────────────────────────────────────┘
                                     │
                                     ▼
             ┌───────────────────────────────────────────────┐
             │           IntegratorList_Setup()              │
             │  • Iterates over integratorList(:)            │
             │  • Calls element%Setup()                      │
             └───────────────────────────────────────────────┘
                                     │
                                     ▼
             ┌───────────────────────────────────────────────┐
             │    Propagator calls element%Integrate(...)    │
             │  • state is updated in-place                  │
             │  • TimeDerivative callback computes d/dt      │
             └───────────────────────────────────────────────┘
```

---

## Key Interfaces

### `I_TimeDerivative`

```fortran
subroutine I_TimeDerivative(dState, state, time)
  complex(R64), intent(out), contiguous, target :: dState(:)
  complex(R64), intent(in),  contiguous, target :: state(:)
  real(R64),    intent(in)                      :: time
end subroutine
```

The user-supplied callback that computes **d(state)/dt**. In quantum dynamics
this is typically **−i Ĥ |ψ⟩**, but the interface is general for any first-order
ODE system.

### `I_Integrate`

```fortran
subroutine I_Integrate(this, state, t0, t1)
  class(T_IntegratorList_E), intent(inout) :: this
  complex(R64), intent(inout), contiguous  :: state(:)
  real(R64),    intent(in)                 :: t0, t1
end subroutine
```

Advances `state` from `t0` to `t1` **in-place**.

---

## JSON Configuration Examples

### Runge–Kutta (classical RK4)

```json
{
  "integratorList": [
    { "rk": { "o4Expl": {} } }
  ]
}
```

### Expokit Krylov exponential

```json
{
  "integratorList": [
    {
      "expokit": {
        "krylov_dim": 30,
        "tolerance": 1.0e-7,
        "max_steps": 1000
      }
    }
  ]
}
```

### GSL adaptive stepper

```json
{
  "integratorList": [
    { "gslOdeiv2": { "stepperType": "rkf45" } }
  ]
}
```

Supported GSL stepper types: `rk4`, `rkf45`, `rkck`, `rk8pd`.

### Short Iterative Lanczos (SIL)

```json
{
  "integratorList": [
    {
      "sil": {
        "krylovDim": 32,
        "tolerance": 1.0e-10,
        "maxRestarts": 8
      }
    }
  ]
}
```

The Krylov dimension adapts per step: the Lanczos basis stops growing as
soon as the residual error estimate drops below `tolerance` (up to
`krylovDim` vectors). If the tolerance cannot be met, the step is halved
and retried up to `maxRestarts` times. Requires a Hermitian generator
(TimeDerivative = −i Ĥ ψ with Ĥ Hermitian).

---

## Adding a New Integrator

1. **Create directory** `src/IntegratorList/MyMethod/`.
2. **Interface module** `M_IntegratorList_MyMethod.f90`:
   - Define `T_IntegratorList_E_MyMethod` extending `T_IntegratorList_E`.
   - Declare `IntegratorList_MyMethod_Allocate`.
3. **Submodule** `S_IntegratorList_MyMethod.f90`:
   - Implement `Allocate`, `Fabricate`, `Setup`, `Integrate`.
4. **Register in** `S_IntegratorList.f90`:
   - Add `use M_IntegratorList_MyMethod`.
   - Extend the `if/elseif` chain for child-name dispatch.
5. **Update CMakeLists** to include new source files.

---

## Numerical Method Summary

| Method           | Order | Explicit | Stability               | Best For                        |
|------------------|:-----:|:--------:|-------------------------|---------------------------------|
| O1Expl (Euler)   |   1   |    ✓     | Conditionally stable    | Baseline / debugging            |
| O2Expl (midpt)   |   2   |    ✓     | Conditionally stable    | Smooth, non-stiff dynamics      |
| O2Impl (midpt)   |   2   |          | A-stable                | Mildly stiff / oscillatory      |
| O4Expl (RK4)     |   4   |    ✓     | Conditionally stable    | General purpose, smooth RHS     |
| Crank–Nicolson   |   2   |          | A-stable, symplectic    | Unitary propagation (H indep.)  |
| Expokit (Krylov) |  var  |    ✓*    | Exact for linear        | Large sparse Hamiltonians       |
| GSL odeiv2       |  var  |   both   | Adaptive error control  | Production runs, stiff/non-stiff|
| SIL (Lanczos)    |  var  |    ✓*    | Exact for Hermitian     | Memory-limited, Hermitian H     |

*Krylov/Lanczos methods are "explicit" in the sense that they do not require
solving linear systems, but they approximate the matrix exponential directly.

---

## Tips for Contributors

- **Preserve style**: 2-space indentation, `!>` Doxygen comments, `M_`/`S_`
  naming, `procedure(I_...)` abstract interfaces.
- **Lifecycle discipline**: Never skip `Fabricate` → `Setup` → `Integrate`.
- **Avoid global state leaks**: All mutable data lives in the type instance or
  the `integratorList` registry.
- **Testing**: Add a JSON + expected-output pair under `test/` when introducing
  a new backend.

---

## See Also

- `M_Propagator` / `S_Propagator` — higher-level time-loop orchestration.
- `M_Method` — supplies the `TimeDerivative` callback for various physics.
- `M_Utils_ExpokitLib`, `M_Utils_Odeiv2GslLib` — low-level library wrappers.
