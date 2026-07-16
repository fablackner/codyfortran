# Mixing Module

## Purpose

Abstracts self-consistency mixers (convergence accelerators/dampers) behind a
small interface. Given the currently accepted iterate `x` and a new candidate
`xNew` produced by one fixed-point step (e.g., one SCF cycle), a mixer
computes the next accepted iterate. Mixers are agnostic to what the vector
represents — consumers decide whether to feed orbitals, a density/potential,
or any other flattened complex quantity.

## Architecture Overview

Standard CodyFortranRDM fabrication pattern:

```
src/Mixing/
├── M_Mixing.f90          # Interface: T_Mixing handle + procedure pointers
├── S_Mixing.f90          # Dispatch: JSON → backend selection
├── Linear/
│   ├── M_Mixing_Linear.f90
│   └── S_Mixing_Linear.f90
└── Diis/
    ├── M_Mixing_Diis.f90
    └── S_Mixing_Diis.f90
```

The mixing *algorithm* is selected once per run via JSON (`mixing.*`), but
consumers may hold several independent `T_Mixing` instances (one per mixed
quantity), each with its own dimension and history. `Mixing_Fabricate` is
called by the program alongside the other module fabrications (before
`GroundSolver_Fabricate`, which consumes the mixer).

## Key Interfaces

```fortran
call Mixing_Fabricate()                      ! bind backend from JSON
call Mixing_SetupMixer(mixer, dim)           ! per mixed quantity
call Mixing_Mix(mixer, x, xNew, alpha)       ! x ← next accepted iterate
call Mixing_ResetMixer(mixer)                ! drop history
```

Contract of `Mixing_Mix`: `x` holds the accepted iterate on entry and the
next accepted iterate on exit; `xNew` is the candidate produced by one
fixed-point step applied to `x`; the residual is `r = xNew − x`;
`alpha ∈ (0, 1]` is the damping.

## Backends

### linear (default)

`x ← (1−α)·x + α·xNew`. Stateless. Selected automatically when the JSON has
no `mixing` block, which keeps legacy configurations bit-compatible.

```json
"mixing": { "linear": { } }
```

### diis (Pulay/Anderson)

Stores (iterate, residual) pairs in a circular buffer and minimizes
`‖Σ c_i·r_i‖²` subject to `Σ c_i = 1` (bordered system, SVD pseudo-inverse,
real coefficients), then steps `x ← Σ c_i·(x_i + α·r_i)`. With one stored
pair this is exactly linear mixing.

Safeguards (all trigger a history restart + plain damped step):
- **Gating**: extrapolation starts only once the residual norm has dropped
  below `startThreshold` × (first residual norm) — the linear model is only
  trusted in the quasi-linear regime
- **Growth**: residual norm grows by more than 3× in one step
- **Coefficients**: `Σ|c_i|` exceeds 1e4 (near-dependent residuals)

```json
"mixing": { "diis": { "nHistory": 8, "startThreshold": 0.1 } }
```

| Parameter        | Type | Default | Description                              |
|------------------|------|---------|------------------------------------------|
| `nHistory`       | int  | 8       | Maximum stored (iterate, residual) pairs |
| `startThreshold` | real | 0.1     | Relative residual drop enabling extrapolation |

## Consumers (GroundSolver)

All `GroundSolver_Approach` backends route their orbital update through a
mixer instance. `groundSolver.scf.stdImpl` additionally offers `mixTarget`:

- `"orbitals"` (default): mix the gauge-aligned Fock eigenvectors into the
  old orbitals, then orthonormalize
- `"potential"`: mix the Hartree potential (equivalent to density mixing —
  the interaction potential is linear in the source density) and replace the
  orbitals outright

The MCSCF CI coefficients are *not* routed through Mixing: an eigenvector is
not a fixed-point iterate in a linear space, so residual extrapolation does
not apply; they keep phase-aligned linear mixing.

## Common Pitfalls

- **Never feed raw eigenvectors**: remove gauge freedom first
  (`Orbs_AlignOnReference`), otherwise the residual `xNew − x` is
  dominated by arbitrary phases/rotations and any mixing diverges.
- **DIIS prefers gauge-invariant, linear-space targets**: densities and
  potentials. Orbital sets live on a manifold (orthonormalization after the
  mix breaks the linear model); DIIS on orbitals is allowed but usually no
  better than linear.
- **DIIS with very small α is counterproductive**: the exchange operator is
  built from the (lagging) orbitals, so with heavy damping the residual is
  not a clean function of the mixed potential and extrapolation misfires.
  Use α ≳ 0.5 and let the extrapolation provide the stability.
- One `T_Mixing` instance per mixed quantity — histories must not be shared
  between quantities of different meaning or dimension.

## Testing

Exercised through the ground-solver simulation tests
(`test/simulationTests/Ne3d/T_Ne3d_02_GroundSolverScf`,
`T_Ne3d_03_GroundSolverScfYlm`); the default (no `mixing` block) reproduces
the legacy linear behavior exactly.
