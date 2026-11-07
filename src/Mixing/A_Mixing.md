# Mixing Module

## Purpose

Abstracts self-consistency mixers (convergence accelerators/dampers) behind a
small interface. Given the currently accepted iterate `x` and a new candidate
`xNew` produced by one fixed-point step (e.g., one SCF cycle), a mixer
computes the next accepted iterate. Mixers are agnostic to what the vector
represents — consumers decide whether to feed orbitals, a density/potential,
or any other flattened complex quantity.

## Architecture Overview

Standard CodyFortran fabrication pattern:

```
src/Mixing/
├── M_Mixing.f90          # Interface: procedure pointers Mixing_Mix, Mixing_Reset
├── S_Mixing.f90          # Dispatch: JSON → backend selection
├── Linear/
│   ├── M_Mixing_Linear.f90
│   └── S_Mixing_Linear.f90
└── Diis/
    ├── M_Mixing_Diis.f90
    └── S_Mixing_Diis.f90
```

There is one mixer per run. Its state (vector length, history) lives in the
backend submodule; `Mixing_Fabricate` starts it empty (also when called
again), a new vector length restarts it, and `Mixing_Reset` empties it when
the fixed-point map changes while the iteration continues (another
functional, another system). `Mixing_Fabricate` is called by the program
alongside the other module fabrications (before `GroundSolver_Fabricate`,
whose backends call `Mixing_Mix`).

## Key Interfaces

```fortran
call Mixing_Fabricate()          ! bind the backend from JSON, empty mixer
call Mixing_Mix(xMixed, xRaw)    ! xMixed ← next accepted iterate
call Mixing_Reset()              ! drop the history (new fixed-point map)
```

Contract of `Mixing_Mix`: `xMixed` holds the accepted iterate on entry and
the next accepted iterate on exit; `xRaw` is the candidate produced by one
fixed-point step applied to `xMixed`; the residual is `r = xRaw − xMixed`.
The damping α is a backend parameter (`mixing.<backend>.alpha`).

## Backends

### linear (default)

`x ← (1−α)·x + α·xNew`. Stateless (`Mixing_Reset` is a no-op). Selected
automatically when the JSON has no `mixing` block, which keeps legacy
configurations bit-compatible.

```json
"mixing": { "linear": { "alpha": 0.7 } }
```

### diis (Pulay/Anderson)

Stores per step the residual r_i and the damped candidate y_i = x_i + α·r_i
in a circular buffer, minimizes ‖Σ c_i·r_i‖² subject to Σ c_i = 1 (bordered
system, pseudo-inverse of its symmetric eigen-decomposition, real
coefficients), then steps `x ← Σ c_i·y_i`. With one stored pair this is
exactly linear mixing with the same α.

Safeguards (each keeps only the newest pair, i.e. takes a plain damped step):
- **Start gate**: plain steps while the residual norm is above
  `startThreshold` × the first residual norm (since fabrication or the last
  `Mixing_Reset`): the linear model is trusted only in the quasi-linear
  regime. `startThreshold: 0` switches the gate off
- **Growth**: the residual norm grows by more than 3× in one step
- **Coefficients**: the coefficient sum degenerates or Σ|c_i| exceeds 1e4
  (near-dependent residuals)

```json
"mixing": { "diis": { "nHistory": 8, "startThreshold": 0.0, "alpha": 1.0 } }
```

| Parameter        | Type | Default | Description                              |
|------------------|------|---------|------------------------------------------|
| `nHistory`       | int  | 8       | Maximum stored pairs                     |
| `startThreshold` | real | 0.1     | Relative residual drop below which the extrapolation acts (0: no gate) |
| `alpha`          | real | 1.0     | Damping 0 < α ≤ 1                         |
| `printLevel`     | int  | 0       | 1: one line per step (residual norm, stored pairs, restarts) |

Memory: 2·`nHistory` vectors of the mixed length. For orbitals that is
2·nHistory·nPoints·nOrbs complex numbers: 1.4 GB for the 5 orbitals of H2O
on the accurate 1.1 M-node grid.fem mesh at nHistory 8. A shorter history
converges more slowly (H2O Hartree–Fock of T_Fem_04: nHistory 3, 4, 6, 8 →
17, 15, 14, 12 SCF steps).

## What to mix

DIIS models the residual as a linear function of the mixed quantity, so the
mixed quantity must be a *complete iterate*: the next Fock operator, and
with it the next candidate, must be a function of it alone.

- **Orbitals** are complete for every SCF-type map (Hartree–Fock, Kohn–Sham
  with any functional, the Fock operator of MCSCF stdImpl/ylmOpt) and for
  the MCSCF gradient step: every mean field is rebuilt from the mixed
  orbitals.
- **The Kohn–Sham potential** v_KS = J + v_xc of a local functional (LDA,
  GGA) is complete.
- **The Hartree potential** of Hartree–Fock (`mixTarget: "potential"`), v_KS
  of a hybrid and the Hartree potential of the MCSCF Fock operator are
  *not*: the exchange (and the correlated part) is rebuilt from the orbitals
  and lags the mixed potential, so the history pairs do not belong to one
  map. Linear mixing tolerates this (the lag acts as extra damping); DIIS
  extrapolates over inconsistent pairs and is slower than linear mixing.

Measured SCF steps to the test's convergence criterion:

| System (test setup) | linear on the potential | DIIS on the potential | DIIS on the orbitals (no gate, α 1) |
|---|---|---|---|
| Ne Hartree–Fock, Ylm (T_Ne3d_02) | 12 (α 0.85) | 87, 53, 17, 39 (α 0.3, 0.5, 0.8, 1) | 13 (linear on the orbitals: 20 at α 0.85, 28 at α 1) |
| Ne Hartree–Fock, Ylm per l (T_Ne3d_03, scf.ylmOpt) | 14 (α 0.85) | 104, 41, 22, 66 | 12 |
| Ne MCSCF Fock (mcscf.stdImpl, T_Ne3d_04) | 14 (α 0.85) | 68, 43, 17, 39 | 12 |
| He exact exchange only, Ylm (T_He3d_03) | 14 (α 0.8) | not converged in 40 (α 0.3, 0.5) | 9 |
| H2O Hartree–Fock, grid.fem (T_Fem_04) | 22 (α 0.7) | 51 (erratic) | 12 |
| Ne LDA, Ylm (T_Ne3d_18) | oscillates (α 0.8) | 13 (α 0.3) | 11 |
| He PBE, Ylm (T_He3d_02) | 8 (α 0.8) | 8 (α 0.5) | 5 |
| H2O LDA, grid.fem (app) | — | 14 (α 0.3; 564 Kohn–Sham products, 82 s) | 12 (629 products, 86 s) |

## Recommended settings

| Consumer | `mixTarget` | `mixing` |
|---|---|---|
| `groundSolver.scf.stdImpl` / `scf.ylmOpt` (Hartree–Fock, any functional), `mcscf.stdImpl` / `mcscf.ylmOpt` | `"orbitals"` | `"diis": { "startThreshold": 0.0 }` |
| `groundSolver.mcscf.gradient` | (always the orbitals) | `"diis": { }` (the default gate 0.1) |

For a local functional, DIIS on `"potential"` with α ≈ 0.3 works as well
(H2O LDA: 82 s against 86 s). The orbitals are the input that stays
correct when the functional becomes a hybrid.

The start gate is what differs, and either choice fails one of the two map
classes:

- **SCF-type maps** (eigenvectors of the Fock operator of the current
  orbitals) need DIIS from the second step on, because a gate can block it
  altogether:
  - Ne LDA at α 1: the plain orbital iteration is a 2-cycle whose residual
    hardly decays (39 → 25.6, then flat), so with `startThreshold` 0.1 DIIS
    never starts (not converged in 60 steps);
  - H2O LDA: the plain step overshoots, so the residual exceeds the first
    one; even `startThreshold` 1 blocks DIIS (a 2-cycle of E_el −78.4 /
    −73.8 for 60 steps). Without the gate it converges in 12 steps.
- **The MCSCF gradient step** is a preconditioned, Newton-like step whose
  norm need not decrease while the energy gradient does (weakly occupied
  orbitals). Early extrapolation leaves the basin. Be 1s 2s 2p 3s 3p
  (T_Be3d_02):
  - linear, α 1: 27 steps;
  - DIIS, α 1, gate 0.1: 19;
  - DIIS, α 1, no gate: not converged in 100;
  - DIIS, α 0.5, no gate: 23.

  With the default gate, DIIS is never worse than linear mixing with the
  same α: before the gate opens it *is* linear mixing.

## Implementation notes (performance)

- **Incremental Pulay matrix.** B_ij = Re⟨r_i|r_j⟩ is kept across calls:
  one new pair adds one row, ⟨r_i|r_new⟩ for all stored i, from a single
  ZGEMV over the stored residuals. Rebuilding B every step costs m² dot
  products.
- **One-pass step.** Storing y_i = x_i + α·r_i instead of x_i makes the step
  one ZGEMV, x = Y·c.
- **Measured.** One step with 8 stored pairs on 5.5 M complex entries (the
  H2O orbitals on the accurate mesh) takes 0.24 s on 1 thread and 0.19 s on
  2 threads; the previous implementation, which rebuilt B, took 0.78 s. The
  iteration count is the same, and the arithmetic differs only in round-off.
- **Euclidean metric.** Weighting the inner product with the grid metric
  (FEDVR or finite-element weights) does not speed up orbital DIIS: H2O
  Hartree–Fock 12 → 13, Ne Hartree–Fock 12 → 17, He exact exchange 9 → 7
  SCF steps. So the mixer stays Grid-free.

## Common Pitfalls

- **Never feed raw eigenvectors**: remove gauge freedom first
  (`Orbs_AlignOnReference`), otherwise the residual `xNew − x` is
  dominated by arbitrary phases/rotations and any mixing diverges.
- **Mix a complete iterate** (see "What to mix"): DIIS on the Hartree
  potential of Hartree–Fock, a hybrid or MCSCF is slower than linear mixing.
- **Reset between maps**: a program that switches the functional (or the
  system) and continues from the converged orbitals calls `Mixing_Reset`
  (T_FemDft_01, T_FemDft_02); without it the extrapolation mixes pairs of
  two maps (the growth safeguard catches only a converged previous map).
- **Energy-only convergence tests**: a stalled DIIS iteration also has
  |ΔE| → 0. Check the orbital energies too (FemScf.inc: |Δε|), or watch the
  residual with `printLevel: 1`.
- One mixer per run: two quantities of different meaning must not be mixed
  alternately (their pairs would share one history).

## Testing

- `test/componentTests/Mixing/T_Mixing_Diis`:
  - Anderson finite termination: on a 6-dimensional affine map whose plain
    iteration diverges (eigenvalues −1.4 and 1.2), DIIS converges within
    n + 1 = 7 steps;
  - the step after `Mixing_Reset` is the plain step;
  - a new vector length restarts the mixer;
  - `nHistory` 1 reproduces linear mixing (difference 9e-16).
- Simulation tests:
  - DIIS on the orbitals: T_FemDft_02 (Kohn–Sham LDA and PBE0 on grid.fem,
    `Mixing_Reset` between its three SCFs) and T_Be3d_01 (MCSCF gradient);
  - DIIS on the Kohn–Sham potential: T_Ne3d_18;
  - linear mixing: the other ground-solver tests (the default without a
    `mixing` block reproduces the legacy behavior exactly).
