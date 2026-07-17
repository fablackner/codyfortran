# SysGauge Module — Onboarding Guide

## Overview

The **SysGauge** module provides the light–matter coupling term of the
Hamiltonian in **velocity gauge** and dipole approximation:

```
Ĥ_A = A(t)·p̂ / m + A(t)² / (2m)
```

for an electron (charge −1, atomic units) in a linearly z-polarized laser
pulse. It complements the other Hamiltonian-term modules:

```
Ĥ = T̂ (SysKinetic) + V̂_ext (SysPotential) + Ĥ_A (SysGauge) + Ŵ (SysInteraction)
```

**Why a separate module?** The gauge coupling contains the momentum operator
and is therefore *not diagonal in position space*. It cannot satisfy the
SysPotential contract (`FillExternalPotential` + point-wise
`MultiplyWithExternalPotential`), which is the contract of multiplicative
potentials. Instead SysGauge follows the SysKinetic operator contract: a
single application routine `SysGauge_MultiplyWithGaugeOp(dOrb, orb, time, bt_)`.
Length-gauge couplings (E(t)·z) remain multiplicative and belong to
SysPotential (see `SysPotential/Ylm/CoulombLaser`).

The uniform A²/(2m) term only contributes a global time-dependent phase; it
is included so that the operator is exactly (p̂ + A)²/2m − p̂²/2m.

## Architecture

The module follows the CodyFortranRDM fabrication pattern:

```
┌────────────────────────────────────────────────────────────────┐
│                     M_SysGauge (Facade)                        │
│  - SysGauge_Fabricate            → runtime configuration       │
│  - SysGauge_Setup                → capture grid data (pointer) │
│  - SysGauge_MultiplyWithGaugeOp  → apply Ĥ_A·ψ (pointer)       │
├────────────────────────────────────────────────────────────────┤
│                          Ylm                                   │
│            (spherical harmonics, z-polarized)                  │
│                           │                                    │
│                       Velocity                                 │
│         (sin²-envelope pulse, Δl = ±1 channel coupling)        │
│                      │           │                             │
│                    Fedvr      FedvrEcs                         │
│              (radial d/dr)  (d/dz on ECS contour)              │
└────────────────────────────────────────────────────────────────┘
```

## File Structure

```
SysGauge/
├── M_SysGauge.f90            # Public facade: pointers, flags, Fabricate
├── S_SysGauge.f90            # Dispatch to grid backends (Ylm)
├── A_SysGauge.md             # This file
└── Ylm/
    ├── M_SysGauge_Ylm.f90    # Radial building blocks: first-derivative
    │                         #   pointer + (complex) radial coordinates
    ├── S_SysGauge_Ylm.f90    # Dispatch to gauge models (velocity)
    └── Velocity/
        ├── M_SysGauge_Ylm_Velocity.f90   # Pulse parameters, A(t)
        ├── S_SysGauge_Ylm_Velocity.f90   # Channel-coupling operator
        ├── Fedvr/                        # FEDVR radial derivative
        └── FedvrEcs/                     # FEDVR-ECS contour derivative
```

## Physics Notes

### Pulse Definition

The velocity-gauge pulse is specified through the **vector potential**
(the natural quantity in this gauge):

```
A(t) = -(E0/ω) · sin²(π (t−t0)/T) · sin(ω (t−t0−T/2) + cep),  T = 2π nCycles/ω
```

zero outside [t0, t0+T]. A(t) vanishes identically at both pulse edges. The
physical electric field E(t) = −dA/dt equals the CoulombLaser length-gauge
field `E0 sin² cos` **plus an envelope-derivative correction of relative size
1/nCycles** — observables from the two gauges therefore agree only in the
many-cycle limit; they are not term-by-term identical pulses.

### Ylm Channel Coupling

With ψ = Σ_{lm} f_{lm}(r) Y_l^m, a z-polarized A couples Δl = ±1, Δm = 0:

```
∂_z [f_l Y_{lm}] = c_{l,m} (∂_r − l/r) f_l · Y_{l+1,m}
                 + c_{l−1,m} (∂_r + (l+1)/r) f_l · Y_{l−1,m}

c_{l,m} = sqrt( ((l+1)² − m²) / ((2l+1)(2l+3)) )
```

**Implementation trick** (same as SysKinetic): with g = r·f,

```
∂_r − l/r      → (1/r)(d/dr − (l+1)/r) g
∂_r + (l+1)/r  → (1/r)(d/dr + l/r) g
```

The r² metric factor is absorbed by the transform, the weak-form FEDVR first
derivative acting on g is exactly antisymmetric under the plain radial
weights, and the ±k/r pieces become exactly Hermitian diagonal terms. The
discrete Ĥ_A is Hermitian in the FEDVR-weighted metric to machine precision
— required for norm conservation with the SIL integrator in the split-step
propagator — up to a boundary term ∝ g_bra(rmax)·g_ket(rmax) (the radial
grid includes the rmax endpoint), which is negligible for bounded or
absorbed wavefunctions. One radial derivative per channel serves both the
raising and lowering terms.

On the ECS contour, derivatives and the ±k/z terms use the complex contour
points, matching the analytic continuation of the ECS kinetic operator
(c-product metric).

## Usage

### Initialization Sequence

```fortran
call SysGauge_Fabricate    ! Reads JSON, binds procedure pointers
...
call Grid_Setup            ! Grid first: backends capture radial data
call SysGauge_Setup        ! Copies radial/contour coordinates
```

### JSON Configuration

```json
{
  "sysGauge": {
    "ylm": {
      "velocity": {
        "fieldStrength": 0.1,
        "omega": 0.5,
        "nCycles": 1.0,
        "cep": 0.0,
        "tStart": 0.0,
        "bodyMass": [ 1.0 ],
        "fedvr": { }
      }
    }
  }
}
```

`fedvrEcs: { }` selects the ECS contour derivative instead (requires
`grid.ylm.fedvrEcs`).

### Consumption by Method

`Method_Mb_OrbBased` binds `Method_Mb_OrbBased_ApplyGaugeOp` when the JSON
contains `sysGauge` (no-op otherwise) and applies it inside
`ApplySingleBodyOp` and `FillH1`, so the coupling automatically enters:

- TDHX/MCTDHX time derivatives and the h1 matrix (energies, CI part)
- the *linear* part of the split-step propagator (integrated by SIL — the
  gauge term is stiff like the kinetic term and belongs there)
- ground-state solvers via h1 (a pulse with tStart ≥ 0 has A(0) = 0, so
  ground states are unaffected)

## Exported Symbols

| Symbol | Type | Description |
|--------|------|-------------|
| `SysGauge_timeIndependentQ` | `logical` | Always `.false.` for a pulse |
| `SysGauge_bodyTypeIndependentQ` | `logical` | Mass-independence flag |
| `SysGauge_Setup` | pointer | Capture grid-dependent data |
| `SysGauge_MultiplyWithGaugeOp` | pointer | `(dOrb, orb, time, bt_)`: dOrb = Ĥ_A·orb |
| `SysGauge_Ylm_ApplyRadialFirstDerivative` | pointer | `(dfLm, fLm)` per-channel d/dr |
| `SysGauge_Ylm_radialCoordinates` | `complex(:)` | r (Fedvr) or z(r) (ECS) |
| `SysGauge_Ylm_Velocity_VectorPotentialAmplitude` | function | A(t) of the pulse |

## Tests

- `test/componentTests/SysGauge/Ylm/Velocity/Fedvr/T_SysGauge_Ylm_Velocity_Fedvr`:
  analytic action on a (0,0) channel, selection rules, Hermiticity in the
  weighted metric (machine precision).
- `test/simulationTests/Ne3d/T_Ne3d_08_TimePropagationMctdhxLaserVelocity`:
  velocity-gauge companion of the length-gauge T_Ne3d_07; checks the
  ground state is unaffected (A(0) = 0), norm conservation through SIL,
  and a dipole regression value.

## Adding a New Backend

- **New grid** (e.g., Linear): add `SysGauge/Linear/...` with the 1D coupling
  A(t) p_x = −i A ∂_x, dispatch in `S_SysGauge.f90`. The 1D case needs no
  channel coupling — only a first-derivative application.
- **New radial discretization**: mirror `Ylm/Velocity/Fedvr/`, bind
  `SysGauge_Ylm_ApplyRadialFirstDerivative` and fill
  `SysGauge_Ylm_radialCoordinates` in Setup.
- **New pulse shape**: extend the Velocity model or add a sibling model under
  `Ylm/` and dispatch in `S_SysGauge_Ylm.f90`.

## Dependencies

- `M_Grid_Ylm`, `M_Grid_Ylm_Fedvr`, `M_Grid_Ylm_FedvrEcs`: channel layout,
  derivative contexts, radial/contour points
- `M_Utils_DerivativeFedvr`, `M_Utils_DerivativeFedvrEcs`: first derivative
- `M_Utils_Json`, `M_Utils_Say`, `M_Utils_Types`, `M_Utils_Constants`
