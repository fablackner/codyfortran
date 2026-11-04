# Absorber Module – Agent Onboarding Guide

## Purpose

The **Absorber** module implements boundary-damping schemes for quantum mechanical simulations on finite grids. When propagating wave packets, hard domain boundaries cause spurious reflections that contaminate the physical results. Absorbers mitigate this by applying a smooth mask function M(x) that attenuates wavefunctions near the boundaries:

```
ψ(x) ← M(x) · ψ(x)
```

where M(x) = 1 in the interior (physical region) and M(x) → 0 at the domain edges.

---

## Architecture Overview

The Absorber follows the **Interface Module + Implementation Submodule** pattern used throughout CodyFortran:

```
src/Absorber/
├── M_Absorber.f90              # Public interface (procedure pointers)
├── S_Absorber.f90              # Top-level fabrication dispatcher
├── AGENTS.md                   # This file
├── Linear/
│   ├── M_Absorber_Linear.f90   # Linear-grid family interface
│   ├── S_Absorber_Linear.f90   # Linear-grid variant dispatcher
│   └── Cosinus/
│       ├── M_Absorber_Linear_Cosinus.f90  # Cosine absorber interface
│       └── S_Absorber_Linear_Cosinus.f90  # Cosine absorber implementation
├── Ylm/
│   ├── M_Absorber_Ylm.f90      # Ylm-grid family interface
│   ├── S_Absorber_Ylm.f90      # Ylm-grid variant dispatcher
│   └── Cosinus/
│       ├── M_Absorber_Ylm_Cosinus.f90  # Radial cosine absorber interface
│       └── S_Absorber_Ylm_Cosinus.f90  # Radial cosine absorber implementation
└── Fem/
    ├── M_Absorber_Fem.f90      # FEM-grid family interface
    ├── S_Absorber_Fem.f90      # checks grid.fem, variant dispatcher
    └── Cosinus/
        ├── M_Absorber_Fem_Cosinus.f90  # width, order (exported)
        └── S_Absorber_Fem_Cosinus.f90  # product of per-axis layer profiles
```

### Key Concepts

| Component | Role |
|-----------|------|
| **M_Absorber** | Public API: exports `Absorber_Fabricate`, `Absorber_Setup`, `Absorber_ApplyAbsorber` |
| **Procedure Pointers** | Runtime polymorphism – concrete implementations bind these at startup |
| **Fabricate** | Reads JSON config, selects implementation, binds procedure pointers |
| **Setup** | Precomputes mask arrays (called once after fabrication) |
| **ApplyAbsorber** | Applies mask to orbitals (called each timestep) |

---

## Data Flow

```
┌──────────────────────────────────────────────────────────────────┐
│                        JSON Configuration                         │
│  { "absorber": { "linear": { "cosinus": { ... } } } }            │
└──────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌──────────────────────────────────────────────────────────────────┐
│                    Absorber_Fabricate                             │
│  1. Parse JSON tree                                               │
│  2. Dispatch to Absorber_Linear_Fabricate                        │
│     → Dispatch to Absorber_Linear_Cosinus_Fabricate              │
│  3. Bind: Absorber_Setup → Setup                                 │
│           Absorber_ApplyAbsorber → ApplyAbsorber                 │
└──────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌──────────────────────────────────────────────────────────────────┐
│                       Absorber_Setup                              │
│  1. Allocate maskFunction(Grid_nPoints)                          │
│  2. Evaluate M(x) = cos^(1/n)(π/2 · ξ) at each grid point        │
└──────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌──────────────────────────────────────────────────────────────────┐
│                  Absorber_ApplyAbsorber(orbs)                     │
│  For each orbital: orbs(:,i) ← orbs(:,i) * maskFunction(:)       │
└──────────────────────────────────────────────────────────────────┘
```

---

## JSON Configuration

### Example: Cosine Absorber on Linear Grid

```json
{
  "absorber": {
    "linear": {
      "cosinus": {
        "onset": 80.0,
        "order": 6
      }
    }
  }
}
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `onset` | real | 100.0 | Coordinate where absorption begins (`|x| ≥ onset`) |
| `order` | int | 6 | Exponent denominator (larger = sharper transition) |

### Example: Radial Cosine Absorber on Ylm Grid

```json
{
  "absorber": {
    "ylm": {
      "cosinus": {
        "onset": 16.0,
        "order": 6
      }
    }
  }
}
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `onset` | real | 0.8 · rmax | Radius where absorption begins (`r ≥ onset`) |
| `order` | int | 6 | Exponent denominator (larger = sharper transition) |

The mask depends only on r, so it acts identically on every spherical-harmonic (l,m) channel of the flattened Ylm grid (via `Grid_Ylm_rCoord`). Do not combine with an ECS radial back-end (`fedvrEcs`) — complex scaling already absorbs outgoing flux.

### Example: Boundary-Layer Cosine Absorber on the FEM Grid

```json
{
  "absorber": {
    "fem": {
      "cosinus": {
        "width": 4.0,
        "order": 8
      }
    }
  }
}
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `width` | real | 0.2 · smallest box side | Thickness of the absorbing layer along every face; 0 < width ≤ (smallest side)/2 |
| `order` | int | 6 | Exponent denominator (larger = sharper transition), ≥ 1 |

The `grid.fem` box has homogeneous Dirichlet faces, so without an absorber
outgoing waves reflect there. The mask is a product of one-dimensional layer
profiles (a layer along each of the six faces; edges and corners get the
product):

```
M(x, y, z) = m_x(x) m_y(y) m_z(z),
m_a(s) = cos^(1/order)(π/2 · (width − d_a)/width)   for d_a < width,
       = 1                                           otherwise,
d_a    = min(s − lower_a, upper_a − s)   (distance to the nearer face of axis a)
```

It is 1 in the inner box, decays monotonically toward each face and would
vanish on the faces themselves (Dirichlet nodes, not unknowns). It is built
once in Setup (`Grid_nPoints` values; the default width needs the box, so
`Grid_Fabricate` must run before `Absorber_Fabricate`), and `ApplyAbsorber`
multiplies every orbital column (`size(orbs, 2)`). The family checks that
`grid.fem` exists. `Absorber_Fem_Cosinus_width` and
`Absorber_Fem_Cosinus_order` are exported for printers.

### Orthonormal Orbitals after the Mask (`keepOrthonormalQ`)

```json
{ "absorber": { "keepOrthonormalQ": true, "ylm": { "cosinus": { "onset": 120.0, "order": 4 } } } }
```

The mask on the orbitals is the product mask on the many-body state, so the
masked orbitals lose norm and are no longer orthonormal, which the orbital
equations of motion assume. With `keepOrthonormalQ` (default `false`, family
independent) the time loop calls `Method_Mb_OrbBased_OrthonormalizeKeepingState`
after every mask: Loewdin orthonormalization of the orbitals with the
compensating CI transformation (MCTDHX only, exact in the full CI space;
error stop for the other methods). The many-body state is unchanged; its lost
norm moves into the CI vector.

The 1-RDM is then that of the unnormalized masked state, so one-body
observables are ⟨Ψ|Ô|Ψ⟩ with ⟨Ψ|Ψ⟩ = 1 − P (P: probability that an electron
was absorbed), and N_box = Tr ρ = N ⟨Ψ|Ψ⟩ instead of N minus the absorbed
electrons. Without the option the observables are Tr(ρ O) of the
non-orthonormal orbitals (N_box = Tr(ρ S)). For Ne at t = 198 a.u.
(P ≈ 1.2e-3) the two differ by about that factor: energy and N_box by
⟨Ψ|Ψ⟩, ⟨z⟩ by 5e-4 and the dipole acceleration by 1e-3 relative.

### No Absorber (Passthrough)

If no `absorber` block is present, a no-op implementation is installed that leaves wavefunctions unchanged.

---

## Mathematical Details: Cosine Mask

The cosine-profile mask is defined as:

```
M(x) = 1                              if |x| < onset
     = cos^(1/n)(π/2 · ξ(x))          if |x| ≥ onset
```

where:
- **ξ(x)** = (|x| − onset) / (x_max − onset) ∈ [0, 1]
- **n** = `order` parameter

Properties:
- M = 1 in the interior (no damping)
- M → 0 as x → ±x_max (full absorption at boundaries)
- Smooth transition minimizes spurious reflections
- Larger `order` produces **sharper** transitions (the exponent is 1/n, not n)

---

## Typical Usage in Client Code

```fortran
! Initialization sequence (once at startup)
call Absorber_Fabricate   ! Reads JSON, binds implementation
call Absorber_Setup       ! Precomputes mask array

! Time-stepping loop
do iStep = 1, nSteps
  call Propagator_Propagate(state, t0, t1)
  call Absorber_ApplyAbsorber(Orbs_orbs)  ! Damp outgoing waves
  if (Absorber_keepOrthonormalQ) call Method_Mb_OrbBased_OrthonormalizeKeepingState(Method_state)
end do
```

---

## Adding a New Absorber Implementation

1. **Create directory**: `src/Absorber/Linear/MyAbsorber/` (or new family if not linear)

2. **Interface module** (`M_Absorber_Linear_MyAbsorber.f90`):
   ```fortran
   module M_Absorber_Linear_MyAbsorber
     interface
       module subroutine Absorber_Linear_MyAbsorber_Fabricate
       end subroutine
     end interface
   end module
   ```

3. **Implementation submodule** (`S_Absorber_Linear_MyAbsorber.f90`):
   - Implement local `Setup` and `ApplyAbsorber` subroutines
   - In `Fabricate`, bind to `M_Absorber` procedure pointers

4. **Register in factory**: Add branch in `S_Absorber_Linear.f90`:
   ```fortran
   if (Json_GetExistence("absorber.linear.myabsorber")) then
     call Absorber_Linear_MyAbsorber_Fabricate
   ```

5. **Update CMakeLists.txt** to include new source files

---

## Dependencies

| Depends On | Purpose |
|------------|---------|
| `M_Utils_Types` | R64, I32 type kinds |
| `M_Utils_Json` | Configuration parsing |
| `M_Utils_Say` | Logging/diagnostics |
| `M_Utils_Constants` | PI |
| `M_Utils_NoOpProcedures` | Default no-op implementations |
| `M_Grid`, `M_Grid_Linear`, `M_Grid_Ylm`, `M_Grid_Fem` | Spatial grid info (nPoints, xCoord, rCoord/radialPoints, Fem axes and box) |
| `M_Orbs` | Orbital count (`Orbs_nOrbsInState`) |

---

## Testing Considerations

- Verify mask values: M(x) = 1 for |x| < onset, M(x_max) ≈ 0
- Check norm reduction: after applying absorber, norm should decrease for wavefunctions with amplitude in absorbing region
- Test with different `order` values to verify transition sharpness
- Benchmark ApplyAbsorber performance (should be O(nGrid × nOrbitals))
- `test/componentTests/Absorber/Fem/Cosinus/T_Absorber_Fem_Cosinus`: asymmetric
  box [−6, 6] × [−5, 7] × [−4, 4], width 2, order 4; constant orbitals become
  the mask: exactly 1 at the 27195 inner nodes, the product formula at all
  nodes (difference 0), the explicit profile at the last x node, monotone
  decay toward every face, 0 < M ≤ 1

---

## Common Pitfalls

1. **Forgetting to call Setup**: `ApplyAbsorber` will crash if `maskFunction` is not allocated
2. **Onset > x_max**: Absorber has no effect (entire domain is interior)
3. **Order confusion**: Larger `order` = sharper transition (exponent is 1/n)
4. **Grid mismatch**: Must use same grid type that absorber was configured for
