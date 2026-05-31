# Absorber Module – Agent Onboarding Guide

## Purpose

The **Absorber** module implements boundary-damping schemes for quantum mechanical simulations on finite grids. When propagating wave packets, hard domain boundaries cause spurious reflections that contaminate the physical results. Absorbers mitigate this by applying a smooth mask function M(x) that attenuates wavefunctions near the boundaries:

```
ψ(x) ← M(x) · ψ(x)
```

where M(x) = 1 in the interior (physical region) and M(x) → 0 at the domain edges.

---

## Architecture Overview

The Absorber follows the **Interface Module + Implementation Submodule** pattern used throughout CodyFortranRDM:

```
src/Absorber/
├── M_Absorber.f90              # Public interface (procedure pointers)
├── S_Absorber.f90              # Top-level fabrication dispatcher
├── AGENTS.md                   # This file
└── Linear/
    ├── M_Absorber_Linear.f90   # Linear-grid family interface
    ├── S_Absorber_Linear.f90   # Linear-grid variant dispatcher
    └── Cosinus/
        ├── M_Absorber_Linear_Cosinus.f90  # Cosine absorber interface
        └── S_Absorber_Linear_Cosinus.f90  # Cosine absorber implementation
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
| `M_Grid`, `M_Grid_Linear` | Spatial grid info (nPoints, xCoord) |
| `M_Orbs` | Orbital count (`Orbs_nOrbsInState`) |

---

## Testing Considerations

- Verify mask values: M(x) = 1 for |x| < onset, M(x_max) ≈ 0
- Check norm reduction: after applying absorber, norm should decrease for wavefunctions with amplitude in absorbing region
- Test with different `order` values to verify transition sharpness
- Benchmark ApplyAbsorber performance (should be O(nGrid × nOrbitals))

---

## Common Pitfalls

1. **Forgetting to call Setup**: `ApplyAbsorber` will crash if `maskFunction` is not allocated
2. **Onset > x_max**: Absorber has no effect (entire domain is interior)
3. **Order confusion**: Larger `order` = sharper transition (exponent is 1/n)
4. **Grid mismatch**: Must use same grid type that absorber was configured for
