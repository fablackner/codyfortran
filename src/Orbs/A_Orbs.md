# Orbs Module — Agent Onboarding Guide

## Overview

The **Orbs** module manages single-particle orbital coefficients in CodyFortranRDM.
It provides the orbital part of the many-body wavefunction ansatz used in MCTDH-type
methods (MCTDHF, MCTDHB, TD-2RDM, etc.).

**Key responsibility**: Store and manipulate the orbital coefficient matrix
`Orbs_orbs(nBasis, nOrbs)` where each column represents one orbital expanded in
the spatial basis provided by the Grid module.

---

## Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                        M_Orbs.f90                               │
│  Public API (interface module)                                  │
│  • Orbs_orbs            — orbital coefficient matrix (pointer)  │
│  • Orbs_nOrbsInState    — number of active orbitals             │
│  • Orbs_restrictedQ     — spin-restricted flag                  │
│  • Orbs_Setup           — procedure pointer → backend           │
│  • Orbs_Orthonormalize  — procedure pointer → backend           │
│  • Orbs_ProjectOnSubspace — procedure pointer → backend         │
│  • Orbs_SaveOrbs        — procedure pointer → backend           │
└─────────────────────────────────────────────────────────────────┘
                              │
                              │ submodule implements
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                        S_Orbs.f90                               │
│  Implementation (submodule)                                     │
│  • Orbs_Fabricate — binds procedure pointers, reads JSON        │
│  • Orthonormalize — per-body-type Gram-Schmidt via Grid         │
│  • ProjectOnSubspace — gauge enforcement via Grid               │
│  • SaveOrbs — binary file I/O per orbital                       │
└─────────────────────────────────────────────────────────────────┘
                              │
                              │ delegates to
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                        M_Grid                                   │
│  Spatial discretization layer                                   │
│  • Grid_Orthonormalize  — metric-aware Gram-Schmidt             │
│  • Grid_ProjectOnSubspace — metric-aware projection             │
│  • Grid_nPoints         — basis dimension                       │
└─────────────────────────────────────────────────────────────────┘
```

---

## Data Layout

### Orbital Matrix

```fortran
complex(R64), pointer :: Orbs_orbs(:, :)
! Shape: (Grid_nPoints, Orbs_nOrbsInState)
! Convention: Orbs_orbs(i, j) = j-th orbital's coefficient at grid point i
```

### Body-Type Partitioning

Orbitals are partitioned by **body type** (particle species):

```
Orbital index:   1  2  3  4  5  6  7  8
                 |--BT1--|  |---BT2---|
                 nOrbs(1)=3   nOrbs(2)=5
```

Key variables from `M_Method_Mb_OrbBased`:
- `nOrbs(ibt)` — number of orbitals for body type `ibt`
- `nOrbsSum` — total orbitals = sum(nOrbs)
- `nOrbsStart(ibt)`, `nOrbsEnd(ibt)` — index range for body type

### Restricted vs Unrestricted

| Mode         | `Orbs_restrictedQ` | `Orbs_nOrbsInState`  | Use case               |
|--------------|-------------------|----------------------|------------------------|
| Unrestricted | `.false.`         | `nOrbsSum`           | General (default)      |
| Restricted   | `.true.`          | `nOrbsSum / 2`       | Closed-shell systems   |

In restricted mode, spin-up and spin-down electrons share the same spatial
orbitals, reducing computational cost and exploiting spin symmetry.

---

## Lifecycle

```fortran
! 1. Fabrication phase (bind procedure pointers)
call Grid_Fabricate()
call Method_Mb_OrbBased_Fabricate()
call Orbs_Fabricate()           ! Reads orbs.restrictedQ, sets nOrbsInState
call OrbsInit_Fabricate()

! 2. Setup phase (allocate and initialize)
call Grid_Setup()
call Method_Mb_OrbBased_Setup()
call OrbsInit_Setup()           ! Allocates Orbs_orbs, calls OrbsInit_Initialize

! 3. Runtime operations
call Orbs_Orthonormalize(Orbs_orbs)           ! Ensure orthonormality
call Orbs_ProjectOnSubspace(dOrbs, Orbs_orbs) ! Gauge enforcement
call Orbs_SaveOrbs(Orbs_orbs)                 ! Checkpoint to disk
```

---

## Key Operations

### Orthonormalization

**Purpose**: Ensure `<φ_i | φ_j> = δ_ij` within each body type.

**Algorithm**: Modified Gram-Schmidt with grid metric (handles non-uniform
grids, FEDVR elements, spherical coordinates, etc.).

**Why per-body-type?** Different particle species are distinguishable; their
orbitals need not be mutually orthogonal.

### Subspace Projection

**Purpose**: Remove components of orbital derivatives parallel to the current
orbital space.

**Math**: `dφ' = dφ - Σ_j <φ_j|dφ> φ_j`

**Physical meaning**: Enforces the MCTDH gauge condition `<φ_j|dφ_i/dt> = 0`,
which eliminates redundant parametrization in the variational equations.

### File I/O

**Format**: One binary file per orbital, named `orbBB_II.in`
- `BB`: body type (01, 02, ...)
- `II`: orbital index within body type (01, 02, ...)

**Content**: Raw `complex(R64)` array of length `Grid_nPoints`.

---

## JSON Configuration

```json
{
  "orbs": {
    "restrictedQ": false
  }
}
```

| Key               | Type    | Default | Description                          |
|-------------------|---------|---------|--------------------------------------|
| `orbs.restrictedQ`| logical | false   | Use spin-restricted representation   |

---

## Integration Points

### Upstream Dependencies
- **M_Grid**: Spatial basis and metric for inner products
- **M_Method_Mb**: Number of body types (`nBodyTypes`)
- **M_Method_Mb_OrbBased**: Orbital partitioning (`nOrbs`, `nOrbsSum`)
- **M_Utils_Json**: Configuration parsing

### Downstream Consumers
- **M_OrbsInit**: Initializes orbital coefficients
- **M_Method_Mb_OrbBased**: Computes Hamiltonian matrix elements using orbitals
- **M_Propagator**: Time-evolves the orbital+coefficient state

---

## Common Patterns

### Body-Type Loop

All orbital operations use this pattern to process each species independently:

```fortran
integer(I32) :: ibt, startOrb, endOrb

startOrb = 1
do ibt = 1, Method_Mb_nBodyTypes
  endOrb = startOrb + Method_Mb_OrbBased_nOrbs(ibt) - 1
  
  ! Process orbs(:, startOrb:endOrb) for body type ibt
  call SomeOperation(orbs(:, startOrb:endOrb))
  
  startOrb = endOrb + 1
end do
```

### State Packing

In propagation methods, orbitals are packed into a flat state vector:

```fortran
state(1 : nBasis*nOrbsInState) = reshape(Orbs_orbs, [nBasis*nOrbsInState])
```

---

## Testing

Run Orbs-related tests:
```bash
cd build && ctest -R T_.*Orb
```

Relevant test files demonstrate:
- Orbital initialization (Harmonic, Hydrogen-like, Lattice)
- MCTDH propagation with orbital dynamics
- Restricted vs unrestricted configurations

---

## Troubleshooting

| Symptom | Likely Cause | Solution |
|---------|--------------|----------|
| Orbitals not orthonormal | Missing `Orbs_Orthonormalize` call | Call after initialization and periodically during propagation |
| Gauge drift in MCTDH | Missing `ProjectOnSubspace` | Ensure projection is called in time derivative |
| Wrong `nOrbsInState` | `restrictedQ` mismatch | Check JSON config matches physical system |
| File not found on restart | Wrong body type indexing | Verify `orbBB_II.in` naming convention |

---

## Extension Points

To add a new orbital backend:
1. Create a new submodule `S_Orbs_MyBackend.f90`
2. Implement the procedure interfaces from `M_Orbs`
3. Add fabrication logic to select the backend based on JSON configuration
4. Update `Orbs_Fabricate` to bind to the new implementation

The module/submodule architecture allows multiple backends to coexist without
modifying the public API.
