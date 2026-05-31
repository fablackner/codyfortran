# GroundSolver Module

## Overview

The GroundSolver module implements **self-consistent field (SCF) ground-state solvers** for quantum many-body systems. It finds the ground state of a Hamiltonian by iteratively diagonalizing an effective single-particle (mean-field) operator until convergence.

Currently, the only available method is **TDHx** (Time-Dependent Hartree–Exchange), a Hartree–Fock mean-field approach.

## Architecture

```
M_GroundSolver.f90              ← Interface: procedure pointers + abstract interfaces
└── S_GroundSolver.f90          ← Dispatch: JSON routing to backends
    └── Tdhx/
        ├── M_GroundSolver_Tdhx.f90      ← TDHx interface
        ├── S_GroundSolver_Tdhx.f90      ← TDHx dispatch
        ├── StdImpl/
        │   ├── M_GroundSolver_Tdhx_StdImpl.f90  ← Standard impl interface
        │   └── S_GroundSolver_Tdhx_StdImpl.f90  ← Standard impl (general grids)
        └── Ylm/
            ├── M_GroundSolver_Tdhx_Ylm.f90     ← Ylm-optimized interface
            └── S_GroundSolver_Tdhx_Ylm.f90     ← Ylm-optimized impl (spherical grids)
```

The module follows the **CodyFortranRDM fabrication pattern**:
- **M_* files**: Declare procedure pointers and abstract interfaces (contracts)
- **S_* files**: Provide implementations and JSON-based dispatch logic

## Key Concepts

### SCF Iteration

Each `GroundSolver_Approach` call performs one SCF cycle:

1. **Build mean-field potentials** from current orbitals
   - Hartree (direct): Ĵ = ∑_j ∫|φ_j|²/|r−r'| dr'
   - External: V̂_ext (e.g., nuclear Coulomb −Z/r)

2. **Diagonalize the Fock operator**
   - F̂ = T̂ + V̂_ext + Ĵ − K̂
   - Uses iterative ARPACK solver via DiagonalizerList

3. **Mix eigenvectors with old orbitals**
   - φ_new = (1−α)·φ_old + α·evec
   - α < 1 damps oscillations for difficult convergence

4. **Orthonormalize** via Gram–Schmidt

The caller controls the convergence loop (checking energy change).

### TDHx Backends

| Backend | Use Case | Grid | Diagonalizer |
|---------|----------|------|--------------|
| **stdImpl** | General 1D/2D grids | Any | Single nPoints×nPoints |
| **ylmOpt** | Spherical atoms | Ylm | (lmax+1) × nRadial×nRadial |

**ylmOpt** exploits angular symmetry: instead of one large matrix, it diagonalizes smaller per-l-channel matrices. This requires the user to set up multiple diagonalizers with l-specific wrappers.

## Public Interface

```fortran
! Bind backend from JSON configuration
call GroundSolver_Fabricate()

! Allocate working arrays
call GroundSolver_Setup()

! Perform one SCF iteration
call GroundSolver_Approach(state, alpha, time)
```

### Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `state` | `complex(R64)(:)` | Packed orbital coefficients (inout) |
| `alpha` | `real(R64)` | Mixing parameter (0 < α ≤ 1) |
| `time` | `real(R64)` | Time for potentials (usually 0) |

## JSON Configuration

### Standard Implementation (stdImpl)

```json
{
    "groundSolver": {
        "tdhx": {
            "stdImpl": { }
        }
    },
    "diagonalizerList": {
        "arpack": {
            "nEvals": 5,
            "which": "SR",
            "nKry": 15
        }
    }
}
```

### Ylm-Optimized Implementation (ylmOpt)

```json
{
    "groundSolver": {
        "tdhx": {
            "ylmOpt": { }
        }
    },
    "diagonalizerList": {
        "arpack1": { "nEvals": 2, "which": "SR" },
        "arpack2": { "nEvals": 1, "which": "SR" }
    }
}
```

For ylmOpt, you need one diagonalizer entry per l-channel (l=0, l=1, ...).

## Usage Example

```fortran
program scf_example
  use M_GroundSolver
  use M_Method
  
  real(R64) :: energy_old, energy_new, conv
  real(R64) :: alpha = 0.5_R64
  real(R64) :: conv_thresh = 1e-12_R64
  
  ! ... setup Grid, SysKinetic, SysPotential, etc. ...
  
  call GroundSolver_Fabricate()
  call GroundSolver_Setup()
  
  ! SCF convergence loop
  energy_new = 1e10_R64
  do while (abs(conv) > conv_thresh)
    call GroundSolver_Approach(Method_state, alpha, time=0.0_R64)
    energy_old = energy_new
    energy_new = Method_GetEnergy(0.0_R64)
    conv = energy_new - energy_old
  end do
end program
```

## Dependencies

### Required Modules

- **Grid**: Spatial discretization, `Grid_nPoints`
- **SysKinetic**: Kinetic energy operator, `MultiplyWithKineticOp`
- **SysPotential**: External potential, `FillExternalPotential`
- **SysInteraction**: Two-body interaction, `FillInteractionSrc/Potential`
- **Orbs**: Orbital container, `Orbs_Orthonormalize`
- **DiagonalizerList**: Iterative eigensolvers (ARPACK)

### Ylm-Specific Dependencies

- **Grid_Ylm**: `nRadial`, `lmax`, `SetLmComponent`, `GetLmComponent`
- **SysKinetic_Ylm**: `MultiplyWithRadialKineticOp`
- **SysPotential_Ylm**: `FillExternalPotentialRadial`
- **SysInteraction_Ylm**: `FillInteractionPotentialRadial`
- **OrbsInit_Ylm_HydrogenLike**: Quantum numbers `n(j)`, `l(j)`, `m(j)`

## Extending the Module

### Adding a New Backend

1. Create `M_GroundSolver_NewBackend.f90` with interface + fabricate declaration
2. Create `S_GroundSolver_NewBackend.f90` with implementation
3. Add dispatch branch in `S_GroundSolver.f90`:

```fortran
if (Json_GetExistence("groundSolver.newBackend")) then
  call GroundSolver_NewBackend_Fabricate
```

### Implementation Checklist

Your implementation must bind these procedure pointers:
- `GroundSolver_Setup` → Allocate working arrays
- `GroundSolver_Approach` → One SCF iteration
- (Optional) Backend-specific action pointer for the diagonalizer

## Tests

| Test | Description |
|------|-------------|
| `T_Ne3d_02_GroundSolverTdhx` | Neon atom (stdImpl), E ≈ −127.74 Ha |
| `T_Ne3d_03_GroundSolverTdhxYlm` | Neon atom (ylmOpt), E ≈ −128.55 Ha |

Run with: `ctest -R T_Ne3d_0[23]`

## Physical Background

The TDHx method solves the **Hartree–Fock equations**:

```
F̂ φ_i = ε_i φ_i
```

where the Fock operator is:

```
F̂ = T̂ + V̂_ext + Ĵ − K̂
```

- **T̂**: Kinetic energy (−½∇²)
- **V̂_ext**: External potential (e.g., −Z/r for atoms)
- **Ĵ**: Hartree (direct) potential, electron repulsion mean-field
- **K̂**: Fock (exchange) operator, antisymmetry requirement

The self-consistency arises because Ĵ and K̂ depend on the orbitals {φ_i} being solved for.

## Known Limitations

1. **Restricted calculations only**: Spin-up and spin-down orbitals are assumed identical
2. **No DIIS acceleration**: Simple linear mixing only; may converge slowly for some systems
3. **ylmOpt assumes m=0 trials**: The Fock action uses m=0 in the diagonalizer callbacks

## Troubleshooting

| Issue | Solution |
|-------|----------|
| Slow/no convergence | Reduce α (e.g., 0.3–0.5), increase `nTimeSteps` |
| ARPACK convergence failure | Increase `nKry`, check `nEvals` ≤ dim/2 |
| Wrong energy (ylmOpt) | Verify DiagonalizerList order matches l=0,1,2,... |
