# Reconstruction Module

> **Note:** The neural-network-based closures documented below (`NnCholesky`,
> `MlpFullCop`, `CnnFullCop`, `CnnCorrelCop`) and `M_Utils_FtorchInference` depend on
> FTorch/libtorch, which is not enabled on `main`. That code lives on the
> `neural-recon` branch until FTorch is enabled here. On `main`, only the
> `Valdemoro` and `Ny` (Nakatsuji-Yasuda) closures are available.

## Purpose

The Reconstruction module implements **BBGKY hierarchy closure** for time-dependent
two-body reduced density matrix (TD-2RDM) dynamics. During propagation of the 2-RDM,
the equation of motion contains terms involving the **3-RDM** (collision operator).
Since tracking the full 3-RDM is computationally expensive, reconstruction schemes
**approximate the 3-RDM from lower-order RDMs** (1-RDM and 2-RDM), enabling tractable
quantum many-body simulations.

Mathematically, reconstruction computes:

```
d(D²)/dt ← ... + Tr₃[Ŵ, D³]
                 ↑
         reconstructed from D¹, D²
```

## Architecture Overview

```
M_Reconstruction (root interface)
│
├── Generic/           ← spin-unresolved implementations
│   └── Valdemoro/     ← free-cumulant approximation C³=0
│
├── Singlet/           ← spin-singlet (S=0) implementations
│   ├── Valdemoro/     ← Valdemoro for α-β blocks
│   ├── Ny/            ← Nakatsuji-Yasuda (cumulant products via R = 2 P_occ - 1)
│   └── NnCholesky/    ← neural-network via Cholesky features
│
└── SingletTr0/        ← trace-free singlet state representation
    ├── Valdemoro/     ← Valdemoro on trace-free 2-RDM
    ├── Ny/            ← Nakatsuji-Yasuda (cumulant products via R = 2 P_occ - 1)
    ├── MlpFullCop/    ← neural-network on packed features
    ├── CnnFullCop/    ← convolutional NN on 4D tensors
    └── CnnCorrelCop/  ← CNN correction to Valdemoro collision operator
```

## Key Concepts

### Procedure Pointer Binding

The module follows the CodyFortran **fabrication pattern**:

1. `M_Reconstruction.f90` declares the **public interface** with procedure pointers
2. `S_Reconstruction.f90` reads JSON config and delegates to a family (Generic/Singlet/SingletTr0)
3. Each family reads its JSON section and binds the final implementation

```fortran
! At runtime, JSON config selects which implementation backs this pointer:
Reconstruction_AddCollisionOp => AddCollisionOp   ! bound in submodule
```

### State Representations

| Family      | 2-RDM Layout                | 3-RDM Blocks          |
|-------------|-----------------------------|-----------------------|
| Generic     | Full nOrbs⁴ tensor          | Full nOrbs⁶ tensor    |
| Singlet     | α-β block rdm2ab(nα,nα,nβ,nβ) | rdm3aab(nα³×nβ³)    |
| SingletTr0  | Trace-free α-β + explicit 1-RDM | Same as Singlet    |

### Reconstruction Approximations

| Scheme        | Method | Description |
|---------------|--------|-------------|
| **Valdemoro** | Classical | Sets 3-cumulant C³=0, reconstructs D³ from free-particle products |
| **Nakatsuji-Yasuda (NY)** | Classical | Adds particle-number-dependent corrections to C³ |
| **NnCholesky** | ML | Neural network maps Cholesky-decomposed 2-RDM → 3-RDM Cholesky vectors |
| **MlpFullCop** | ML | Neural network maps packed trace-free 2-RDM → collision operator |
| **CnnFullCop** | ML | Convolutional NN processes local periodic blocks with real/imag channel axis |
| **CnnCorrelCop** | ML | CNN predicts correction to Valdemoro collision operator from the 2-RDM cumulant |

## Public Interface

```fortran
! From M_Reconstruction:
call Reconstruction_Fabricate       ! Bind implementation from JSON
call Reconstruction_Setup           ! Allocate workspace (if needed)
call Reconstruction_AddCollisionOp(dD2, D2, h2)  ! Add [Ŵ,D³] to dD²/dt
```

### `Reconstruction_AddCollisionOp(dTwoRdmState, twoRdmState, h2)`

**Core routine** called during each time derivative evaluation:

| Argument        | Type                  | Description |
|-----------------|-----------------------|-------------|
| `dTwoRdmState`  | `complex(R64)(:)` INOUT | Packed derivative vector; collision term **added** in-place |
| `twoRdmState`   | `complex(R64)(:)` IN   | Packed 2-RDM state vector |
| `h2`            | `complex(R64)(:,:,:,:)` IN | Two-body Hamiltonian tensor |

## JSON Configuration

```json
{
  "reconstruction": {
    "generic": {
      "genericValdemoro": {}
    }
  }
}
```

```json
{
  "reconstruction": {
    "singlet": {
      "valdemoro": { "ccQ": false },
      "singletNy": { "ccQ": true },
      "nnCholesky": { "modelFile": "model.pt" }
    }
  }
}
```

```json
{
  "reconstruction": {
    "singletTr0": {
      "valdemoro": {},
      "mlpFullCop": { "modelFile": "model.pt", "natOrbTransQ": true },
      "cnnFullCop": {
        "modelFile": "model.pt",
        "r": 3
      },
      "cnnCorrelCop": {
        "modelFile": "model.pt",
        "r": 3
      }
    }
  }
}
```

### Configuration Options

| Path | Type | Default | Description |
|------|------|---------|-------------|
| `*.valdemoro.ccQ` | bool | false | Apply contraction-consistent trace constraints |
| `*.nnCholesky.modelFile` | string | "model.pt" | PyTorch model file path |
| `*.mlpFullCop.natOrbTransQ` | bool | false | Transform to natural orbital basis before inference |
| `*.cnnFullCop.modelFile` | string | "model.pt" | Single TorchScript model for real/imag channel input |
| `*.cnnFullCop.r` | int | 3 | Local periodic block size (must be odd) |
| `*.cnnCorrelCop.modelFile` | string | "model.pt" | TorchScript model for Valdemoro correction |
| `*.cnnCorrelCop.r` | int | 3 | Local periodic block size (must be odd) |

## Data Flow in `AddCollisionOp`

```
twoRdmState (packed)
       │
       ▼
┌──────────────────┐
│ Unpack to rdm2ab │
└──────────────────┘
       │
       ▼
┌──────────────────┐
│ Trace → rdm1a    │
└──────────────────┘
       │
       ▼
┌──────────────────────────────────┐
│ Reconstruct rdm3aab from         │
│   rdm1a, rdm2ab (or via NN)      │
└──────────────────────────────────┘
       │
       ▼
┌───────────────────────────────────┐
│ RdmCollisionOp_FillCollisionOp2ab │
│   collOp2ab = Tr₃[h2, rdm3aab]    │
└───────────────────────────────────┘
       │
       ▼
┌──────────────────────────────────┐
│ (Optional) Enforce trace / project│
│   for trace-free representations  │
└──────────────────────────────────┘
       │
       ▼
dRdm2ab += collOp2ab
       │
       ▼
Pack back → dTwoRdmState
```

## Adding a New Reconstruction Scheme

1. **Create directory** under appropriate family:
   ```
   src/Reconstruction/{Family}/MyScheme/
   ```

2. **Create module file** `M_Reconstruction_{Family}_MyScheme.f90`:
   ```fortran
   module M_Reconstruction_{Family}_MyScheme
     interface
       module subroutine Reconstruction_{Family}_MyScheme_Fabricate
       end subroutine
     end interface
   end module
   ```

3. **Create submodule** `S_Reconstruction_{Family}_MyScheme.f90`:
   ```fortran
   submodule(M_Reconstruction_{Family}_MyScheme) S_Reconstruction_{Family}_MyScheme
   contains
     module subroutine Reconstruction_{Family}_MyScheme_Fabricate
       use M_Reconstruction
       Reconstruction_AddCollisionOp => AddCollisionOp
     end subroutine

     subroutine AddCollisionOp(dTwoRdmState, twoRdmState, h2)
       ! Your implementation
     end subroutine
   end submodule
   ```

4. **Register in parent factory** (e.g., `S_Reconstruction_Singlet.f90`):
   ```fortran
   else if (Json_GetExistence("reconstruction.singlet.myScheme")) then
     call Reconstruction_Singlet_MyScheme_Fabricate
   ```

5. **Update CMakeLists.txt** to include new files

## Dependencies

| Module | Purpose |
|--------|---------|
| `M_Utils_RdmCumulants` | Cumulant expansion (free-part, NY corrections) |
| `M_Utils_RdmCollisionOp` | Compute Tr₃[Ŵ, D³] |
| `M_Utils_RdmTrace` | Extract 1-RDM from 2-RDM |
| `M_Utils_RdmEnforceTraceM2ab` | Project 2-RDM to trace-free subspace |
| `M_Utils_RdmNatOrbTransform` | Natural orbital transformation |
| `M_Utils_FtorchInference` | PyTorch model inference via FTorch |
| `M_TwoRdm_*` | Pack/unpack state vectors |
| `M_Method_Mb` | Access particle numbers |
| `M_Method_Mb_OrbBased` | Access orbital counts |

## Testing

Reconstruction is exercised in TD-2RDM propagation tests:
```
test/T_He1d_06_TimePropagationTd2rdm.f90
test/T_He1d_06_TimePropagationTd2rdm.json
```

## References

- **Valdemoro approximation**: C. Valdemoro, Phys. Rev. A 45, 4462 (1992)
- **Nakatsuji-Yasuda**: H. Nakatsuji & K. Yasuda, Phys. Rev. Lett. 76, 1039 (1996)
- **Neural network closures**: Project-specific ML models trained on exact RDM data
