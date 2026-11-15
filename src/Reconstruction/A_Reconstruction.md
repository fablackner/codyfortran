# Reconstruction Module

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
│   ├── Triples/       ← optional extension: connected triples (own A_*.md)
│   └── Anchor/        ← optional extension of the closure (own A_*.md)
│
└── SingletTr0/        ← trace-free singlet state representation
    ├── Valdemoro/     ← Valdemoro on trace-free 2-RDM
    └── Ny/            ← Nakatsuji-Yasuda (cumulant products via R = 2 P_occ - 1)
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
| **Nakatsuji-Yasuda (NY)** | Classical | 3-cumulant from products of two 2-cumulants contracted through R = 2 P_occ - 1 |

### Singlet Driver

A singlet scheme fabricate binds only the closure
`Reconstruction_Singlet_FillRdm3aab` (`M_Reconstruction_Singlet`: every
element of the ket-first rdm3aab from rdm1a and rdm2ab, no contraction
consistency). The family (`S_Reconstruction_Singlet`) reads `<scheme>.ccQ`
and binds `Reconstruction_AddCollisionOp`, which

1. unpacks rdm2ab and forms rdm1a;
2. reconstructs rdm3aab = R(rdm2ab) with the closure of the scheme;
3. enforces the contraction consistency (ccQ, `RdmEnforceTraceM3aab_Apply`);
4. adds Tr3[W, D3] (`RdmCollisionOp_FillCollisionOp2ab`, one ZGEMM);
5. packs the derivative.

Steps 2 and 3 are `Reconstruction_Singlet_ReconstructRdm3aab`. Extensions
wrap the closure pointer (and, if needed, the collision operator) instead of
changing the driver, each one wrapping the closure bound before it:

- `Singlet/Triples/` (`triples`): the connected triples block of the
  3-cumulant from the 2-cumulant, the interaction and the orbital energies
  (`moellerPlesset`). It needs the Hamiltonian
  (`Reconstruction_SetHamiltonian`). It improves NY for C2H4 and not for Ne
  and N2.
- `Singlet/Anchor/` (`anchor`): the closure error at a reference state.

The Generic and SingletTr0 schemes bind their own AddCollisionOp. Cost per
derivative (NY with ccQ, static orbitals, ethylene, 14 orbitals, 8 threads):
0.23 s, 0.18 s without ccQ. SingletTr0 has no ccQ; for its accuracy see the
ethylene section of `app/README.md`.

## Public Interface

```fortran
! From M_Reconstruction:
call Reconstruction_Fabricate       ! Bind implementation from JSON
call Reconstruction_Setup           ! Workspace, reference data (after Method_Setup)
call Reconstruction_AddCollisionOp(dD2, D2, h2)  ! Add [Ŵ,D³] to dD²/dt
call Reconstruction_SetOrbitals(orbs)  ! orbitals of the 2-RDM (moving-orbital methods; default no-op)
call Reconstruction_SetHamiltonian(h1, h2)  ! Hamiltonian of the next AddCollisionOp (td2rdmStaticOrbs; default no-op)
```

### `Reconstruction_AddCollisionOp(dTwoRdmState, twoRdmState, h2)`

**Core routine** called during each time derivative evaluation:

| Argument        | Type                  | Description |
|-----------------|-----------------------|-------------|
| `dTwoRdmState`  | `complex(R64)(:)` INOUT | Packed derivative vector; collision term **added** in-place |
| `twoRdmState`   | `complex(R64)(:)` IN   | Packed 2-RDM state vector |
| `h2`            | `complex(R64)(:,:,:,:)` IN | Two-body Hamiltonian tensor |

### `Reconstruction_SetHamiltonian(h1, h2)`

For closures that depend on the Hamiltonian (`Singlet/Triples/`): the
one-body matrix of all spin orbitals and the h2 of the next
`Reconstruction_AddCollisionOp`, with every factor h2 carries in the
Hamiltonian. `td2rdmStaticOrbs` calls it before every AddCollisionOp and
passes λ h2 to both during adiabatic switching. The TD-2RDM method with
moving orbitals does not call it yet.

## JSON Configuration

A family block takes exactly one scheme, for example:

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
      "singletNy": { "ccQ": true }
    }
  }
}
```

```json
{
  "reconstruction": {
    "singletTr0": {
      "ny": {}
    }
  }
}
```

### Configuration Options

| Path | Type | Default | Description |
|------|------|---------|-------------|
| `reconstruction.singlet.<scheme>.ccQ` | bool | false | Contraction consistency of D3 (valdemoro, singletNy) |
| `reconstruction.singlet.triples.moellerPlesset` | object | absent | Connected triples (`Singlet/Triples/A_Reconstruction_Singlet_Triples.md`) |
| `reconstruction.singlet.anchor` | object | absent | Anchored closure (`Singlet/Anchor/A_Reconstruction_Singlet_Anchor.md`) |

## Adding a New Reconstruction Scheme

A singlet closure: an M_/S_ pair under `src/Reconstruction/Singlet/MyScheme/`
whose fabricate binds only `Reconstruction_Singlet_FillRdm3aab` (the contract
in `M_Reconstruction_Singlet`), never AddCollisionOp or Setup, and a branch
in `Reconstruction_Singlet_Fabricate` that sets `scheme` and calls it; the
family then reads `<scheme>.ccQ`. CMake picks up new
files by its glob.

Any other scheme (Generic, SingletTr0, or a collision operator of its own):

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

4. **Register in parent factory** (e.g., `S_Reconstruction_SingletTr0.f90`):
   ```fortran
   else if (Json_GetExistence("reconstruction.singletTr0.myScheme")) then
     call Reconstruction_SingletTr0_MyScheme_Fabricate
   ```

## Dependencies

| Module | Purpose |
|--------|---------|
| `M_Utils_RdmCumulants` | Cumulant expansion (free-part, NY corrections) |
| `M_Utils_RdmCollisionOp` | Compute Tr₃[Ŵ, D³] |
| `M_Utils_RdmTrace` | Extract 1-RDM from 2-RDM |
| `M_Utils_RdmEnforceTraceM2ab` | Project 2-RDM to trace-free subspace |
| `M_Utils_RdmEnforceTraceM3aab` | Contraction consistency of rdm3aab (ccQ) |
| `M_Utils_RdmSingletRelations` | aa block from the ab block of a singlet |
| `M_TwoRdm_*` | Pack/unpack state vectors |
| `M_Method_Mb` | Access particle numbers |
| `M_Method_Mb_OrbBased` | Access orbital counts |

## Testing

- `test/simulationTests/He1d/T_He1d_06_TimePropagationTd2rdm*`: TD-2RDM
  propagation (generic and singlet, weak orbitals, split step).
- `test/simulationTests/Be3d/T_Be3d_04_Td2rdmSpectralAnchor` pins the NY
  closure at a correlated Be state (the anchor tests: `Singlet/Anchor/`).
- `test/componentTests/Reconstruction/Singlet/Triples/T_Reconstruction_Singlet_Triples`:
  the triples cumulant (symmetries, vanishing contractions, covariance, pinned norm).
- `test/componentTests/Reconstruction/SingletTr0/T_Reconstruction_SingletTr0_Ny`,
  `test/componentTests/Utils/T_Utils_RdmEnforceTrace` (idempotent and
  contraction-free projections).

## References

- **Valdemoro approximation**: C. Valdemoro, Phys. Rev. A 45, 4462 (1992)
- **Nakatsuji-Yasuda**: H. Nakatsuji & K. Yasuda, Phys. Rev. Lett. 76, 1039 (1996)
