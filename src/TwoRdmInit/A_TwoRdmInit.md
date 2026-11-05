# TwoRdmInit Module — Agent Onboarding Guide

## Purpose

The **TwoRdmInit** module provides initialization strategies for the two-body reduced density matrix (2-RDM) in quantum many-body simulations. The 2-RDM encodes all two-particle correlations:

```
D²(i,j;k,l) = ⟨Ψ| a†_i a†_j a_l a_k |Ψ⟩
```

This quantity is central to:
- **TD-2RDM dynamics**: Time propagation of the 2-RDM using BBGKY hierarchy
- **Energy evaluation**: E = Tr[h₁·D¹] + ½Tr[W·D²]
- **Correlation analysis**: Beyond mean-field two-particle correlations

---

## Architecture Overview

```
src/TwoRdmInit/
├── M_TwoRdmInit.f90              # Public interface (procedure pointers)
├── S_TwoRdmInit.f90              # Factory dispatcher + block assembly
├── AGENTS.md                     # This file
├── Excited/
│   ├── M_TwoRdmInit_Excited.f90  # Excited-state interface
│   └── S_TwoRdmInit_Excited.f90  # Build D² from 1-RDM products
└── Load/
    ├── M_TwoRdmInit_Load.f90     # File-load interface
    └── S_TwoRdmInit_Load.f90     # Read D² from binary files
```

### Design Pattern: Factory + Strategy

| Component | Role |
|-----------|------|
| **M_TwoRdmInit** | Public API: `TwoRdmInit_Fabricate`, `TwoRdmInit_Setup`, `TwoRdmInit_Initialize`, `TwoRdmInit_InitializeBlock` |
| **S_TwoRdmInit** | Factory logic: reads JSON, dispatches to implementation |
| **Excited** | Strategy: construct D² from 1-RDM products (Hartree-Fock-like) |
| **Load** | Strategy: read D² blocks from binary files |

---

## Data Flow

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                           JSON Configuration                                 │
│  { "twoRdmInit": { "excited": { ... } } }                                   │
│  { "twoRdmInit": { "load": { } } }                                          │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                      TwoRdmInit_Fabricate                                    │
│  1. Parse JSON "twoRdmInit" block                                           │
│  2. Dispatch to Excited or Load Fabricate                                   │
│  3. Bind: TwoRdmInit_Initialize → Initialize (block assembly)               │
│           TwoRdmInit_InitializeBlock → implementation-specific              │
└─────────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                      TwoRdmInit_Initialize(twoRdmState)                      │
│  1. Allocate full tensor rdm2(nO, nO, nO, nO)                               │
│  2. For each (bt1, bt2) ∈ {1,2}×{1,2}:                                      │
│     a. Allocate block-sized tensor                                          │
│     b. Call TwoRdmInit_InitializeBlock(block, bt1, bt2)                     │
│     c. Copy block into full tensor at correct slice                         │
│  3. Pack tensor → linearized state via TwoRdm_FillTwoRdmStateFromRdm2       │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## JSON Configuration

### Excited: Build from 1-RDM Products

```json
{
  "twoRdmInit": {
    "excited": {
      "creates": [3],
      "destroys": [1],
      "bodyType1": 1,
      "bodyType2": 2
    }
  }
}
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `creates` | int[] | [1] | Orbital indices to occupy |
| `destroys` | int[] | [1] | Orbital indices to vacate |
| `bodyType1` | int | 1 | Body type for creates/destroys |
| `bodyType2` | int | 2 | Second body type |

**Note**: The creates/destroys parameters are currently read but not fully wired up in the implementation. See "Known Issues" below.

### Load: Read from Files

```json
{
  "twoRdmInit": {
    "load": { }
  }
}
```

Optional `"dir"`: directory of the files (default: the working directory). Requires the files:
- `rdm2Bt01_01.in` — α-α block
- `rdm2Bt01_02.in` — α-β block
- `rdm2Bt02_01.in` — β-α block
- `rdm2Bt02_02.in` — β-β block

---

## Mathematical Details

### 2-RDM Block Structure

For a two-species fermionic system (e.g., spin-up α and spin-down β):

```
        ┌─────────────────┬─────────────────┐
        │     bt2 = 1     │     bt2 = 2     │
        │    (α orbitals) │    (β orbitals) │
┌───────┼─────────────────┼─────────────────┤
│bt1 = 1│   D²_αα(i,j;k,l)│   D²_αβ(i,j;k,l)│
│  (α)  │  antisymmetric  │   no symmetry   │
├───────┼─────────────────┼─────────────────┤
│bt1 = 2│   D²_βα(i,j;k,l)│   D²_ββ(i,j;k,l)│
│  (β)  │   no symmetry   │  antisymmetric  │
└───────┴─────────────────┴─────────────────┘
```

### Excited Implementation: Mean-Field Factorization

For same-spin blocks (bt1 == bt2), fermionic antisymmetry is enforced:
```
D²(i₁,i₂;j₁,j₂) = D¹(i₁,j₁)·D¹(i₂,j₂) - D¹(i₂,j₁)·D¹(i₁,j₂)
```

For opposite-spin blocks (bt1 ≠ bt2), simple product:
```
D²(i₁,i₂;j₁,j₂) = D¹(i₁,j₁)·D¹(i₂,j₂)
```

The 1-RDM is constructed as a diagonal matrix with occupation numbers:
- n_i = 1 for i ≤ nBodies (occupied orbitals)
- n_i = 0 otherwise (virtual orbitals)

---

## Typical Usage

```fortran
! Initialization sequence
call TwoRdm_Fabricate           ! Set up 2-RDM module first
call TwoRdmInit_Fabricate       ! Reads JSON, binds implementation

call TwoRdm_Setup
call TwoRdmInit_Setup           ! Currently no-op

! Initialize the 2-RDM state
complex(R64), allocatable :: twoRdmState(:)
allocate(twoRdmState(TwoRdm_nTwoRdmState))
call TwoRdmInit_Initialize(twoRdmState)

! Now twoRdmState is ready for propagation
call Propagator_Propagate(state, t0, t1)
```

---

## Dependencies

| Depends On | Purpose |
|------------|---------|
| `M_Utils_Types` | R64, I32 type kinds |
| `M_Utils_Json` | JSON configuration parsing |
| `M_Utils_Say` | Logging/diagnostics |
| `M_Utils_NoOpProcedures` | Default no-op for Setup pointer |
| `M_Utils_DataStorage` | Binary file I/O (Load implementation) |
| `M_TwoRdm` | 2-RDM state packing (`FillTwoRdmStateFromRdm2`) |
| `M_Method_Mb` | Particle counts (`nBodies`) |
| `M_Method_Mb_OrbBased` | Orbital indices (`nOrbs`, `nOrbsStart`, `nOrbsEnd`) |

---

## Adding a New Initialization Strategy

1. **Create directory**: `src/TwoRdmInit/MyStrategy/`

2. **Interface module** (`M_TwoRdmInit_MyStrategy.f90`):
   ```fortran
   module M_TwoRdmInit_MyStrategy
     use M_Utils_Types
     implicit none

     interface
       module subroutine TwoRdmInit_MyStrategy_Fabricate
       end subroutine
     end interface
   end module
   ```

3. **Implementation submodule** (`S_TwoRdmInit_MyStrategy.f90`):
   ```fortran
   submodule(M_TwoRdmInit_MyStrategy) S_TwoRdmInit_MyStrategy
     implicit none
   contains

     module subroutine TwoRdmInit_MyStrategy_Fabricate
       use M_TwoRdmInit
       TwoRdmInit_InitializeBlock => MyInitializeBlock
     end subroutine

     subroutine MyInitializeBlock(rdm2Bt, bt1, bt2)
       complex(R64), intent(out), contiguous :: rdm2Bt(:,:,:,:)
       integer(I32), intent(in) :: bt1, bt2
       ! ... your implementation
     end subroutine

   end submodule
   ```

4. **Register in factory** (`S_TwoRdmInit.f90`):
   ```fortran
   use M_TwoRdmInit_MyStrategy

   ! In TwoRdmInit_Fabricate:
   else if (Json_GetExistence("twoRdmInit.mystrategy")) then
     call TwoRdmInit_MyStrategy_Fabricate
   ```

5. **Update CMakeLists.txt** to include new source files

---

## Known Issues and TODOs

### 1. Unused Configuration Parameters (Excited)

The `creates` and `destroys` arrays are read from JSON but the local variables `initCreatesBt1`, `initDestroysBt1`, etc. are hardcoded to `[1]`:

```fortran
! In S_TwoRdmInit_Excited.f90
creates = Json_Get("twoRdmInit.excited.creates", [1])    ! Read but unused!
destroys = Json_Get("twoRdmInit.excited.destroys", [1])  ! Read but unused!
...
allocate(initDestroysBt1, source=[1])  ! Hardcoded!
```

**Suggested fix**: Wire up the module-level `creates`/`destroys` to the block initialization logic, respecting `bodyType1`/`bodyType2`.

### 2. Memory Efficiency

The `Initialize` routine allocates a full O(nOrbs⁴) tensor before packing. For large orbital spaces, consider direct packing strategies that work block-by-block without materializing the full tensor.

### 3. Missing File Existence Check (Load)

The Load implementation will crash if files don't exist. Consider adding existence checks with informative error messages:

```fortran
inquire(file=trim(filename), exist=file_exists)
if (.not. file_exists) then
  error stop "TwoRdmInit_Load: File not found: "//trim(filename)
end if
```

---

## Testing

The module is exercised by TD-2RDM propagation tests:
- `test/T_He1d_06_TimePropagationTd2rdm.f90` — Uses Load strategy

Run relevant tests:
```bash
cd build && ctest -R Td2rdm
```

---

## Related Modules

| Module | Relationship |
|--------|--------------|
| **TwoRdm** | Defines 2-RDM state vector layout and operations |
| **Method** | TD-2RDM dynamics methods that consume initialized 2-RDMs |
| **Reconstruction** | BBGKY closure approximations for 3-RDM |
| **Purification** | N-representability enforcement on 2-RDMs |
| **CoeffsInit** | Analogous initializer for CI coefficient vectors |
| **OrbsInit** | Analogous initializer for orbital coefficients |
