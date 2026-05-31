# CoeffsInit Module — Agent Context

## Purpose

**CoeffsInit** initializes the Configuration Interaction (CI) coefficient vector `c` that defines the many-body quantum state:

```
|Ψ⟩ = Σᵢ cᵢ |Φᵢ⟩
```

where `|Φᵢ⟩` are Fock-space configurations enumerated by the `ConfigList` module. This module controls **what state the simulation starts from**.

---

## Architecture

```
M_CoeffsInit.f90          ← Public interface (procedure pointers)
S_CoeffsInit.f90          ← Fabrication dispatcher

Unary/                    ← c = [1, 0, 0, ...] (reference state)
├── M_CoeffsInit_Unary.f90
└── S_CoeffsInit_Unary.f90

Load/                     ← Read from coeffs.in
├── M_CoeffsInit_Load.f90
└── S_CoeffsInit_Load.f90

Excited/                  ← Apply â†/â operators to |Φ₀⟩
├── M_CoeffsInit_Excited.f90
└── S_CoeffsInit_Excited.f90
```

**Pattern:** Interface module (`M_`) exports procedure pointers; submodule (`S_`) provides implementation. Fabrication reads JSON and binds pointers at runtime.

---

## Available Backends

| Backend | JSON Key | Description | Typical Use |
|---------|----------|-------------|-------------|
| **Unary** | `coeffsInit.unary` | Sets `c(1)=1`, all others zero | Ground-state search, mean-field |
| **Load** | `coeffsInit.load` | Reads from `coeffs.in` | Restarts, checkpoints |
| **Excited** | `coeffsInit.excited` | Applies â†/â to build excitations | Response, dynamics |

---

## JSON Configuration Examples

### Unary (most common)
```json
"coeffsInit": {
  "unary": { }
}
```

### Load from file
```json
"coeffsInit": {
  "load": { }
}
```
Expects binary file `coeffs.in` in working directory.

### Excited state
```json
"coeffsInit": {
  "excited": {
    "creates":   [3],
    "destroys":  [1],
    "bodyType1": 1,
    "bodyType2": 2
  }
}
```
Applies `â†₃ â₁` to both body types (e.g., spin-up and spin-down).

---

## Call Sequence

```fortran
! 1. Fabrication (binds procedure pointers)
call CoeffsInit_Fabricate

! 2. Setup (optional preprocessing — usually no-op)
call CoeffsInit_Setup

! 3. Initialize (fills coefficient vector)
call CoeffsInit_Initialize(coeffs)
```

All calls after fabrication go through the bound procedure pointers in `M_CoeffsInit`.

---

## Key Interfaces

```fortran
! From M_CoeffsInit:
procedure(I_CoeffsInit_Setup),      pointer :: CoeffsInit_Setup
procedure(I_CoeffsInit_Initialize), pointer :: CoeffsInit_Initialize

! Initialize signature:
subroutine I_CoeffsInit_Initialize(coeffs)
  complex(R64), intent(out), contiguous :: coeffs(:)
end subroutine
```

---

## Dependencies

| Module | Role |
|--------|------|
| `M_Utils_Types` | Type kinds (R64, I32) |
| `M_Utils_Json` | JSON configuration parsing |
| `M_Utils_DataStorage` | Binary I/O for Load backend |
| `M_Coeffs` | `Coeffs_ApplyExcitation` for Excited backend |
| `M_Utils_NoOpProcedures` | Default no-op for Setup pointer |

---

## Adding a New Backend

1. Create directory `src/CoeffsInit/MyInit/`
2. Add interface module `M_CoeffsInit_MyInit.f90`:
   - Declare `CoeffsInit_MyInit_Fabricate` in interface block
3. Add submodule `S_CoeffsInit_MyInit.f90`:
   - Implement fabrication (bind procedure pointers)
   - Implement local `Initialize(coeffs)` procedure
4. Update `S_CoeffsInit.f90`:
   - Add `use M_CoeffsInit_MyInit`
   - Add branch: `if (Json_GetExistence("coeffsInit.myInit")) then ...`
5. Update CMakeLists.txt to include new files

---

## Common Pitfalls

- **Dimension mismatch:** `coeffs(:)` must match the CI space size from `ConfigList`. The Load backend assumes the file has exactly the right size.
- **Body types in Excited:** For single-species systems, set `bodyType1 = bodyType2 = 1`.
- **Branching order:** In `S_CoeffsInit.f90`, Load takes precedence over Unary over Excited. Only one backend is activated.

---

## Related Modules

- **Coeffs** — CI state container and Hamiltonian application
- **ConfigList** — Fock-space configuration enumeration
- **OrbsInit** — Analogous initializer for orbital coefficients
- **TwoRdmInit** — Analogous initializer for 2-RDM state
