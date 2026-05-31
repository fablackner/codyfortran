# DiagonalizerList Module — Agent Onboarding Guide

## Purpose

The **DiagonalizerList** module provides a unified interface for eigenvalue problem solvers used in quantum many-body simulations. It abstracts away the differences between dense (LAPACK) and iterative (ARPACK) eigensolvers, allowing the rest of the codebase to request eigenpairs through a **matrix-free** interface.

**Core problem solved:** Given a linear operator Ĥ (typically a Hamiltonian), find eigenvalues λᵢ and eigenvectors |ψᵢ⟩ satisfying:

```
Ĥ |ψᵢ⟩ = λᵢ |ψᵢ⟩
```

The operator is accessed only through a callback `ApplyMatOnVec(y, x, t)` that computes y = Ĥ(t)·x.

---

## Architecture Overview

```
src/DiagonalizerList/
├── M_DiagonalizerList.f90        # Public interface (abstract type + procedure pointers)
├── S_DiagonalizerList.f90        # Factory logic (JSON dispatch)
├── AGENTS.md                     # This file
├── Lapack/
│   ├── M_DiagonalizerList_Lapack.f90   # LAPACK backend interface
│   └── S_DiagonalizerList_Lapack.f90   # LAPACK implementation
└── Arpack/
    ├── M_DiagonalizerList_Arpack.f90   # ARPACK backend interface
    └── S_DiagonalizerList_Arpack.f90   # ARPACK implementation
```

### Design Pattern: Strategy + Factory

| Component | Role |
|-----------|------|
| **M_DiagonalizerList** | Abstract interface `T_DiagonalizerList_E` with deferred methods |
| **S_DiagonalizerList** | Factory that reads JSON and instantiates backends by name |
| **T_DiagonalizerList_E_Lapack** | Dense solver: assembles explicit matrix, calls LAPACK |
| **T_DiagonalizerList_E_Arpack** | Iterative solver: Arnoldi/Lanczos via ARPACK |

---

## Data Flow

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                           Client Code                                        │
│  1. Define callback: ApplyMatOnVec(dState, state, time)                     │
│  2. Create input: T_DiagonalizerList_FabricateInput                         │
└─────────────────────────────────────────────────────────────────────────────┘
                                   │
                                   ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                    DiagonalizerList_Fabricate(input)                         │
│  1. Read JSON "diagonalizerList" block                                       │
│  2. For each child: allocate backend by name (lapack/arpack)                │
│  3. Bind callback, set dim, call backend%Fabricate()                        │
└─────────────────────────────────────────────────────────────────────────────┘
                                   │
                                   ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                       DiagonalizerList_Setup                                 │
│  Calls Setup() on each backend (allocate workspaces)                        │
└─────────────────────────────────────────────────────────────────────────────┘
                                   │
                                   ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│              DiagonalizerList(i)%e%Diagonalize(time, evecsQ)                │
│  Compute eigenpairs → store in %evals(:), %evecs(:,:), %nFound              │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## JSON Configuration

### Example: Multiple Solvers

```json
{
  "diagonalizerList": {
    "arpack_ground": {
      "nEvals": 5,
      "which": "SR",
      "nKry": 20,
      "checkConvergenceQ": true,
      "printLevel": 0
    },
    "lapack_full": {
      "nEvals": -1,
      "printLevel": 0
    }
  }
}
```

**Backend selection:** The JSON child name must contain `"lapack"` or `"arpack"` as a substring.

### LAPACK Parameters

| Parameter    | Type | Default | Description                       |
|--------------|------|---------|-----------------------------------|
| `nEvals`     | int  | 1       | Number of eigenvalues (-1 = all)  |
| `printLevel` | int  | 0       | Verbosity (0=silent, 1=progress)  |

### ARPACK Parameters

| Parameter           | Type   | Default      | Description                              |
|---------------------|--------|--------------|------------------------------------------|
| `nEvals`            | int    | 1            | Number of eigenvalues (-1 = all)         |
| `which`             | string | "SR"         | Spectrum region (see below)              |
| `bmat`              | string | "I"          | "I" = standard, "G" = generalized        |
| `nKry`              | int    | 2*nEvals+1   | Krylov subspace dimension                |
| `checkConvergenceQ` | bool   | true         | Verify residuals after solve             |
| `printLevel`        | int    | 0            | Verbosity level                          |

**Spectrum selection (`which`):**
- `"LM"` / `"SM"` : Largest / Smallest Magnitude
- `"LR"` / `"SR"` : Largest / Smallest Real part ← **use for ground states**
- `"LI"` / `"SI"` : Largest / Smallest Imaginary part

---

## Key Interfaces

### ApplyMatOnVec Callback

```fortran
subroutine ApplyMatOnVec(dState, state, time)
  complex(R64), intent(out), contiguous :: dState(:)  ! Result: y = H·x
  complex(R64), intent(in), contiguous  :: state(:)   ! Input: x
  real(R64), intent(in)                 :: time       ! For time-dependent H
end subroutine
```

This is the **only** interface the eigensolver uses to access the operator.

### Fabrication Input

```fortran
type(T_DiagonalizerList_FabricateInput) :: inputs(1)
inputs(1)%ApplyMatOnVec => MyHamiltonianCallback
inputs(1)%dim = hilbert_space_size
call DiagonalizerList_Fabricate(inputs)
```

### Accessing Results

```fortran
call DiagonalizerList(1)%e%Diagonalize(0.0_R64, .true.)

! Results:
eigenvalues  => DiagonalizerList(1)%e%evals(1:nFound)   ! real(R64)
eigenvectors => DiagonalizerList(1)%e%evecs(:, 1:nFound) ! complex(R64)
num_found    =  DiagonalizerList(1)%e%nFound
```

---

## Typical Usage in CodyFortranRDM

```fortran
! 1. Setup fabrication input (callback + dimension)
type(T_DiagonalizerList_FabricateInput) :: DiagonalizerListInput(1)
DiagonalizerListInput(1)%ApplyMatOnVec => Actions_ApplyHamiltonian
DiagonalizerListInput(1)%dim = Coeffs_nCoeffs

! 2. Fabricate (reads JSON, creates backends)
call DiagonalizerList_Fabricate(DiagonalizerListInput)

! 3. Setup (allocate workspaces)
call DiagonalizerList_Setup

! 4. Solve eigenvalue problem
call DiagonalizerList(1)%e%Diagonalize(0.0_R64, .true.)

! 5. Use results
if (DiagonalizerList(1)%e%nFound > 0) then
  ground_state_energy = DiagonalizerList(1)%e%evals(1)
  ground_state_vector = DiagonalizerList(1)%e%evecs(:, 1)
end if
```

---

## Backend Comparison

| Aspect | LAPACK | ARPACK |
|--------|--------|--------|
| **Matrix access** | Assembles explicit dim×dim matrix | Matrix-free (callback only) |
| **Complexity** | O(dim³) | O(nKry × dim × T_matvec) |
| **Memory** | O(dim²) | O(nKry × dim) |
| **Best for** | dim ≲ 1000, full spectrum | dim > 1000, partial spectrum |
| **Convergence** | Deterministic | Iterative (may fail) |

**Rule of thumb:** Use ARPACK for production (large systems), LAPACK for validation/testing.

---

## Adding a New Backend

1. **Create directory:** `src/DiagonalizerList/MyBackend/`

2. **Interface module** (`M_DiagonalizerList_MyBackend.f90`):
   ```fortran
   module M_DiagonalizerList_MyBackend
     use M_DiagonalizerList, only: T_DiagonalizerList_E

     type, extends(T_DiagonalizerList_E) :: T_DiagonalizerList_E_MyBackend
       ! backend-specific fields
     contains
       procedure :: Fabricate
       procedure :: Setup
       procedure :: Diagonalize
     end type

     interface
       module subroutine DiagonalizerList_MyBackend_Allocate(e, path)
         class(T_DiagonalizerList_E), allocatable, intent(out) :: e
         character(len=*), intent(in) :: path
       end subroutine
     end interface
   end module
   ```

3. **Implementation submodule** (`S_DiagonalizerList_MyBackend.f90`):
   - Implement `Allocate`, `Fabricate`, `Setup`, `Diagonalize`

4. **Register in factory** (`S_DiagonalizerList.f90`):
   ```fortran
   use M_DiagonalizerList_MyBackend, only: DiagonalizerList_MyBackend_Allocate
   ...
   else if (index(childName, "mybackend") .ne. 0) then
     call DiagonalizerList_MyBackend_Allocate(DiagonalizerList(i)%e, "diagonalizerList."//childName)
   ```

5. **Update CMakeLists.txt** to include new source files

---

## Dependencies

| Depends On | Purpose |
|------------|---------|
| `M_Utils_Types` | R64, I32 type kinds |
| `M_Utils_Json` | JSON configuration parsing |
| `M_Utils_Say` | Logging/diagnostics |
| `M_Utils_LapackLib` | LAPACK wrapper (ZHEEVR/DSYEVR) |
| `M_Utils_ArpackLib` | ARPACK wrapper (znaupd/zneupd) |
| `M_Utils_BlasLib` | BLAS norms for residual checks |
| `M_Utils_NoOpProcedures` | Default no-op for Setup pointer |

---

## Common Pitfalls

1. **Input array size mismatch:** `size(input)` must equal number of JSON children under `"diagonalizerList"`.

2. **Callback dimension mismatch:** `input(i)%dim` must match the vector size expected by `ApplyMatOnVec`.

3. **ARPACK convergence failure:** If `nFound < nEvals`, increase `nKry` or check operator conditioning.

4. **Forgetting Setup:** Always call `DiagonalizerList_Setup` before `Diagonalize`.

5. **Time parameter:** For time-independent Hamiltonians, pass any value (e.g., `0.0_R64`).

6. **Name convention:** JSON child names must contain `"lapack"` or `"arpack"` substring.

---

## Testing

Relevant test files in `test/`:
- `T_FermiHubbard_01_ExactDiagTdciLapack.f90` — LAPACK on Fermi-Hubbard
- `T_FermiHubbard_02_ExactDiagTdciArpack.f90` — ARPACK on Fermi-Hubbard
- `T_BoseHubbard_01_ExactDiagTdciLapack.f90` — LAPACK on Bose-Hubbard
- `T_BoseHubbard_02_ExactDiagTdciArpack.f90` — ARPACK on Bose-Hubbard

Run tests:
```bash
cd build && ctest -R ExactDiag
```

---

## Performance Tuning

### ARPACK
- **Increase `nKry`:** Better convergence, more memory. Rule: `nKry ≈ 2-3 × nEvals`.
- **Decrease `nKry`:** Faster iterations, may fail to converge.
- **`checkConvergenceQ=false`:** Skip residual verification (faster, less safe).

### LAPACK
- **`nEvals=-1`:** Compute all eigenvalues (O(dim³) cost).
- **`nEvals=k`:** Compute only lowest k (uses ZHEEVR index range).

---

## Related Modules

- **Method** — Dynamics methods that use DiagonalizerList for ground-state initialization
- **Propagator/EigenExpansion** — Time evolution via diagonalization
- **GroundSolver** — SCF iteration using eigensolvers
