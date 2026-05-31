# Utils Module — Developer Onboarding Guide

## Overview

The **Utils** module is a collection of 68 utility modules providing the numerical
and algebraic infrastructure for the CodyFortranRDM framework. These modules are
**stateless** helper functions—they do not follow the Fabricate/Setup pattern and
have no JSON configuration. They are imported directly by other modules as needed.

This guide focuses on the **RDM algebra utilities**, which provide tensor operations
essential for reduced density matrix (RDM) theory and time-dependent simulations.

---

## Architecture Principles

### Stateless Design

All Utils modules are **pure utility libraries** with no internal state:
- No module-level variables storing configuration
- No Fabricate/Setup lifecycle (unlike core physics modules)
- Thread-safe by design (suitable for OpenMP parallelization)
- Directly importable via `use M_Utils_*`

### Naming Conventions

```
M_Utils_{Category}{Operation}.f90
       │         │
       │         └─ What it does (Algebra, Trace, Hole, etc.)
       └─ Domain (Rdm, Blas, Lapack, Fftw3, etc.)
```

Subroutine naming: `{Category}_{Operation}_{Details}`
```fortran
RdmAlgebra_FillSingletPart   ! Category=RdmAlgebra, Op=Fill, What=SingletPart
RdmTrace_FillRdm1aFromRdm2ab ! Category=RdmTrace, Op=Fill, What=Rdm1aFromRdm2ab
```

### Memory Management Patterns

**Allocatable outputs** — routines allocate if needed:
```fortran
complex(R64), intent(out), allocatable :: result(:,:,:,:)
! ...
if (.not. allocated(result)) allocate(result(n,n,n,n))
```

**In-place operations** — modify input directly:
```fortran
complex(R64), intent(inout), contiguous :: A2(:,:,:,:)
```

---

## RDM Utility Modules

| Module                          | Purpose                                                   |
|---------------------------------|-----------------------------------------------------------|
| `M_Utils_RdmAlgebra`            | Tensor symmetrization, antisymmetrization, spin-block decomposition |
| `M_Utils_RdmTrace`              | Partial trace operations (extract 1-RDM from 2-RDM, etc.) |
| `M_Utils_RdmCumulants`          | Cumulant expansions and decompositions                    |
| `M_Utils_RdmHole`               | Q-RDM (hole-RDM) and G-RDM construction                   |
| `M_Utils_RdmDefective`          | Extract negative eigenspace (defective part) of matrices  |
| `M_Utils_RdmPacking`            | Pack/unpack tensors to/from 1D state vectors              |
| `M_Utils_RdmProjection`         | Project corrections onto energy-conserving subspaces      |
| `M_Utils_RdmDiagonalize`        | Diagonalize RDMs to obtain natural orbitals/occupations   |
| `M_Utils_RdmObservables`        | Compute expectation values from RDMs                      |
| `M_Utils_RdmNatOrbTransform`    | Transform RDMs to natural orbital basis                   |
| `M_Utils_RdmCollisionOp`        | Compute collision operator Tr₃[Ŵ, D³]                     |
| `M_Utils_RdmCholeskyExtract`    | Cholesky decomposition of 2-RDM for ML features           |
| `M_Utils_RdmSingletRelations`   | Singlet spin-symmetry relations between RDM blocks        |
| `M_Utils_RdmEnforceTraceM2ab`   | Enforce trace constraints on α-β 2-RDM                    |
| `M_Utils_RdmEnforceTraceM2sym`  | Enforce trace constraints on symmetric 2-RDM              |
| `M_Utils_RdmEnforceTraceM2anti` | Enforce trace constraints on antisymmetric 2-RDM          |
| `M_Utils_RdmEnforceTraceM3aab`  | Enforce trace constraints on 3-RDM (αα-β block)           |

---

## M_Utils_RdmAlgebra Deep Dive

### Purpose

Provides fundamental tensor algebra operations for 4-index (rank-4) tensors
representing two-body reduced density matrices. These operations are building
blocks used throughout the framework for:

- **Cumulant decomposition**: Separating connected from disconnected correlations
- **Spin adaptation**: Decomposing RDMs into singlet/triplet irreducible components
- **N-representability**: Enforcing fermionic exchange symmetry

### Tensor Indexing Convention

All routines follow the physicist's convention for RDM indices:

```
A2(i1, i2, j1, j2) ≡ ⟨Ψ| a†_{i1} a†_{i2} a_{j2} a_{j1} |Ψ⟩
                   ≡ ⟨i1, i2 | Â | j1, j2⟩
                      └─────┘   └─────┘
                      creation  annihilation
                      (bra)     (ket)
```

- Indices 1,2 refer to particles 1,2 respectively
- α-β blocks: `rdm2ab(iα, iβ, jα, jβ)` — α-orbitals in dimensions 1,3

### Public Subroutines

| Subroutine                       | Description                                           |
|----------------------------------|-------------------------------------------------------|
| `RdmAlgebra_Antisymmetrize(A2)`  | Project tensor onto antisymmetric subspace (in-place) |
| `RdmAlgebra_Symmetrize(A2)`      | Project tensor onto symmetric subspace (in-place)     |
| `RdmAlgebra_FillAntisymProduct`  | Compute A0 ∧ A1 (fermionic wedge product)             |
| `RdmAlgebra_FillSymProduct`      | Compute A0 ⊙ A1 (bosonic symmetric product)           |
| `RdmAlgebra_FillSingletPart`     | Extract S=0 (singlet) component of α-β 2-RDM          |
| `RdmAlgebra_FillTripletPart`     | Extract S=1 (triplet) component of α-β 2-RDM          |

### Mathematical Background

#### Antisymmetrization (Fermionic)

For fermions, the 2-RDM must be antisymmetric under particle exchange:

```
D²(i1,i2,j1,j2) = -D²(i2,i1,j1,j2) = -D²(i1,i2,j2,j1) = +D²(i2,i1,j2,j1)
```

The projection operator P_A (idempotent: P_A² = P_A) is:

```
P_A[A] = (1/4) * [A(i1,i2,j1,j2) - A(i2,i1,j1,j2) - A(i1,i2,j2,j1) + A(i2,i1,j2,j1)]
```

#### Wedge Product (Exterior Product)

The antisymmetric product of two 1-RDMs gives the "free-particle" 2-RDM:

```
D²_free = D¹ ∧ D¹

D²_free(i1,i2,j1,j2) = (1/4) * [D¹(i1,j1)*D¹(i2,j2) - D¹(i2,j1)*D¹(i1,j2)
                               - D¹(i1,j2)*D¹(i2,j1) + D¹(i2,j2)*D¹(i1,j1)]
```

This is the mean-field (Hartree-Fock) approximation to the 2-RDM.
The 2-body cumulant Δ² = D² - D²_free measures electron correlation.

#### Singlet/Triplet Decomposition

For spin-1/2 fermions, the α-β 2-RDM decomposes into spin irreducible representations
(Clebsch-Gordan: ½ ⊗ ½ = 0 ⊕ 1):

```
D²_αβ = D²_singlet + D²_triplet
```

where:
- **Singlet (S=0)**: Symmetric under combined spatial+spin exchange
- **Triplet (S=1)**: Antisymmetric under combined exchange (vanishes on diagonal: i1=i2 or j1=j2)

This decomposition is crucial because:
1. Each block has independent N-representability constraints (D-condition, Q-condition)
2. Purification algorithms correct singlet/triplet blocks separately
3. For true singlet states, D²_triplet = 0 identically (useful diagnostic)

---

## Usage Examples

### Basic Symmetrization

```fortran
use M_Utils_RdmAlgebra
use M_Utils_Types

complex(R64), allocatable :: rdm2(:,:,:,:)

allocate(rdm2(nOrbs, nOrbs, nOrbs, nOrbs))
! ... fill rdm2 ...

! Project onto antisymmetric subspace (fermionic constraint)
call RdmAlgebra_Antisymmetrize(rdm2)
```

### Cumulant Decomposition

```fortran
use M_Utils_RdmAlgebra
use M_Utils_RdmCumulants

complex(R64), allocatable :: rdm1(:,:), rdm2(:,:,:,:)
complex(R64), allocatable :: rdm2_free(:,:,:,:), cumulant2(:,:,:,:)

! Build free-particle 2-RDM (disconnected part)
call RdmAlgebra_FillAntisymProduct(rdm2_free, rdm1, rdm1)

! Extract connected 2-body correlation (cumulant)
call RdmCumulants_FillCumul2(cumulant2, rdm1, rdm2)
! cumulant2 = rdm2 - rdm2_free
```

### Spin-Block Analysis

```fortran
use M_Utils_RdmAlgebra

complex(R64), allocatable :: rdm2ab(:,:,:,:)     ! α-β block
complex(R64), allocatable :: singlet(:,:,:,:)    ! S=0 part
complex(R64), allocatable :: triplet(:,:,:,:)    ! S=1 part

! Decompose into spin irreps
call RdmAlgebra_FillSingletPart(singlet, rdm2ab)
call RdmAlgebra_FillTripletPart(triplet, rdm2ab)

! Verify decomposition: rdm2ab = singlet + triplet
! Check singlet character: norm(triplet) ≈ 0 for S=0 states
```

---

## Related Modules

### M_Utils_RdmHole

Constructs the Q-RDM (hole-RDM) and G-RDM from the particle 2-RDM:

```
Q²(i,j,k,l) = δ_ik*δ_jl - δ_il*δ_jk - δ_jl*D¹_ki - δ_ik*D¹_lj + ... + D²(k,l,i,j)

G²(i,j,k,l) = δ_kl*D¹_ij + D²(i,l,j,k)
```

Key subroutines:
- `RdmHole_FillQrdm2ab` — Q-RDM for α-β block
- `RdmHole_FillQrdm2` — Q-RDM for full tensor
- `RdmHole_FillGrdm2` — G-RDM (particle-hole RDM)
- `RdmHole_FillQrdm1a` — 1-body hole density (q = 1 - γ)

### M_Utils_RdmCumulants

Implements cumulant expansions for RDM decomposition:

```
D² = D¹ ∧ D¹ + Δ²        (cumulant = Δ²)
D³ = D¹ ∧ D¹ ∧ D¹ + (Δ² ∧ D¹ + perms) + Δ³
```

Key subroutines:
- `RdmCumulants_FillCumul2` — Extract 2-body cumulant
- `RdmCumulants_FillCumul2ab` — Extract α-β 2-body cumulant
- `RdmCumulants_FillCumul3` — Extract 3-body cumulant
- `RdmCumulants_FillCumulFreePartM2` — Free-particle 2-RDM
- `RdmCumulants_FillCumulFreePartM3` — Free-particle 3-RDM
- `RdmCumulants_AddNyCumulantM3aab` — Nakatsuji-Yasuda 3-RDM correction

---

## Data Flow in TD-2RDM Propagation

The RDM algebra utilities are used throughout the TD-2RDM propagation pipeline:

```
┌─────────────────────────────────────────────────────────────────────┐
│                     Method_TimeDerivative                           │
└─────────────────────────────────────────────────────────────────────┘
                                │
                                ▼
┌─────────────────────────────────────────────────────────────────────┐
│ 1. Unpack twoRdmState → rdm2ab                                      │
│    (M_Utils_RdmPacking or TwoRdm module)                            │
└─────────────────────────────────────────────────────────────────────┘
                                │
                                ▼
┌─────────────────────────────────────────────────────────────────────┐
│ 2. Extract rdm1a from rdm2ab via partial trace                      │
│    (M_Utils_RdmTrace::RdmTrace_FillRdm1aFromRdm2ab)                 │
└─────────────────────────────────────────────────────────────────────┘
                                │
                                ▼
┌─────────────────────────────────────────────────────────────────────┐
│ 3. Reconstruct rdm3aab from rdm1a, rdm2ab                           │
│    (M_Utils_RdmCumulants::FillCumulFreePartM3aab + corrections)     │
└─────────────────────────────────────────────────────────────────────┘
                                │
                                ▼
┌─────────────────────────────────────────────────────────────────────┐
│ 4. Compute collision operator: collOp = Tr₃[Ŵ, D³]                  │
│    (M_Utils_RdmCollisionOp::RdmCollisionOp_FillCollisionOp2ab)      │
└─────────────────────────────────────────────────────────────────────┘
                                │
                                ▼
┌─────────────────────────────────────────────────────────────────────┐
│ 5. Purify: decompose into singlet/triplet, remove negative evals    │
│    (M_Utils_RdmAlgebra::FillSingletPart, FillTripletPart)           │
│    (M_Utils_RdmDefective::RdmDefective_FillDefectivePart)           │
│    (M_Utils_RdmHole::RdmHole_FillQrdm2ab for Q-condition)           │
└─────────────────────────────────────────────────────────────────────┘
```

---

## Performance Considerations

1. **O(nOrbs⁴) complexity**: Most RDM routines iterate over all tensor elements.
   Memory scales as O(nOrbs⁴) which limits practical orbital counts to ~50-100.

2. **Temporary allocations**: In-place operations (`Antisymmetrize`, `Symmetrize`)
   allocate a temporary copy. For repeated calls in tight loops, consider
   workspace-based variants if available.

3. **Spin symmetry exploitation**: For singlet states, only the singlet block
   carries information—use `FillSingletPart` to avoid processing zero elements.
   The triplet extraction already exploits diagonal zeros (i1≠i2, j1≠j2).

4. **Loop ordering**: Fortran is column-major. Current loop orderings (j2,j1,i2,i1
   outer-to-inner) are memory-friendly for the result tensor but access input
   in non-contiguous patterns. Profile before optimizing.

---

## Testing

The RDM algebra utilities are tested indirectly through higher-level tests:

```
test/T_He1d_06_TimePropagationTd2rdm.f90  — TD-2RDM propagation
test/T_Rdm_*.f90                          — RDM-specific unit tests (if present)
```

Run tests with:
```bash
cd build && ctest -R T_Rdm --output-on-failure
```

---

## Common Pitfalls

1. **Index convention mismatch**: Ensure calling code uses the same physicist's
   convention (creation before annihilation). Chemist's convention swaps indices.

2. **Missing antisymmetrization**: If constructing 2-RDMs from other sources,
   always call `RdmAlgebra_Antisymmetrize` before using in fermionic calculations.

3. **Singlet assumption**: Code paths assuming `rdm2ab ≡ rdm2singlet` will fail
   silently for non-singlet states. Always decompose and check if spin state
   is not guaranteed.

4. **Complex vs Real**: All routines use `complex(R64)` even for ground states.
   For real-valued RDMs, the imaginary part should be numerically zero.

---

## References

- **Cumulant theory**: W. Kutzelnigg & D. Mukherjee, J. Chem. Phys. 110, 2800 (1999)
- **Spin-adapted RDMs**: A.J. Coleman, Rev. Mod. Phys. 35, 668 (1963)
- **N-representability**: D.A. Mazziotti, Chem. Rev. 112, 244 (2012)
- **TD-2RDM theory**: K. Giesbertz & R. van Leeuwen, J. Chem. Phys. 139, 104109 (2013)
