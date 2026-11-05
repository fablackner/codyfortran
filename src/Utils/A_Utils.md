# Utils Module — Developer Onboarding Guide

## Overview

The **Utils** module is a collection of 72 utility modules providing the numerical
and algebraic infrastructure for the CodyFortran framework. These modules are
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
| `M_Utils_RdmDiagonalize`        | Diagonalize RDMs to obtain natural orbitals/occupations; `RdmDiagonalize_Rdm1(..., blockOf_)` diagonalizes per block (e.g. `Orbs_blockOfOrb`) |
| `M_Utils_RdmObservables`        | Compute expectation values from RDMs                      |
| `M_Utils_RdmNatOrbTransform`    | Transform RDMs to natural orbital basis                   |
| `M_Utils_RdmCollisionOp`        | Compute collision operator Tr₃[Ŵ, D³]                     |
| `M_Utils_RdmCholeskyExtract`    | Cholesky decomposition of 2-RDM for ML features           |
| `M_Utils_RdmSingletRelations`   | Singlet spin-symmetry relations between RDM blocks        |
| `M_Utils_RdmEnforceTraceM2ab`   | Enforce trace constraints on α-β 2-RDM; coefficients cached per nDim and recomputed when nDim changes |
| `M_Utils_RdmEnforceTraceM2sym`  | Enforce trace constraints on symmetric 2-RDM; coefficients cached per nDim and recomputed when nDim changes |
| `M_Utils_RdmEnforceTraceM2anti` | Enforce trace constraints on antisymmetric 2-RDM; coefficients cached per nDim and recomputed when nDim changes |
| `M_Utils_RdmEnforceTraceM3aab`  | Enforce trace constraints on 3-RDM (αα-β block); coefficients cached per nDim and recomputed when nDim changes |

---

## Observable Printers

| Module | Public routines | Notes |
|--------|-----------------|-------|
| `M_Utils_PrinterObservable` | `PrinterObservable_DumpNorm`, `_DumpOrthonormality`, `_DumpOrbitalEnergyMatrix`, `_DumpEnergy` | Method-aware: uses `M_Method`/`M_Method_Mb_OrbBased` (Apply*Op, `Method_Mb_OrbBased_FillInverseRdm1` for the orbital energy matrix), so not a stateless utility; one inner product per orbital (norms) and per pair (overlap matrix) |
| `M_Utils_PrinterObservableLinear` | `PrinterObservableLinear_DumpDipole`, `_DumpPotentialFullExpansion` | Linear grids |
| `M_Utils_PrinterObservableLattice` | `PrinterObservableLattice_DumpDensityCenter`, `_DumpDensityOnSite` | Lattice grids |
| `M_Utils_PrinterObservableYlm` | `PrinterObservableYlm_DumpDipole` (optional `dipole_` out), `PrinterObservableYlm_DumpDipoleAcceleration` (optional `acceleration_` out) | Ylm grids; z = r√(4π/3)Y₁₀ applied as spatial product; acceleration form −Z⟨cosθ/r²⟩ |
| `M_Utils_PrinterObservableFem` | `PrinterObservableFem_DumpDipole(rdm1, orbs, filename, toScreenQ, time, dipole_)` (line "time ⟨x⟩ ⟨y⟩ ⟨z⟩", optional `dipole_(3)` out), `PrinterObservableFem_DumpNorm(..., norm_)` (line "time N", N = Tr(ρ S)), `PrinterObservableFem_DumpDipoleVelocity(..., velocity_)` (line "time v_x v_y v_z", the velocity form v_a = i Tr(ρ_s [T, r_a]) of the dipole, d⟨r_a⟩/dt wherever the equations of motion obey Ehrenfest's theorem for r_a, e.g. TDHF; optional `velocity_(3)` out) | grid.fem; Method-aware like `M_Utils_PrinterObservable`: ρ_s = `Method_Mb_OrbBased_FillSpatialRdm1` of rdm1 (rdm1 itself, or the sum of the two spin blocks for restricted storage), body type of stored orbital i = `Method_Mb_OrbBased_bodyTypeOfOrb(i)`. Dipole and N: one pass over the nodes with the density ρ = Σ_ij ρ_s(i,j) φ_j* φ_i and the nodal weights (no orthonormality assumed), node coordinates from `Grid_Fem_FillNodeCoordinate`. Velocity: Re(i Σ_ij ρ_s(i,j)(⟨Tφ_j\|r_aφ_i⟩ − ⟨r_aφ_j\|Tφ_i⟩)) with `SysKinetic_MultiplyWithKineticOp` (one T application per orbital, 2 nOS² inner products per axis); the potential, the length gauge and the interaction commute with r_a, so after a kick e^{iκr_a} it gives the TRK value Nκ/m on the grid (P_KickTdhfFem and FemKick.inc take their initial slope from it); early return without output request |

---

## Binary Data I/O: `M_Utils_DataStorage`

`SaveData(data, filename, elementBits, startPos_)` and `LoadData(...)` move
the raw bytes of an array of any type and rank (assumed rank, `c_loc`) to and
from an unformatted stream file (orbitals `orbXX_YY.in`, `coeffs.in`,
`rdm2BtXX_YY.in`). Both use one stream transfer of all bytes, as 8- or 4-byte
words when the size allows (libgfortran moves 1-byte items one by one):
0.05 s per 16 MB (1e6-point) orbital instead of 2.7 s (save) and 1.8 s (load)
with the former one I/O statement per byte; the files are byte-identical.
Byte counts and positions are 64-bit. Behaviour unchanged: a file shorter than requested
stops LoadData with "more data requested than stored in ...", a file read up
to its end prints "reached end of ...", a longer one is read in part without
a note.

---

## Block-Diagonal Eigenproblems: `M_Utils_LapackLib`

`LapackLib_DiagonalizeHermitianBlocks(evals, evecs, a, blockOf)` diagonalizes
a Hermitian matrix that is block diagonal in the labels `blockOf` block by
block (ZHEEV): the eigenvectors are exactly zero outside their block also
where eigenvalues of different blocks are degenerate (a diagonalization of
the whole matrix mixes such blocks at round-off). Eigenvalues ascending over
all blocks; with a single block the result equals
`LapackLib_DiagonalizeGeneric` bitwise. Used for orbital-space matrices with
the blocks `Orbs_blockOfOrb` (body type and declared m): the inverse 1-RDM of
the correlation operator and the Loewdin steps of the Be/Ne programs.
Test: `test/componentTests/Utils/T_Utils_LapackLib.f90`.

---

## GPU Linear Algebra: `M_Utils_CublasLib`

Minimal double-precision cuBLAS layer (`src/Utils/CublasLib/`): device memory
as `type(c_ptr)` plus element offsets, `CublasLib_AvailableQ`,
`CublasLib_Allocate`/`Deallocate`, `CublasLib_CopyToDevice`/`CopyToHost` and
`CublasLib_Dgemm` (column-major, 'N'/'T'). The CUDA runtime and cuBLAS are
called through `bind(C)`, so no CUDA compiler is needed. CMake compiles
exactly one implementation submodule: `S_Utils_CublasLib_Cuda` with
`-DCODY_USE_CUDA=ON` (links `CUDA::cublas`, `CUDA::cudart`), otherwise
`S_Utils_CublasLib_Stub` (no device; every call except `AvailableQ` stops).
Used by the exact field-free propagation (`method.mb.orbBased.exactFieldFreeGpuQ`).

---

## Spherical Harmonics: `M_Utils_SphericalHarmonics`

Gaunt coefficients (direct via GSL Wigner-3j, or from the table built by
`SphericalHarmonics_EnsureGauntTable` and read by
`SphericalHarmonics_GauntTabulated`). The module also owns the compact (l,m)
channel layout, 0 ≤ l ≤ lmax, |m| ≤ min(l, mmax), l-major with m ascending:
`SphericalHarmonics_ChannelOffset`, `SphericalHarmonics_ChannelIndex` and
`SphericalHarmonics_NChannels`. Both the Gaunt table and the Grid_Ylm storage
(`Grid_Ylm_NChannels`) use it.

---

## Solid Harmonics: `M_Utils_SolidHarmonics`

Real regular solid harmonics S_lm(x, y, z) in Racah normalization,
S_lm = √(4π/(2l+1)) r^l Y_lm (real Y_lm; m > 0 cos-type, m < 0 sin-type),
from the Helgaker/Jørgensen/Olsen recursion (S₀₀ = 1, S₁₀ = z, S₁₁ = x,
S₁,₋₁ = y, S₂₀ = (3z² − r²)/2, S₂₂ = √3/2 (x² − y²), S₂,₋₂ = √3 xy, …), so
that 1/|r₁ − r₂| = Σ_lm S_lm(r₁) S_lm(r₂)/|r₂|^(2l+1) for |r₁| < |r₂|.
Stateless, thread-safe; tables are indexed by
`SolidHarmonics_LmIndex(l, m)` = l² + l + m + 1.

| Routine | Purpose |
|---------|---------|
| `SolidHarmonics_FillCoefficients(coeffs, lmax)` | Polynomial coefficients: S_lm = Σ_(a+b+c=l) coeffs(a, b, c, lm) x^a y^b z^c (multipole moments and boundary values of `Grid_Fem_SolvePoisson`) |
| `SolidHarmonics_FillValues(values, r, lmax)` | Values of all S_lm, l ≤ lmax, at the point r (pure; no division by r, regular at r = 0) — used by the Fem STO initializer |

Test: `test/componentTests/Utils/T_Utils_SolidHarmonics` (closed forms for
l ≤ 3, values vs. polynomial coefficients for l ≤ 8 to 7e-16, addition
theorem to 3e-16, regular at the origin).

---

## HHG Spectrum: `M_Utils_HhgSpectrum`

Stateless spectrum of a dipole time series, written at the end of the laser
programs in `app/` (`hhg.dat`).

| Routine | Purpose |
|---------|---------|
| `HhgSpectrum_Write(outFile, time, dipole, omega0, acceleration_, maxOrder_, dOrder_)` | Hann-windowed Fourier transforms over the whole equidistant series: acceleration form \|FT[a w]\|² (if `acceleration_` is given) and length form ω⁴ \|FT[(z − ⟨z⟩) w]\|², versus harmonic order ω/ω0 from `dOrder_` (0.1) to `maxOrder_` (300), stopping at the Nyquist frequency |
| `HhgSpectrum_WriteFromTable(outFile, inFile, omega0, colTime, colDipole, colAcceleration_, maxOrder_, dOrder_)` | Same, with the series read from columns of a whitespace table (`#` lines skipped) |

---

## Laser Pulses: `M_Utils_Pulse`

Stateless, pure pulse shapes shared by the laser couplings.

| Routine | Purpose |
|---------|---------|
| `Pulse_Sin2Field(time, amplitude, omega, nCycles, cep, tStart)` | E(t) = E0 sin²(π t_rel/T) cos(ω (t_rel − T/2) + cep), t_rel = t − tStart, T = 2π nCycles/ω; zero outside [tStart, tStart + T] and for ω ≤ 0 (carrier centered on the envelope maximum) |

Used by the length-gauge couplings `sysPotential.ylm.coulombLaser`
(`SysPotential_Ylm_CoulombLaser_FieldAmplitude`; test T_Ne3d_07) and
`sysGauge.fem.length` (`SysGauge_Fem_Length_FieldAmplitude`; tests
T_SysGauge_Fem_Length and the Fem length-gauge tests of Method), with the
same floating-point operations as their former private copies (bit-identical
outputs). The velocity gauge `sysGauge.ylm.velocity` keeps its own shapes:
its envelope sits on A(t) = −(E0/ω) sin²(π t_rel/T) sin(ω (t_rel − T/2) + cep),
or, with envelope `field`, on E(t) with the uncentered carrier
cos(ω t_rel + cep).

---

## Linear Systems: `M_Utils_Gmres`

Restarted, flexible, right-preconditioned GMRES(m) for a complex system
A y = b in an inner product of the caller, matrix-free: the operator, the
preconditioner and the inner product are callbacks, the context is owned by
the caller (the pattern of `M_Utils_ConvolutionFftw`).

| Routine / type | Purpose |
|---------|---------|
| `T_Gmres_Ctx` | Controls `krylovDim` (m), `maxRestarts`, `tolerance`, and the workspace: Arnoldi vectors (n × m+1), their preconditioned images (n × m) and the residual — allocated by the first solve for size(b) and kept (a solve of another size reallocates) |
| `Gmres_CreateCtx(ctx, krylovDim, maxRestarts, tolerance)` | Fresh context with these controls (releases an earlier workspace; error stop for m < 1, restarts < 0, tolerance ≤ 0) |
| `Gmres_DestroyCtx(ctx)` | Releases the workspace |
| `Gmres_Solve(y, nIterations, relativeResidual, convergedQ, ctx, b, ApplyPrecOnVec, ApplyMatOnPrecVec, ApplyMatOnVec, InnerProduct)` | y ≈ A⁻¹ b; returns the Arnoldi iterations, the last true residual ‖b − A y‖/‖b‖ and whether it meets the tolerance (neither prints nor stops) |

Callbacks (abstract interfaces, public): `I_Gmres_ApplyPrecOnVec(z, v)`,
z = M⁻¹ v; `I_Gmres_ApplyMatOnPrecVec(w, v, z)`, w = A z for z = M⁻¹ v (v is
given so that A = M + N can form v + N z without a second solve);
`I_Gmres_ApplyMatOnVec(w, y)`, A applied explicitly for the true residuals;
`I_Gmres_InnerProduct(fConjg, f)`, pure, the characteristics of
`I_Grid_InnerProduct` (pass `Grid_InnerProduct` directly). Pass module
procedures: internal procedures as actual arguments need trampolines (an
executable stack). Callbacks may read data the caller prepared for the solve.

Algorithm: initial guess y0 = M⁻¹ b; every convergence decision, the first
one included, on the true residual (y0 and the end of every cycle); modified
Gram–Schmidt Arnoldi of A M⁻¹ with the preconditioned vectors kept (flexible
form, no extra solve for the update); complex Givens rotations; a cycle ends
after m iterations, once the estimated residual meets the tolerance, or when
the new Arnoldi vector vanishes (happy breakdown); at most 1 + maxRestarts
cycles; b = 0 gives y = 0. A solve touches only its context and arguments
(solves with different contexts may run concurrently); the vector passes are
OpenMP parallel outside parallel regions.

User: the Cayley step of the exact field-free part of Method on grid.fem
(`S_Method_Mb_OrbBased_FieldFreeCayley`: P⁻¹ = (1 + i Δt T/2)⁻¹ by
`SysKinetic_SolveShifted`, `Grid_InnerProduct`). Test:
`test/componentTests/Utils/T_Utils_Gmres` (dense non-Hermitian n = 40 system
in a diagonal metric against ZGESV with one cycle and with m = 2 restarts,
not converged at maxRestarts = 0, exact preconditioner, happy breakdown,
b = 0, workspace life cycle).

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
