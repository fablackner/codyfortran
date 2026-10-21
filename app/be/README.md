# Be: high-harmonic generation with TD-2RDM vs MCTDHF

Recreation of the beryllium results (Fig. 4, and the time-frequency analysis of
Fig. 6) of

> F. Lackner, I. Březinová, T. Sato, K. L. Ishikawa, and J. Burgdörfer,
> *High-harmonic spectra from time-dependent two-particle reduced density
> matrix theory*, Phys. Rev. A **95**, 033414 (2017),
> [arXiv:1611.00888](https://arxiv.org/abs/1611.00888).

A Be atom in a z-polarized 800 nm few-cycle pulse is propagated with
MCTDHF and with the TD-2RDM method (time-dependent orbitals, propagation of
the (↑↓) block of the 2-RDM, contraction-consistent Nakatsuji–Yasuda
reconstruction of the 3-RDM, D/Q purification for the strong pulse). The
induced dipole moment and the HHG spectrum (from the dipole acceleration via
the Larmor formula) of both methods are compared.

## Workflow

```bash
source setBuildVars.sh -c gnu -t release
cody build app                         # builds the programs below
cd app/be
python3 make_inputs.py                 # writes inputs/*.json (one parameter set)
./run.sh inputs runs                   # all steps, or e.g. ./run.sh inputs runs relaxA relaxB
python3 analysis/plot_fig4.py runs --out fig4.png --spectrogram fig6.png
```

| step | program | what it does |
|---|---|---|
| `relaxA` | `relax/P_ImagTimeMctdhf` | imaginary-time MCTDHF ground state at small `lmax` |
| `relaxB` | `relax/P_ImagTimeMctdhf` | embed into the propagation grid, relax again |
| `average` | `td2rdmInit/P_Td2rdmTimeAverage` | TD-2RDM initial state: field-free TD-2RDM propagation and time average over 40 a.u. (appendix C of the paper) |
| `mctdhf_<pulse>` | `laser/P_LaserMctdhf` | MCTDHF in the pulse |
| `td2rdm_<pulse>` | `laser/P_LaserTd2rdm` | TD-2RDM in the pulse |

Pulses (`weak4c`, `strong2c`): F(t) = F₀ cos(ωt) sin²(ωt/2N_c), ω = 0.057 a.u.,
velocity gauge (`sysGauge.ylm.velocity` with `"envelope": "field"`, i.e. the
sin² envelope on the field as in the paper).

| pulse | I [W/cm²] | F₀ [a.u.] | N_c | purification |
|---|---|---|---|---|
| weak4c | 0.5×10¹⁴ | 0.038 | 4 | none |
| strong2c | 4×10¹⁴ | 0.107 | 2 | 40 iterations per step |

`make_inputs.py --test DIR` writes a tiny configuration (small box, few
orbitals, half-cycle pulses) to check the pipeline in minutes. `--norbs N`
keeps the first N orbitals, `--nweak-orbs NW` switches on the TD-2RDM
active/weak natural-orbital split (NW weak spatial natural orbitals of the
spin-summed 1-RDM, see `src/Method/A_Method.md`; off by default) for all
TD-2RDM inputs, i.e. the time average and the laser runs.

**Time step.** Each output step of 0.2 a.u. has 80 split steps (dt = 0.0025)
for TD-2RDM and MCTDHF alike. With 6 or more orbitals dt = 0.01 is unstable for
TD-2RDM (the inverse 1-RDM makes the orbital equations stiff, ~1/n_min) and
lets the MCTDHF energy drift (1.8e-4 in 10 a.u. vs. 9e-7 at dt = 0.0025).

**Exact field-free propagation (default; `--no-exact-field-free` to switch
off).** The one-body split step is stiff (centrifugal term l(l+1)/2r² ≈ 4e6
a.u. at the first FEDVR point for l = 32), so SIL needs ~250 Krylov vectors
per half step. With `method.mb.orbBased.exactFieldFreeQ` the field-free h0 is
propagated exactly per l (diagonalized once) and SIL only sees the gauge
coupling A·p. A plain Strang splitting of h0 and A·p adds a large error near
the nucleus ([h0, A·p] = iA ∂zV); the fourth-order Blanes–Moan composition
(`propagator.splitStep.exactPartScheme = "blanesMoan6"`, 7 exact and 6 SIL
substeps per linear step; now the default is its swapped variant
`"blanesMoan6Swap"` with 6 exact and 7 SIL substeps, see "Performance" below)
removes it. MCTDHF, strong2c pulse, 0.2 a.u. at the
peak vector potential, maximal error against the Richardson-extrapolated
reference, one thread:

| linear step | dt | ⟨z⟩ | −Z⟨z/r³⟩ (acceleration) | energy | wall time |
|---|---|---|---|---|---|
| SIL (previous default) | 0.0025 | 1.0e-6 | 6.7e-4 | 8.9e-6 | 3014 s |
| SIL | 0.00125 | 2.5e-7 | 1.7e-4 | 2.2e-6 | 2030 s |
| exact h0, Strang | 0.0025 | 7.5e-6 | 2.3e-2 | 2.0e-4 | 366 s |
| exact h0, Strang | 0.000625 | 4.7e-7 | 1.3e-3 | 1.2e-5 | 1648 s* |
| exact h0, Yoshida-4 | 0.0025 | 1.1e-6 | 1.0e-2 | 1.1e-4 | 579 s |
| **exact h0, Blanes–Moan** | 0.0025 | 1.0e-6 | 6.7e-4 | 9.4e-6 | 678 s |
| **exact h0, Blanes–Moan** | 0.00125 | 2.5e-7 | 1.7e-4 | 2.2e-6 | 1247 s |

(*shared the two cores with other runs part of the time.) With Blanes–Moan the
errors coincide with SIL at the same dt: what remains is the outer
(one-body/interaction) Strang error common to both. TD-2RDM agrees with the
SIL run to 7e-10 in ⟨z⟩ (627 s vs. 3185 s). Per step early in the pulse (weak
field, one thread): 17.4 s with SIL, 4.9 s with Blanes–Moan, 3.1 s with the
GPU. The imaginary-time relaxation (`relax/P_ImagTimeMctdhf`, no field) uses
exp(−h0 dτ) directly: one check of stage B takes 157 s instead of 278 s
(including ~45 s one-time diagonalization) with the same energy to 3e-14.

**GPU (`--gpu`).** The exact h0 transformations (two dense GEMMs per l, 14 per
time step) can run on an NVIDIA GPU with cuBLAS; the eigenvectors stay on the
device (227 MiB). Build once with the cuBLAS backend (no CUDA compiler needed,
the toolkit libraries suffice):

```bash
cmake -B build -DCODY_USE_CUDA=ON      # then cody build app as usual
python3 make_inputs.py --gpu
```

On a GTX 1050 Ti (FP64 at 1/32 rate) this gives 1.6× per step over the CPU
path; results agree with the CPU to round-off.

**Performance (MCTDHF, 5 orbitals, strong2c, at the peak vector potential).**
Wall time per split step (dt = 0.0025) on a Pentium G4560 (2 cores, no AVX),
OpenBLAS, `OMP_NUM_THREADS` as given:

| | 1 thread | 2 threads |
|---|---|---|
| before (blanesMoan6, SIL tolerance 1e-12) | 4.54 s | 3.17 s |
| now (CPU) | 3.71 s | 2.45 s |
| now (`--gpu`, GTX 1050 Ti) | 2.33 s | 1.74 s |

What changed, measured against a dt/4 reference over 0.1 a.u. (maximal error
of the dipole acceleration 6.1e-4 in all cases):

- `exactPartScheme = "blanesMoan6Swap"`: the Blanes–Moan composition with the
  roles exchanged, 6 instead of 7 exact substeps (the dense transformations
  dominate); the observables change by 2.5e-5 (acceleration).
- SIL tolerance 1e-10 in real time (1e-12 in imaginary time): fewer Krylov
  vectors; the observables change by 1e-10.
- `ApplyCorrelationOpYlmNodes` with blocked ZGEMMs: results unchanged to 1e-13,
  the mean-field part 4× cheaper.
- MRRR (DSYEVR) for the exact field-free blocks: setup 34 s → 13 s.

Two threads pay off (1.5×) when no other job runs; the GPU path only offloads
the dense transformations, so the CPU threads still matter there. Tried and
dropped: SHTns (the m = 0 angular transforms are ~1.5% of the run), the
reversed outer splitting B(½Δt) A(Δt) B(½Δt) (1.2× faster per step but twice
the acceleration error, i.e. slower at equal accuracy), padded or re-laid-out
DGEMMs (≤ 4%; the 999 × 999 × 10 products already run at 60% of peak).

**Spectral cut (optional, `python3 make_inputs.py --spectral-cut 1000`).** The
highest h0 eigenvalues of the FEDVR grid (up to 3.6e5 a.u.) are unphysical.
`method.mb.orbBased.exactFieldFreeCut` drops the eigenvectors above the cut
(1000 a.u. keeps 895–900 of 950 per l) from the exact part and projects the
orbital derivatives onto the rest (src/Method/A_Method.md). Without that
stiffness the plain Strang composition (2 instead of 12 exact transformations)
is as accurate as Blanes–Moan: at the peak 2.12 s vs 4.78 s per step (one
thread), acceleration error 3.9e-3 vs 5.3e-3 for Blanes–Moan with the same cut,
while Strang without cut is off by 1.8e-2. The cut changes the model slightly.
Full strong2c pulse on a small grid (rmax 40, ground states relaxed with the
cut, dt = 0.0025) against the uncut Blanes–Moan reference: the harmonics 21–60
change by 0.02–0.04 decades (2e-3 for the uncut schemes) at 1.7× less wall
time; cuts of 300 and 1000 give the same change. Larger steps do not pay off:
at dt = 0.005 every scheme, with or without cut, is off by ~0.5 decades above
harmonic 20. MCTDHF only; TD-2RDM inputs stay uncut, so use it for MCTDHF
studies, not for the MCTDHF/TD-2RDM comparison.

**Exponential integrators (`propagator.exponential`, not the default).** Same
small grid and pulse, 4 threads, against the uncut reference (Strang with cut
1000 at dt 0.0025: harmonics 21–40 / 41–60 off by 0.019 / 0.036 decades,
2266 s):

| scheme | cut | dt | h21–40 / h41–60 | wall |
|---|---|---|---|---|
| `adamsLawson`, k = 4 | 1000 | 0.0025 | 0.36 / 0.42 | 1298 s |
| `adams`, k = 4 | 300 | 0.0025 | 0.19 / 0.23 | 1450 s |
| `adams`, k = 4 | 100 | 0.005 | 0.70 / 0.52 | 723 s |
| `krogstad4` | 1000 | 0.0025 | 0.28 / 0.26 | 2418 s |

N itself is not stiff (Jacobian spectral radius ~26 a.u. at the ground
state). What limits these schemes is that N couples orbital components
oscillating with h0 frequencies up to the cut (2.5 rad per step at cut 1000,
dt 0.0025), which the interpolation of N in time cannot follow: `adams`
diverges at cut 1000 for k ≥ 3 and runs stably at cut 300. Splitting avoids
the interpolation (N is integrated with h0 frozen).

**Purification.** With `conserveEnergyQ` the D/Q correction is orthogonal to
h2, so `laser/P_LaserTd2rdm` rebuilds h2 from the current orbitals before
every purification (previously h2 was refreshed only at output steps). With
`"purifyAverageQ": True` in `make_inputs.py`, `td2rdmInit/P_Td2rdmTimeAverage`
purifies the time-averaged initial 2-RDM once with the paper's settings (the
average is not N-representable; the minimal D and Q eigenvalues, printed
before and after, drop by one to two orders of magnitude to about the
threshold; without `conserveEnergyQ` the energy shifts by ~1e-3).

## Numerical parameters

| | paper | here |
|---|---|---|
| radial grid | FEDVR, 50 elements × 20 points, 200 a.u. | same |
| absorber | cos^¼ mask from 160 a.u. | same |
| angular basis | spherical harmonics, L_max = 47 | L_max = 32, **m = 0 only** |
| orbitals | not stated | 6 spatial orbitals (m = 0), restricted |
| reconstruction | contraction-consistent Nakatsuji–Yasuda | same |
| propagation | 2nd-order split step | same (one-body part: exact h0 with Blanes–Moan composition of the gauge SIL; RK4 for the rest) |
| initial state | MCTDHF ground state by imaginary time; TD-2RDM: field-free time average (T = 40 a.u.) | same |

Restricting the orbitals to m = 0 is exact for the dynamics of m = 0 orbitals
in a z-polarized field (m is conserved) and reduces the Ylm grid from
(L+1)² to L+1 channels; it restricts, however, the orbital space of the
correlated ground state (no 2p±1 orbitals). MCTDHF and TD-2RDM use the
identical orbital space, so the comparison of the two methods is unaffected.

## Library changes made for this project

Building this comparison exposed and required several changes in `src/`:

- **Coeffs Generic index convention** (bug fix): the CI backend applied the
  transposed one-body/two-body matrices and returned transposed RDMs relative
  to the convention of `Method_Mb_OrbBased_FillH1/H2`, TwoRdm and the orbital
  equations of motion. Invisible for real Hamiltonians and orbitals, wrong for
  complex ones (velocity gauge). Verified against exact propagation (one and
  two electrons, complex basis) to 1e-15; MCTDHF and TD-2RDM now coincide for
  N = 2 (where the 3-RDM vanishes).
- **MCTDHX split step** (bug fix): the non-linear part also applied h¹ to the
  CI coefficients although the linear part already propagates the one-body
  part exactly; the energy drifted by ~0.07 a.u. per a.u. of time in a He test.
- **TD-2RDM split step**: `TimeDerivativeOrbsLin` /
  `TimeDerivativeCoeffsPlusOrbsNonLin` for the TD-2RDM method.
- **SIL imaginary time** (`imagTimeQ`), used for the ground state on the stiff
  FEDVR grid.
- **SysGauge** `"envelope": "field"` (the paper's pulse) and
  `SysGauge_Ylm_Velocity_FieldAmplitude`.
- **Dipole acceleration** observable `PrinterObservableYlm_DumpDipoleAcceleration`.
- **Ylm spatial product**: exact Gauss–Legendre fast path for m = 0 operands
  (3–5× faster propagation).
- **TD-2RDM active/weak natural-orbital split**
  (`method.mb.orbBased.td2rdm.nWeakOrbs`, off by default): stabilizes TD-2RDM
  with weakly occupied natural orbitals; see `src/Method/A_Method.md`.
- **CPU performance** (results unchanged to round-off, 1.8× per step on one
  thread, 2.1× on four): the m = 0 correlation operator
  (now `ApplyCorrelationOpYlmNodes`) was never selected because the Ylm Coulomb
  interaction did not declare itself body-type independent; central
  potentials multiply as a radial scaling instead of the angular quadrature;
  SIL reorthogonalizes with two `ZGEMV` sweeps; the FEDVR element operator
  avoids `matmul` temporaries.
- **Gaunt table** restricted to the m of the grid and operators: setup takes
  milliseconds and ~1 MB instead of ~3.5 min and 2.4 GB per program start.
- **Exact field-free one-body propagation** with a fourth-order composition
  (`method.mb.orbBased.exactFieldFreeQ`, `propagator.splitStep.exactPartScheme`)
  and an optional cuBLAS backend (`exactFieldFreeGpuQ`, CMake `CODY_USE_CUDA`),
  default for this project, see "Time step" above.
- **Contraction-consistent trace constraints** (bug fix): the
  `M_Utils_RdmEnforceTrace*` coefficients were computed for the first
  dimension only and silently reused for other dimensions.

## Results

(see below)
