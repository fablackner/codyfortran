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
orbitals, half-cycle pulses) to check the pipeline in minutes.

## Numerical parameters

| | paper | here |
|---|---|---|
| radial grid | FEDVR, 50 elements × 20 points, 200 a.u. | same |
| absorber | cos^¼ mask from 160 a.u. | same |
| angular basis | spherical harmonics, L_max = 47 | L_max = 32, **m = 0 only** |
| orbitals | not stated | 6 spatial orbitals (m = 0), restricted |
| reconstruction | contraction-consistent Nakatsuji–Yasuda | same |
| propagation | 2nd-order split step | same (SIL for the one-body part, RK4 for the rest) |
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

## Results

(see below)
