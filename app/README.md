# Applications

Production runs built on the CodyFortran library. Every leaf folder holds one
program `P_X.f90` and its input `P_X.json` (and `P_X.exe` after the build).
Programs are named `P_<Stage><Method><Variant>`: stage `Ground` (ground or
initial state), `Laser` (laser pulse) or `Kick` (delta kick); method `Scf`,
`Mcscf`, `Mctdhf`, `Td2rdm` or `Tdhf`; variant the gauge (`Velocity`,
`Length`) or the grid (`Ylm`, `Fem`) where several exist.

```
atoms/
  be/  groundstate/mctdhf                      P_GroundMctdhf         MCTDHF, 5 orbitals
       groundstate/hf                          P_GroundScfYlm         Hartree-Fock (1s, 2s)
       groundstate/td2rdm                      P_GroundTd2rdm         TD-2RDM initial state
       dynamics/{mctdhf,hf}/{weak4c,strong2c}  P_LaserMctdhfVelocity
       dynamics/td2rdm/{weak4c,strong2c}       P_LaserTd2rdmVelocity
       analysis/plot_hhg.py
  ne/  groundstate/mctdhf                      P_GroundMctdhf         MCTDHF, 9 orbitals
       groundstate/hf                          P_GroundScfYlm         Hartree-Fock
       groundstate/mcscf                       P_GroundMcscfYlm       MCSCF, 9 orbitals
       dynamics/mctdhf                         P_LaserMctdhfVelocity  HHG, velocity gauge
       dynamics/mctdhfLengthGauge              P_LaserMctdhfLength    from mcscf
molecules/
  h2/, h2o/  groundstate/scf                   P_GroundScfFem
             dynamics/tdhf                     P_KickTdhfFem          absorption spectrum
             analysis/absorption_spectrum.py
```

## Running

```bash
source setBuildVars.sh -c gnu -t release
cody build app
export CODY_PRETTY_PRINT=OFF OMP_NUM_THREADS=2
cd app/atoms/ne/groundstate/mctdhf && ./P_GroundMctdhf.exe
cd ../../dynamics/mctdhf && ./P_LaserMctdhfVelocity.exe
```

A program reads `P_X.json` next to its executable (or the JSON file given as
the only argument) and writes its output to the working directory, so run it
in its folder. A ground-state run writes the ground state (`orb*.in`,
`coeffs.in`, for MCTDHF also `rdm2Bt*.in`). A dynamics run reads it from the
ground-state folder (`orbsInit.load.dir`, `coeffsInit.load.dir`,
`twoRdmInit.load.dir` in its JSON) and, for the atoms, writes the HHG spectrum
`hhg.dat` at the end (harmonic order, acceleration form |a(ω)|², length form
ω⁴|z(ω)|²; Hann window over the whole pulse, up to `program.hhgMaxOrder`,
default 300). Outputs (`*.dat`, `*.in`) are git-ignored.

The atomic ground states are computed at lmax 1 (Be Hartree-Fock: 0) and
embedded into the lmax of the propagation grid when the dynamics run loads
them (`orbsInit.load.lmaxIn`; ground state and dynamics need the same
`grid.ylm.mmax`). This is exact: the orbitals of a ¹S state with complete p
shells (2p₋₁, 2p₀, 2p₊₁ declared in `method.mb.orbBased.mOfOrbs`) keep a
pure l.

## Be: HHG with MCTDHF, Hartree-Fock and TD-2RDM

Be in a z-polarized 800 nm pulse, velocity gauge. MCTDHF and TD-2RDM use 6
restricted orbitals 1s 2s 2p₋₁ 2p₀ 2p₊₁ (full CI in the total-M = 0 sector
for MCTDHF), Hartree-Fock the orbitals 1s 2s (SCF ground solver, E =
−14.5730236). Ylm grid with
lmax = 32 (mmax 1, Hartree-Fock mmax 0) on FEDVR 50 × 20 points (rmax 200),
cos^¼ absorber from 160 a.u.

Order: `groundstate/mctdhf` and `groundstate/hf`, then `groundstate/td2rdm`
(TD-2RDM initial state from the MCTDHF ground state: field-free TD-2RDM
propagation averaged over 40 a.u.), then any `dynamics/<method>/<pulse>`.

Pulses: `weak4c` (0.5×10¹⁴ W/cm², 4 cycles) and `strong2c` (4×10¹⁴ W/cm²,
2 cycles, D/Q purification for TD-2RDM). From `app/atoms/be`,
`python3 analysis/plot_hhg.py dynamics --out hhg.png --spectrogram spectrogram.png`
plots dipoles and spectra of all methods and the time-frequency spectrum of
the strong pulse; next to each `hhg.dat` it writes `hhg_larmor.dat` with its
own unwindowed Larmor spectrum.

## Ne: ground states and strong-field MCTDHF

- `groundstate/mctdhf` → `dynamics/mctdhf`: 9 restricted orbitals 1s 2s 2p
  3s 3p, full CI in the total-M = 0 sector, 800 nm two-cycle pulse at
  5×10¹⁴ W/cm², lmax = 40, output every 0.2 a.u. (88 000 steps of about
  2 s). With `program.outputStep` 0.05 and `absorber.keepOrthonormalQ: true`
  the acceleration spectrum stays clean up to harmonic ~1100; at 0.2 a.u. it
  shows an aliasing floor above H95.
- `groundstate/hf`, `groundstate/mcscf`: HF and 9-orbital MCSCF ground
  states on a small grid (lmax 1, rmax 8). `dynamics/mctdhfLengthGauge`
  propagates the MCSCF state in a length-gauge pulse on that grid; its
  `hhg.dat` holds the length form only.

## H2 and H2O: RHF and TDHF kick spectra on grid.fem

1. `groundstate/scf` writes the orbitals and `femMesh.dat`.
2. `dynamics/tdhf` reads both from `../../groundstate/scf` (same mesh,
   checked), kicks every orbital with e^{iκ r_a} and propagates with TDHF. It
   writes `dipole.dat` and `observables.dat`.
3. `python3 ../../analysis/absorption_spectrum.py dipole.dat --pade` gives the
   polarizability, the absorption cross section and the TRK sum.

H2: R = 1.4 bohr, SCF in 5 s, 50 a.u. of TDHF in about 15 min. H2O: the SCF
mesh is the coarse one of the kick (the accurate mesh with 1.1 M nodes,
1.7 mHa below the HF limit, is `app/fem_molecules/h2o_scf.json` in the git
history). For the isotropic H2O spectrum run three kicks
(`program.kickDirection` 1, 2, 3, each in its own copy of the folder) and
pass all three `dipole.dat` files to the script. Restarts use
`program.saveStep` and `program.startTime` with `orbsInit.load` reading the
checkpoint of the run directory (no `dir`); keep `meanFieldReference.dat`.

Background on the methods: `src/Method/A_Method.md` (exact field-free part,
Cayley step, frozen mean-field reference) and `src/Propagator/A_Propagator.md`.
Earlier benchmarks and timings are in the git history of `app/be/README.md`,
`app/ne/hhg/README.md` and `app/fem_molecules/README.md`.
