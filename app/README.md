# Applications

Production runs built on the CodyFortran library. Every leaf folder holds one
program `P_X.f90` and its input `P_X.json` (and `P_X.exe` after the build).
Programs are named `P_<Stage><Method><Variant>`: stage `Ground` (ground or
initial state), `Laser` (laser pulse) or `Kick` (delta kick); method `Scf`,
`Mcscf`, `Mctdhf`, `Td2rdm`, `Tdhf`, `Ks` (Kohn–Sham DFT) or `Tddft`; variant the gauge (`Velocity`,
`Length`) or the grid (`Ylm`, `Fem`) where several exist.

```
atoms/
  be/  groundstate/hf                          P_GroundScfYlm         Hartree-Fock (1s, 2s)
       groundstate/mcscf                       P_GroundMcscf          MCSCF, 5 orbitals (ground solver)
       groundstate/mctdhf                      P_GroundMctdhf         MCTDHF, 5 orbitals (imaginary time)
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
             groundstate/dft                   P_GroundKsFem          Kohn-Sham DFT (sysXc, Libxc)
             dynamics/tdhf                     P_KickTdhfFem          absorption spectrum
             dynamics/tddft                    P_KickTddftFem         absorption spectrum, real-time TDDFT
             analysis/absorption_spectrum.py
```

## Running

```bash
source setCodyVars.sh -c gnu -t release
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

Be in a z-polarized 800 nm pulse, velocity gauge. MCTDHF and TD-2RDM use 5
restricted orbitals 1s 2s 2p₋₁ 2p₀ 2p₊₁ (full CI in the total-M = 0 sector
for MCTDHF), Hartree-Fock the orbitals 1s 2s (SCF ground solver, E =
−14.5730236). Ylm grid with
lmax = 32 (mmax 1, Hartree-Fock mmax 0) on FEDVR 50 × 20 points (rmax 200),
cos^¼ absorber from 160 a.u.

Order: `groundstate/hf`, then the correlated ground state with either
`groundstate/mcscf` (ground solver `groundSolver.mcscf.gradient`, 10
iterations, 2 s: E = −14.6168560449) or `groundstate/mctdhf` (imaginary-time
MCTDHF, about 15 min: E = −14.6168560112, 3.4e-8 higher from its dtau = 0.005
step bias), then `groundstate/td2rdm` (TD-2RDM initial state from the MCTDHF
ground state: field-free TD-2RDM propagation averaged over 40 a.u.), then any
`dynamics/<method>/<pulse>`. Both correlated runs start from the
Hartree-Fock 1s 2s of `../hf` plus hydrogen-like 2p (`orbsInit.load.nOrbsIn`)
and the CI eigenvector of these orbitals (`program.ciStartQ` for imaginary
time: from a unary CI vector it ends in a local minimum at E = −14.5956,
where the 2p shell correlates the 1s² core), and both write the same files
(the dynamics load `groundstate/mctdhf`; point `orbsInit.load.dir` and
`coeffsInit.load.dir` at `groundstate/mcscf` to use the MCSCF state).

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

The ground states use DIIS on the orbitals (`mixTarget: "orbitals"`,
`mixing.diis.startThreshold: 0`; `src/Mixing/A_Mixing.md`, "What to mix").
H2: R = 1.4 bohr, SCF in 8 steps (4 s), 50 a.u. of TDHF in about 15 min.
H2O: SCF in 12 steps (22 with linear Hartree-potential mixing), about
1.5 min on two threads; the SCF mesh is the coarse one of the kick (the accurate mesh with 1.1 M nodes,
1.7 mHa below the HF limit, is `app/fem_molecules/h2o_scf.json` in the git
history). For the isotropic H2O spectrum run three kicks
(`program.kickDirection` 1, 2, 3, each in its own copy of the folder) and
pass all three `dipole.dat` files to the script. Restarts use
`program.saveStep` and `program.startTime` with `orbsInit.load` reading the
checkpoint of the run directory (no `dir`); keep `meanFieldReference.dat`.

## H2 and H2O: Kohn–Sham DFT and TDDFT kick spectra on grid.fem

The same steps with a density functional (`sysXc.fem.libxc`, `src/SysXc/A_SysXc.md`):

1. `groundstate/dft` (`P_GroundKsFem`, method `tdks`) writes the Kohn–Sham
   orbitals and `femMesh.dat`; the functional is the `sysXc` block of its
   input (LDA `["lda_x", "lda_c_vwn"]` in the inputs; GGA `["gga_x_pbe",
   "gga_c_pbe"]`, a global hybrid such as `["hyb_gga_xc_pbeh"]` with `aceQ`).
   H2 LDA: E = −1.1375189, HOMO −10.27 eV (7 SCF steps, 6 s on one
   thread); PBE: −1.1667267, HOMO −10.38 eV. H2O LDA (linear mixing
   oscillates): E = −75.9635786, HOMO −7.40 eV, μ = 0.728 a.u. (12 steps,
   about 1.5 min on two threads). That is 48.6 mHa below
   the accurate mesh (nLocals 6, hMin 0.05: −75.9150234), the same offset as
   HF (49.1 mHa, the Galerkin O nucleus); the valence ε agree to 1e-3 Ha.
2. `dynamics/tddft` (`P_KickTddftFem`) reads both from `../../groundstate/dft`
   (same mesh, checked; use the same `sysXc` block), builds the frozen
   Kohn–Sham reference U = v_KS[ρ0] − a K0 of the ground state, kicks every
   orbital with e^{iκ r_a} and propagates with real-time TDDFT. It writes
   `dipole.dat` and `observables.dat` in the format of `dynamics/tdhf`.
3. The spectrum as above.

The time steps of the inputs (H2 0.05, H2O 0.1, as for TDHF) hold for LDA.
With a GGA the split step is stable only if every propagated h0 mode has
E Δt < π, because the kernel of a GGA grows like k², and the stiff nodes at
the nuclei reach E ~ 10³ a.u. H2 PBE is stable at Δt 0.005 but diverges at
0.02 and 0.05 (src/Method/A_Method.md, "Time step with a GGA"). A PBE ground
state is fine at any setting; for PBE dynamics, reduce the time step by
about 10× (more for H2O, whose O core has hMin 0.1).

Background on the methods: `src/Method/A_Method.md` (exact field-free part,
Cayley step, frozen mean-field reference) and `src/Propagator/A_Propagator.md`.
Earlier benchmarks and timings are in the git history of `app/be/README.md`,
`app/ne/hhg/README.md` and `app/fem_molecules/README.md`.
