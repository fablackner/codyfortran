# Applications

Production runs built on the CodyFortran library. Every leaf folder holds one
program `P_X.f90` and its input `P_X.json` (and `P_X.exe` after the build);
programs with a `make_inputs.py` take an input it wrote to `runs/` as the
argument.
Programs are named `P_<Stage><Method><Variant>`: stage `Ground` (ground or
initial state), `Laser` (laser pulse) or `Kick` (delta kick); method `Scf`,
`Mcscf`, `Mctdhf`, `Td2rdm`, `Tdhf`, `Ks` (Kohn–Sham DFT) or `Tddft`; variant the gauge (`Velocity`,
`Length`) or the grid (`Ylm`, `Fem`, `Prolate`) where several exist; analysis
programs `P_<Quantity><Grid>` (`P_IntegralsFem`, `P_StepCostFem`).

```
atoms/
  be/  groundstate/scf                         P_GroundScfYlm         Hartree-Fock (1s, 2s)
       groundstate/mcscf                       P_GroundMcscf          MCSCF, 5 orbitals (ground solver)
       groundstate/mctdhf                      P_GroundMctdhf         MCTDHF, 5 orbitals (imaginary time)
       groundstate/td2rdm                      P_GroundTd2rdm         TD-2RDM initial state
       dynamics/{mctdhf,hf}/{weak4c,strong2c}  P_LaserMctdhfVelocity
       dynamics/td2rdm/{weak4c,strong2c}       P_LaserTd2rdmVelocity
       analysis/plot_hhg.py
  ne/  groundstate/mctdhf                      P_GroundMctdhf         MCTDHF, 9 orbitals
       groundstate/scf                         P_GroundScfYlm         Hartree-Fock
       groundstate/mcscf                       P_GroundMcscf          MCSCF, 9 orbitals (ground solver)
       groundstate/td2rdm                      P_GroundTd2rdm         TD-2RDM by adiabatic switching, no CI
       dynamics/mctdhf                         P_LaserMctdhfVelocity  HHG, velocity gauge
       dynamics/mctdhfLengthGauge              P_LaserMctdhfLength    from mcscf
molecules/
  h2/, h2o/  groundstate/scf                   P_GroundScfFem
             groundstate/dft                   P_GroundKsFem          Kohn-Sham DFT (sysXc, Libxc)
             dynamics/tdhf                     P_KickTdhfFem          absorption spectrum
             dynamics/tddft                    P_KickTddftFem         absorption spectrum, real-time TDDFT
             analysis/absorption_spectrum.py
  h2o/       analysis/stepCost                 P_StepCostFem          cost of MCTDHF vs TD-2RDM per derivative
  h2/        groundstate/scfProlate            P_GroundScfProlate     RHF on the prolate spheroidal grid
  n2/        groundstate/scf                   P_GroundScfProlate     RHF on the prolate spheroidal grid
             groundstate/mcscf                 P_GroundMcscf          MCSCF, 10 orbitals, TD-2RDM anchor
             dynamics/staticOrbs               P_LaserStaticOrbs      TDCI and TD-2RDM in a pulse along the bond
  c2h4/      groundstate/scf                   P_GroundScfFem         RHF on grid.fem
             analysis/integrals                P_IntegralsFem         integrals of HF orbitals + STOs
             dynamics/td2rdm                   P_KickTd2rdmFem        TD-2RDM in 10-16 active orbitals, kick spectrum
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
default 300). Outputs (`*.dat`, `*.in`, `*.bin`, `*.npz`, `*.log`) and the
`runs/` folders of generated inputs are git-ignored.

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

Order: `groundstate/scf`, then the correlated ground state with either
`groundstate/mcscf` (ground solver `groundSolver.mcscf.gradient`, 10
iterations, 2 s: E = −14.6168560449) or `groundstate/mctdhf` (imaginary-time
MCTDHF, about 50 s: E = −14.6168560112, 3.4e-8 higher from its dtau = 0.005
step bias), then `groundstate/td2rdm` (TD-2RDM initial state from the MCTDHF
ground state: field-free TD-2RDM propagation averaged over 40 a.u.), then any
`dynamics/<method>/<pulse>`. Both correlated runs start from the
Hartree-Fock 1s 2s of `../scf` plus hydrogen-like 2p (`orbsInit.load.nOrbsIn`)
and the CI eigenvector of these orbitals (`program.ciStartQ` for imaginary
time: from a unary CI vector it ends in a local minimum at E = −14.5956,
where the 2p shell correlates the 1s² core), and both write the same files
(the dynamics load `groundstate/mctdhf`; point `orbsInit.load.dir` and
`coeffsInit.load.dir` at `groundstate/mcscf` to use the MCSCF state).
Instead of the averaged state of `groundstate/td2rdm`, the TD-2RDM dynamics
can start from the MCSCF state itself with an anchored reconstruction: run
`groundstate/mcscf` with a `program.anchor` block (writes
`rdm3Bt01_01_02.in`), load orbitals and 2-RDM from it, and add
`"anchor": {"dir": "../../../groundstate/mcscf"}` to
`reconstruction.singlet` (`src/Reconstruction/Singlet/Anchor/A_Reconstruction_Singlet_Anchor.md`): the
ground state is then stationary under the TD-2RDM equations.

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
- `groundstate/scf`, `groundstate/mcscf`: HF (E = −128.5470989, 0.7 s) and
  9-orbital MCSCF ground states on a small grid (lmax 1, rmax 8). The MCSCF
  is the variational `groundSolver.mcscf.gradient` (program of
  `atoms/be/groundstate/mcscf`): HF 1s 2s 2p of `../scf` plus hydrogen-like
  3s 3p, 12 iterations, 3.5 s, E = −128.6818315, 0.135 Ha below HF (3s 3p
  occupations 5.7e-3, 3.5e-3). `dynamics/mctdhfLengthGauge` propagates it in
  a length-gauge pulse on that grid (restricted orbitals of declared m,
  total-M = 0 sector) at dt = 0.0025: at dt 0.01 the energy of the
  correlated state drifts by 7 mHa within 2 a.u. before the pulse (splitting
  error; 4th-order splitting is worse), at 0.0025 by 4e-5. Its `hhg.dat`
  holds the length form only.
- `groundstate/td2rdm`: TD-2RDM ground state without CI, moving orbitals
  (`P_GroundTd2rdm`, small grid of `../mcscf`): the Hartree-Fock determinant
  in HF 1s 2s 2p of `../scf` plus hydrogen-like 3s 3p, the interaction
  switched on over 20 a.u. (`method.mb.orbBased.adiabaticSwitching`, frozen
  HF mean field as the reference, `regularizationParameter` 1e-14), singlet
  NY with ccQ, split step dt 0.01, then 10 a.u. of free propagation; 4.5 min
  (8 threads). E(0) = −128.67801 (T = 40: −128.67765), 97% of the MCSCF
  correlation energy in these 9 orbitals; natural occupations of the
  correlating orbitals 5.7e-3 (three) and 3.5e-3, as in MCSCF; min eigenvalues of D and Q about −1e-5; the
  energy drifts by 3e-4 over the 10 a.u. With the default regularization
  1e-10 E(0) is −128.67437 (T = 20) and the occupations oscillate by 10%
  afterwards.
- The MCTDHF grid (rmax 160, 40 uniform elements of 4 bohr) is coarse at the
  nucleus: HF there is 5.8 mHa below the HF limit, MCSCF on that grid has
  the same offset (−128.6876610), so the correlation energy agrees with the
  small grid. Imaginary-time MCTDHF reaches the MCSCF state (same natural
  occupations to 1e-6) with a dtau bias of 2.9 mHa at dtau 0.02, 4.1e-4 at
  0.01 and 3.4e-5 at 0.005.

## H2 and N2: RHF on the prolate spheroidal grid

`P_GroundScfProlate` (`groundSolver.scf.stdImpl`, DIIS on the orbitals) on
`grid.prolate`, nuclei at the foci; σ and π orbitals (azimuthal channels up
to `grid.prolate.mmax`) come out of one Fock eigenproblem, and the program
labels every occupied orbital (σ/π, g/u).

- `h2/groundstate/scfProlate`: R = 1.4, ARPACK, 9 steps, 0.2 s, E =
  −1.133629571, the HF limit to 1e-11 (ximax 20, nEta 10, 5 × 12 ξ points).
- `n2/groundstate/scf`: R = 2.068, LOBPCG preconditioned by the exact
  (T + σ)⁻¹ of the prolate Laplacian with ACE, graded ξ elements
  (`grid.prolate.fedvr.grading` 2.0), 11 steps, 24 s: E = −108.9938256, the
  numerical HF limit −108.993826 to all its digits; orbital energies 1σg
  −15.6819, 3σg −0.6348, πu −0.6156 (πu above 3σg). The π channels use the
  axis-adapted odd-m basis of the grid (`src/Grid/A_Grid.md`, Prolate Grid);
  with the plain basis the same grid gave −108.996274, 2.4 mHa below the limit.

## N2: MCSCF and TD-2RDM against TDCI

- `n2/groundstate/mcscf`: 10 restricted orbitals (the 7 occupied RHF orbitals
  of `../scf` plus LCAO 1πg and 3σu), full CI of 14 electrons (14 400
  configurations), `groundSolver.mcscf.gradient`, 11 iterations, 36 s:
  E_el = −132.8362386, total −109.14184 (0.148 Ha below RHF); natural
  occupations 1πu 0.971, 1πg 0.032, 3σu 0.0096. With `program.anchor` it
  also writes the reference of an anchored TD-2RDM reconstruction
  (`rdm3Bt01_01_02.in`) and, with `nTransitions` 12, the transition RDMs of
  the 12 lowest Ritz states of the CI Hamiltonian in the Krylov space of
  Z|Ψ⟩ (60 Lanczos vectors; excitation energies 0.7651, 1.2135, 1.4375,
  1.9184, … Ha, the dipole-allowed 1Σu+ states): 80 s on three threads,
  208 MB of files.
- `n2/dynamics/staticOrbs`: N2 in the MCSCF orbitals (static), length-gauge
  pulse along the bond (`sysGauge.prolate.length`), TDCI
  (`P_LaserStaticOrbs.json`) or TD-2RDM (`td2rdmNy.json`,
  `td2rdmNyAnchor.json`, `td2rdmNySpectral.json`). `program.kick` replaces
  the pulse by a dipole kick exp(−iκZ) (TDCI and TD-2RDM; remove the
  `sysGauge` block). TD-2RDM costs about 0.05 s per time derivative on one
  thread (RK4, dt = outputStep/nTimeSteps = 0.04): 120 a.u. take about 10 min.

What the anchor of the reconstruction (`src/Reconstruction/Singlet/Anchor/A_Reconstruction_Singlet_Anchor.md`)
does for N2 (NY closure with ccQ, all numbers against the exact TDCI of the
same orbitals):

- Ground state without field: with the plain closure Q gets an eigenvalue of
  −3e-3 within 4 a.u. (the 1σ core holes) and the occupations drift; with the
  anchor the MCSCF state is stationary (energy to 1e-14, eigenvalues of D and
  Q unchanged).
- Linear response (kick κ = 0.01, 60 a.u.): the closure puts the bright
  state at 0.809 Ha instead of 0.769 and adds a spurious peak at 0.633 Ha;
  the dipole error equals the dipole itself (rms 0.0229 vs 0.0228). With the
  spectral anchor (12 transitions) the peaks are those of TDCI (0.769, 1.213,
  1.439, 1.92) and the dipole error is 0.0020; occupations within 3e-4.
- Weak resonant pulse (0.7651 Ha, E0 = 0.01 a.u., 10 cycles, 5% excited):
  TDCI absorbs 0.0380 Ha. Spectral anchor with D/Q purification: 0.0393 Ha,
  dipole within 6% up to t = 50, after that a growing phase error (half the
  dipole at t = 80–90).
- Strong resonant pulse (E0 = 0.04 a.u., two thirds of the population
  excited): TDCI absorbs 0.494 Ha. The NY closure without and with the plain
  anchor absorbs about 0.2 Ha and diverges at t = 65–70 (negative
  occupations); with purification it stays physical (min eig D −1e-5, Q
  −4e-4) but absorbs 0.263 Ha. Spectral anchor with purification: 0.563 Ha,
  dipole close to TDCI up to t = 50; afterwards both differ from TDCI by as
  much as the dipole itself, and the 1πg occupation stays at 0.04 instead of
  0.1–0.17. The spectral anchor is linear in D2 − D2_ref: far from the ground
  state the closure errors of second order remain, and without purification
  the spectral anchor alone diverges after strong driving (occupations
  above 1 that grow exponentially).

## Ethylene: TD-2RDM in a fixed active space beyond CI dynamics

C2H4 (C=C along z, molecule in the yz plane, 16 electrons) on grid.fem with a
fixed orbital basis: the grid enters once, through the integrals, and the
dynamics is a finite-basis problem in which TD-2RDM costs n^6 to n^7 per
derivative independent of the number of electrons, CI binomial(n, 6)^2 x
connections.

1. `groundstate/scf`: RHF in a box of ±16 bohr with a ghost center (charge 0)
   at the C=C midpoint that grades the mesh for the Rydberg functions (685 647
   nodes, 11 steps, 137 s): E = −78.1052105.
2. `analysis/integrals` (`P_IntegralsFem`, 93 min, 3 GB): h1, h2, dipoles and
   overlap of the 8 occupied RHF orbitals plus 39 raw STOs: C 2s 2p (ζ 1.6,
   0.9), C 3d (ζ 2.0), H 1s (ζ 1.2) and 3s 3p 3d (ζ 0.5) on the ghost center
   (`intsDiffuse.bin`).
3. `dynamics/td2rdm/make_inputs.py NAME NVIRT select=cis` writes `runs/NAME/`
   (Hamiltonian and input of `P_KickTd2rdmFem`): it projects the STOs on
   the complement of the occupied orbitals, orthogonalizes them (39
   virtuals), freezes C 1s and writes the active space of the 6 occupied
   valence orbitals and NVIRT virtuals (12 active electrons). `select=cis`
   takes the particle natural transition orbitals of the bright z-polarized
   CIS states first, state by state from the lowest one, and fills up with
   MP2 natural orbitals of the rest; `select=mp2` (default) takes MP2 natural
   orbitals only (MP2 with frozen C 1s: −0.2557 Ha in this basis).
4. `P_KickTd2rdmFem`: TD-2RDM (singlet NY with ccQ, static orbitals) from the
   Hartree–Fock determinant, the correlation switched on adiabatically over
   40 a.u. (`method.mb.orbBased.adiabaticSwitching`, from the Fock operator
   of the determinant; no CI anywhere: NY is exact for a determinant), then
   a z kick (κ = 0.01) and 80 a.u. of free evolution (RK4, dt 0.2).

The rows of the table: `python3 make_inputs.py d10 4 select=cis` (d12 6, d14
8, d16 10), then `cd runs/d10 && ../../P_KickTd2rdmFem.exe
P_KickTd2rdmFem.json`.

The V state (π→π*, B1u, polarized along C=C) is the one to compare with
experiment: band maximum 7.66 eV, vertical excitation about 7.9 to 8.0 eV. CIS
in the full 39-virtual space puts it at 0.2928 Ha = 7.97 eV (oscillator
strength 0.45). In a compact basis (±10 box, 22 STOs without diffuse
functions) it was at 0.3825 Ha = 10.4 eV, and an
exact 10-orbital calculation there gave 11.7 eV: the missing diffuse
functions were the main error. The excitation energies below come from a
matrix-pencil fit of the dipole signal; the median over model orders 8 to 30
is quoted, because a single order can split or shift peaks by 0.1 eV over
80 a.u. (it splits a peak of the exact signal at 40 modes).

| active orbitals | E_corr exact | TD-2RDM error | V state exact | V state TD-2RDM | CIS in the active space | TD-2RDM run |
|---|---|---|---|---|---|---|
| 10 | −22.66 mHa | −0.06 mHa | 0.3232 Ha = 8.79 eV | 8.84 eV | 8.41 eV | 1.8 min |
| 12 | −28.01 mHa | −0.10 mHa | 0.3139 Ha = 8.54 eV | 8.56 eV | 8.40 eV | 5 min |
| 14 | | E_corr −36.0 mHa | | 8.42 eV | 8.35 eV | 12 min |
| 16 | | E_corr −45.5 mHa | | 8.25 eV | 8.33 eV | 31 min |

The exact references are a matrix-free FCI ground state and a Lanczos
spectrum of D_z|0⟩ (Python on `runs/dNN/hamiltonian.npz`, not part of the
repository; 8 s at 10, 5.5 min at 12 orbitals; 14 orbitals would take hours
here). The second bright state agrees as well: 11.62 eV TD-2RDM against
11.67 eV exact at 10 orbitals, 11.45 against 11.38 at 12. The V state comes
down by 0.15 to 0.3 eV per two orbitals and is at 8.25 eV with 16, 0.3 eV
above the experimental vertical value; 18 orbitals would take about 80 min.
The minimal eigenvalue of D stays at −3e-3 to −4e-3 (no purification) and the
signal does not grow.

Orthonormalizing all CIS particle orbitals at once and truncating by singular
value (an earlier version of `select=cis`) lost the π* orbital
of the V state at 10 and 12 orbitals: CIS in those spaces gave 9.18 eV with a
seventh of the strength, and the exact 10-orbital V state was at 8.89 eV
(E_corr −2.6 mHa against −22.7 with the state-ordered selection).

`twoRdm=singletTr0` (explicit 1-RDM and trace-free 2-RDM, NY without ccQ) is
faster per derivative but not usable here: in the compact
basis at 10 orbitals it gives a spurious bright peak at 0.548 Ha (40% of the
strength of the main peak, independent of the switching time) and three
times more negative eigenvalues of D than the singlet storage with ccQ. The
singlet storage costs 0.23 s per derivative at 14 orbitals (8 threads).

## H2O: cost of TD-2RDM and MCTDHF on grid.fem

`h2o/analysis/stepCost/P_StepCostFem` times one time derivative of MCTDHF
(full CI) and TD-2RDM (NY with ccQ) for H2O (10 electrons) on the mesh of
`groundstate/scf` (353 367 nodes) with n moving orbitals (STOs; the cost does
not depend on the values). `python3 make_inputs.py 6 8 10 12 14 16` writes the
inputs to `runs/` (`./P_StepCostFem.exe runs/td2rdm_10.json`). Four threads:

| n | CI coefficients | MCTDHF | TD-2RDM | of which 2-RDM part |
|---|-----------------|--------|---------|---------------------|
| 8  | 3 136      | 5.2 s  | 5.1 s  | 0.01 s |
| 10 | 63 504     | 12.6 s | 12.6 s | 0.06 s |
| 12 | 627 264    | 26.9 s | 22.8 s | 0.16 s |
| 14 | 4 008 004  | 89.9 s | 42.5 s | 0.42 s |
| 16 | 19 079 424 | 497 s  | 65.9 s | 0.88 s |

Both methods share the grid part (one-body operator, n(n+1)/2 Poisson solves
for the pair potentials, mean fields); TD-2RDM replaces the CI vector by the
2-RDM (n^4 entries), and the CI matvec and the RDMs of the CI vector by the
commutators and the reconstructed collision operator. Up to n = 10 the grid
part is everything and both cost the same. From n = 12 the CI part of MCTDHF
grows combinatorially (CI matvec 2.2 s at 12, 23 s at 14, 191 s at 16, plus
the 2-RDM of the CI vector), the 2-RDM part of TD-2RDM stays below a second,
and TD-2RDM is 2.1 times faster at n = 14 and 7.5 times at n = 16 (memory
1.4 GB for both at 16). The grid part (about n^2) then sets the cost of
TD-2RDM: a derivative still takes a minute at n = 16, so a 50 a.u. kick
spectrum (thousands of derivatives) is about a day on this machine for
TD-2RDM and a week for MCTDHF.

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
H2O: SCF in 12 steps (22 with linear Hartree-potential mixing), 30 s on
two threads (466 Fock products; two LOBPCG guard vectors cost 686, tol 1e-6
three more SCF steps); the SCF mesh is the coarse one of the kick (the accurate mesh with 1.1 M nodes,
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
