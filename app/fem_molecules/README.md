# Molecules on the finite-element grid: RHF ground states and TDHF kick spectra

All-electron calculations for small closed-shell molecules on the 3D
finite-element grid `grid.fem`:

1. `P_FemScf` computes the restricted Hartree–Fock ground state.
2. `P_FemKickTdhf` kicks the ground state and propagates it with real-time
   TDHF (linear response).
3. `absorption_spectrum.py` turns the induced dipole into the polarizability,
   the absorption cross section and the dipole strength function, and checks
   the Thomas–Reiche–Kuhn (TRK) sum rule.

The inputs are for H₂ (R = 1.4 bohr) and H₂O (O at the origin, H at
(±1.4305227, 0, −1.1071202): r_OH = 1.80890 bohr, 104.52°, molecular plane xz,
C₂ axis z).

## The grid in short

- **Spectral elements.** The mesh is a tensor product of three 1D element
  partitions of the box, i.e. hexahedral elements. In every element the
  orbitals are Lagrange polynomials on `nLocals` Gauss–Lobatto–Legendre (GLL)
  nodes per edge (degree nLocals − 1), continuous across element faces, and
  zero on the box surface (Dirichlet).
- **GLL lumped mass.** GLL quadrature makes the mass matrix diagonal, so an
  orbital is its vector of nodal values and the inner product is
  Σᵢ conj(fᵢ) gᵢ wᵢ.
- **Graded mesh, nuclei at vertices.** On each axis the element size is
  `hMin` of the nearest nucleus there and grows geometrically (`growth`) up to
  `hMax`. Every nuclear coordinate is an element boundary on every axis, so
  every nucleus sits at an element vertex, where the C⁰ basis can represent
  the cusp.
- **Fast diagonalization (FDM).** The 1D operators of each axis are
  diagonalized once. Shifted kinetic solves (T + s)⁻¹ and the Poisson equation
  (free-space boundary values from a multipole expansion, `multipoleLmax`) are
  then direct and exact at the cost of a few matrix–matrix products. They give
  - the preconditioner of the LOBPCG eigensolver,
  - the Hartree and exchange potentials,
  - the preconditioner of the Cayley time step.
- **Galerkin point nuclei.** The nuclear potential is the Galerkin solution of
  the Poisson equation for point charges. It is finite at the nucleus node
  (roughly −30 … −50 Z/hMin for nLocals 5–6, depending on the mesh; see
  *Nuclear spike* below) and integrates correctly against smooth densities.
  It is **not variational**: on coarse meshes the energy lies below the exact
  (here: Hartree–Fock-limit) value and converges from below as the mesh is
  refined, mostly through the under-resolved core of heavier atoms.
- **Dynamics.** The field-free one-body part h₀ = T + V_nuc is propagated by
  the unitary Cayley (Crank–Nicolson) step, solved with GMRES preconditioned by
  the exact FDM inverse of 1 + iΔt T/2, so the stiff nuclear spike does not
  limit that step. The mean field is integrated by RK4 in a Strang splitting.
  With `exactFieldFreeMeanFieldQ` (set in both kick inputs) the Cayley step
  also carries the frozen Fock mean field of the ground state (built before
  the kick), so RK4 only integrates its change; without it the O core of H₂O
  limits the time step to Δt ≈ 0.0025 (see *H₂O dynamics* below).

Details: `src/Grid/A_Grid.md` (FEM grid), `src/SysPotential/A_SysPotential.md`
(nuclei), `src/GroundSolver/A_GroundSolver.md` (SCF, ACE, LOBPCG
preconditioner), `src/Method/A_Method.md` (Cayley step).

## Files

| File | Purpose |
|---|---|
| `P_FemScf.f90` | RHF ground state. Writes the orbitals `orbXX_YY.in` and the mesh fingerprint `femMesh.dat` |
| `P_FemKickTdhf.f90` | Kick e^{iκ r_a} of every orbital, then TDHF (restricted MCTDHX with nOrbs = nBodies). Writes `dipole.dat` (t, ⟨x⟩, ⟨y⟩, ⟨z⟩), `observables.dat` (t, E_el, E_tot, N), the checkpoint files `orbXX_YY.in`, `coeffs.in`, `checkpoint.dat` and, with the mean-field reference, `meanFieldReference.dat` (the ground state it is built from) |
| `absorption_spectrum.py` | Spectrum from `dipole.dat` (pure Python 3, no numpy) |
| `h2_scf.json`, `h2_kick.json` | H₂ on the mesh of the FEM simulation tests: ±10 box, nLocals 5, hMin 0.3, N = 130 331 nodes |
| `h2o_scf.json` | Accurate H₂O ground state: ±10 box, nLocals 6, hMin 0.05 (O) / 0.2 (H), growth 1.8, hMax 2, N = 1 101 464 |
| `h2o_scf_coarse.json`, `h2o_kick.json` | H₂O dynamics mesh: ±10 box, nLocals 5, hMin 0.1 (O) / 0.3 (H), growth 1.8, hMax 2, N = 353 367 (both files have the same grid block); the kick runs 10 a.u. at Δt = 0.1 with the mean-field reference |

## How to run

Build the programs (release), then run each step in its own directory,
outside the repository or below `app/fem_molecules/runs/` (git-ignored; the
programs and the script, with its default output names, write only `*.dat`
and `*.in` files, and the commands below log to `*.log`: all three are
ignored anywhere in the repository). Use two threads on this two-core
machine.

```bash
source setBuildVars.sh -c gnu -t release
cody build app                    # or: cmake -B build && cmake --build build -j2 --target apps
export OMP_NUM_THREADS=2 CODY_PRETTY_PRINT=OFF
APP=$PWD/app/fem_molecules
RUNS=$APP/runs                    # any directory

# H2: ground state (5 s), 50 a.u. of TDHF after a z kick (about 15 min), spectrum
mkdir -p $RUNS/h2_scf && cd $RUNS/h2_scf
$APP/P_FemScf.exe $APP/h2_scf.json > scf.log
mkdir -p $RUNS/h2_kick && cd $RUNS/h2_kick
cp ../h2_scf/orb*.in ../h2_scf/femMesh.dat .
$APP/P_FemKickTdhf.exe $APP/h2_kick.json > kick.log
python3 $APP/absorption_spectrum.py dipole.dat --pade

# H2O: accurate ground state (about 10 min, 1.1 GB)
mkdir -p $RUNS/h2o_scf && cd $RUNS/h2o_scf
$APP/P_FemScf.exe $APP/h2o_scf.json > scf.log

# H2O dynamics: ground state on the dynamics mesh (about 4 min), then a 10 a.u.
# kick demo (about 50 min; see "H2O dynamics" below before planning longer runs)
mkdir -p $RUNS/h2o_scf_coarse && cd $RUNS/h2o_scf_coarse
$APP/P_FemScf.exe $APP/h2o_scf_coarse.json > scf.log
mkdir -p $RUNS/h2o_kick_z && cd $RUNS/h2o_kick_z
cp ../h2o_scf_coarse/orb*.in ../h2o_scf_coarse/femMesh.dat .
$APP/P_FemKickTdhf.exe $APP/h2o_kick.json > kick.log
```

The JSON file is the only command-line argument, and every file is read from
and written to the working directory.

- **Mesh check.** The kick program loads the orbitals of the SCF run. They
  must come from the identical mesh: the program compares its grid with
  `femMesh.dat` and stops if they differ. A kick input therefore repeats the
  grid block of its SCF input.
- **Spectrum of H₂O.** H₂O is not isotropic. Run the kick three times, with
  `program.kickDirection` 1, 2 and 3 (x, y, z) in separate directories, and
  average: `python3 $APP/absorption_spectrum.py x/dipole.dat y/dipole.dat
  z/dipole.dat` gives (α_xx + α_yy + α_zz)/3.
- **Checkpoints and restarts.** With `program.saveStep` (a multiple of
  `program.outputStep`) the state is saved during the run and at its end; the
  checkpoints overwrite the SCF orbitals in the run directory. To continue,
  set `program.startTime` to the time in `checkpoint.dat` (and `endTime`
  further out) and start the program again in the same directory. There is no
  second kick; `dipole.dat` and `observables.dat` keep their lines up to
  startTime and are continued. The program stops if `checkpoint.dat` does not
  hold startTime. A continued run reproduces the uninterrupted one bit for
  bit (checked on H₂, Δt = 0.05, restart at t = 0.2). With
  `exactFieldFreeMeanFieldQ` (the shipped inputs) the first run saves the
  ground state the mean-field reference is built from
  (`meanFieldReference.dat`), and the continuation rebuilds the same
  reference from it; keep that file in the run directory (the program stops
  without it). Rebuilding the reference from the checkpoint state instead
  would change the splitting: on H₂O the energy then grew quadratically
  (2.7e-6 (t − 5)² Ha after a restart at t = 5).
- **Absorber (optional).** An `absorber.fem.cosinus` block (`width`, `order`;
  see `src/Absorber/A_Absorber.md`) is applied after every time step. With
  `absorber.keepOrthonormalQ` the absorbed norm moves into the CI vector; such
  a run then needs `coeffsInit.load` for continuations. The absorber removes
  oscillator strength, so the TRK sum is no longer N.

### Program keys

| Key | Program | Default | Meaning |
|---|---|---|---|
| `program.convThresh` | SCF | 1e-8 | stop when \|ΔE_el\| between two SCF steps is below this … |
| `program.epsThresh` | SCF | 0 (off) | … and max \|Δε\| of the occupied orbital energies is below this (they converge like √ΔE) |
| `program.maxIter` | SCF | 50 | at most this many SCF steps. Otherwise the program says so, saves the last orbitals (continue with `orbsInit.load`) and stops with an error |
| `program.kappa` | kick | 1e-3 | kick strength κ (linear response: κ ≪ 1) |
| `program.kickDirection` | kick | 3 | 1, 2, 3 = x, y, z |
| `program.timeStep` | kick | 0.05 | Δt |
| `program.outputStep` | kick | 0.1 | output interval, a multiple of Δt; π/outputStep is the Nyquist frequency of the spectrum |
| `program.startTime`, `program.endTime` | kick | 0, 10 | endTime − startTime must be a multiple of outputStep |
| `program.saveStep` | kick | 0 (end only) | checkpoint interval |

The SCF uses `groundSolver.scf.stdImpl` with the Hartree potential of the
total density, adaptively compressed exchange (`aceQ`, needed for LOBPCG on
this grid), linear mixing of the Hartree potential (`mixing.linear.alpha`;
0.7 for H₂O) and `diagonalizerList.lobpcg` preconditioned by
(T + max(−θ, `groundSolver.precShiftMin`))⁻¹. `nEvals` must be at least the
number of doubly occupied orbitals. The start orbitals come from
`orbsInit.fem.sto` (Slater-type functions on the nuclei; for H₂O O 1s, 2s, 2p
mixed with H 1s) or, to continue, from `orbsInit.load`. They are
orthonormalized first.

The kick uses `propagator.splitStep.order2` with `exactPartScheme: "strang"`,
`method.mb.orbBased.exactFieldFreeQ` and `exactFieldFreeMeanFieldQ`: per time
step four Cayley solves per orbital (Δt/4 each, GMRES tolerance
`exactFieldFreeTol` 1e-12) of h₀ + U, U = J₀ − K₀ the Hartree–Fock mean
field of the ground state (built once before the kick from the SCF orbitals,
exchange as ACE: 1 + n(n+1)/2 Poisson solves for n orbitals, and saved for
continued runs), and one RK4 step of the change of the mean field, W − U.

The Cayley step has no imaginary time, no exponential (phi-function)
propagators and no spectral cut (`src/Method/A_Method.md`). A
`sysGauge.fem.length` block (length gauge, sin² pulse) adds a laser; with
`program.kappa` 0 the program then propagates the ground state in the pulse
without a kick (checked on H₂, 1 a.u.).

## Results on this machine

Intel Pentium G4560 (2 cores, 4 hardware threads), `OMP_NUM_THREADS=2`.
Other users' jobs kept the load at 5–10, so the timings vary by a factor of
about 2.

This section is the one home of these measurements: the module docs
(`src/Method/A_Method.md`, `src/Propagator/A_Propagator.md`) keep the
contracts and the values their tests assert, and refer here for drifts,
iteration counts and timings of the app runs.

### H₂ (R = 1.4 bohr, N = 130 331)

| Quantity | Result | Reference |
|---|---|---|
| E_el, E_nn | −1.847953860878, 0.714285714286 Ha | |
| E_tot | −1.133668146592 Ha | HF limit −1.133629571: 3.9e-5 below (Galerkin nuclei) |
| ε(1σg) | −0.5946710914 Ha | HF limit −0.594658 |
| SCF | 13 steps, 186 Fock products, 5.3 s | |
| kick κ = 1e-3 along z: E(0+) − E(0) | 9.99998724e-7 Ha | Nκ²/2 = 1e-6 |
| d⟨z⟩/dt(0+) on the grid | 1.99999747e-3 | Nκ = 2e-3 (TRK) |
| TDHF, Δt = 0.05 with the mean-field reference (`h2_kick.json`), 50 a.u. | \|E_el(t) − E_el(0+)\| ≤ 3.1e-11 Ha (no drift; 3e-5 of the kick energy); \|N − 2\| ≤ 5.8e-11; ⟨x⟩, ⟨y⟩ < 1e-10; induced dipole within 4.5e-4 of its largest value of the plain Δt = 0.025 run below (8.4e-5 over the first 5 a.u.) | |
| cost | 1000 steps in 851 s: 0.85 s per step, 17 s per a.u. (14 min with setup and output); 4 Cayley solves per step of 6 GMRES iterations each; reference built once in 0.03 s. A Cayley solve at Δt/4 = 0.0025–0.005 (5 GMRES iterations, plain split as in the FEM simulation tests) took about 0.12 s, of which 6 FDM solves (the preconditioner) of 14 ms each | |
| for comparison: plain split (no reference), Δt = 0.025, 50 a.u. | \|E_el(t) − E_el(0+)\| ≤ 6.8e-9 Ha (bounded); N − 2 = 3.0e-9 at t = 50; 2000 steps in 1520 s (30 s per a.u.) | |

Time step: without the mean-field reference Δt = 0.05 lets the energy
drift (1.2e-7 Ha in the first a.u. of the ground state, 4e-7 Ha in 5 a.u.
after the kick; 3.5e-9 at Δt = 0.025). With the reference the energy stays
within 3.1e-11 Ha of E(0+) for 50 a.u. at Δt = 0.05 (1.3e-10 with U from the
kicked state), and the induced dipole agrees with the plain Δt = 0.025 run
to 8.4e-5 of its largest value over the first 5 a.u. (the plain Δt = 0.05
run: 1.4e-4), so `h2_kick.json` uses Δt = 0.05.

Spectrum (`absorption_spectrum.py dipole.dat --pade`) of the run of
`h2_kick.json`:

| Quantity | Result |
|---|---|
| lowest Σu⁺ excitation (B state) | **0.4683 Ha = 12.74 eV** (Padé; 0.4691 from the first 45 a.u.). The Fourier transform resolves it with τ = T/2 (0.473 Ha); with the default τ = T/4 it merges with the lines above (one peak at 0.488 Ha) |
| its oscillator strength | f_z ≈ 0.80–0.84 (Padé line, τ = T/2 transform), i.e. f = f_z/3 ≈ 0.27 orientation-averaged |
| next lines | 0.587, 0.751, 0.803 Ha (Padé): box states at and above the TDHF ionization threshold −ε(1σg) = 0.595 Ha; they depend on the box and on the Padé fit (the plain Δt = 0.025 run gives 0.595, 0.725, 0.833) |
| static polarizability α_zz | 6.44 a.u. (Padé), 6.33 a.u. (damped transform) |
| TRK | initial slope of the samples 1.99930 (N = 2); ∫S dω = 0.49 below 0.5 Ha, 1.59 below 1 Ha, 1.87 below 2 Ha, 1.95 below 5 Ha, 1.992 up to the Nyquist frequency 31.4 Ha |

The plain split at Δt = 0.025 (no reference) gives the same damped
transform (peak 0.48784 against 0.48785 Ha, α 6.3299 against 6.3290) and
the lowest Padé line at 0.4678 Ha; U from the kicked state gives the same
lines as U from the ground state to 1e-4 Ha. The Padé peaks are not always stable: the
first 40 a.u. alone give 0.461 and 0.548 Ha (the plain run: 0.437 and 0.527)
instead of one line at 0.468. The script therefore repeats the lowest Padé
peaks from 90 % and 80 % of the samples; trust only lines that stay.

### H₂O ground state

O at the origin, H at (±1.4305227, 0, −1.1071202). The orbital energies are
1a₁, 2a₁, 1b₂, 3a₁, 1b₁ in the usual notation; μ points from O to the H
atoms (−z).

| Quantity | accurate mesh, `h2o_scf.json` | dynamics mesh, `h2o_scf_coarse.json` |
|---|---|---|
| N (nodes per axis) | 1 101 464 (119, 89, 104) | 353 367 (79, 63, 71) |
| E_el | −85.263850350971 Ha | −85.312953232287 Ha |
| E_nn | 9.194689103233 Ha | 9.194689103233 Ha |
| E_tot | −76.069161247738 Ha | −76.118264129054 Ha |
| ε | −20.5646840, −1.3526620, −0.7178110, −0.5849631, −0.5105194 | −20.5856191, −1.3536907, −0.7176764, −0.5849253, −0.5103531 |
| μ_z | −0.7790262 a.u. (1.980 D) | −0.7779986 a.u. (1.977 D) |
| SCF | 20 steps, 576 Fock products, 513 s (26 s per step), 1.1 GB | 22 steps, 735 Fock products, 200 s (9 s per step), 0.38 GB |

Comparison on the accurate mesh:
- **numpy prototype, same mesh:** E_tot −76.0691612, ε −20.5647675,
  −1.35268324, −0.71780699, −0.58500693, −0.51055356. The energy agrees to
  the printed digits. The prototype stopped at \|ΔE\| = 8e-8, while its ε
  still changed by up to 3e-4 per step, so the ε differences of up to 8e-5
  are its convergence error.
- **HF limit:** about −76.0675 Ha. The accurate mesh lies 1.7 mHa below it:
  the O core is resolved to 1–2 mHa, from below, because Galerkin nuclei are
  not variational. The dipole, 0.779 a.u., matches the HF-limit value
  0.78–0.79 a.u.

The dynamics mesh lies 51 mHa below the HF limit, mostly in the O 1s level,
which is 21 mHa too deep. Its valence ε agree with the accurate mesh to 1e-3
Ha and its dipole to 0.13 %.

### H₂O dynamics: time step and cost

The kick works as for H₂ (κ = 1e-3 along z, dynamics mesh):
- E(0+) − E(0) = 4.99990773e-6 against Nκ²/2 = 5e-6;
- d⟨z⟩/dt(0+) on the grid = 9.99981675e-3 against Nκ = 1e-2.

**Why the mean-field reference.** In the O core h₀φ and the mean field Wφ
are both large and cancel only in their sum (⟨1s|h₀|1s⟩ = −33.1 Ha and
⟨1s|J − K|1s⟩ = +12.5 Ha against ε(1s) = −20.6 Ha; valence orbitals: about
−7 and +6.4 Ha). With the plain split (h₀ exact by Cayley, W by RK4) the
ground state itself is not stationary under either part, so the TDHF energy,
which should stay at E(0+), rises steadily, by the same amount with and
without the kick. `method.mb.orbBased.exactFieldFreeMeanFieldQ` (set in
`h2o_kick.json`) moves a frozen Fock mean field U = J₀ − K₀ into the exact
Cayley part (exchange as ACE: 16 Poisson solves once, 1–2 s) and leaves only
W − U to RK4 (`src/Method/A_Method.md`, "Frozen Mean-Field Reference"). The
equations of motion are the same; the ground state is an eigenstate of the
exact part and the RK4 part vanishes on it (2.6e-11 relative at t = 0;
without the reference its orbital derivative has the norm 4.8).

**Which U.** `P_FemKickTdhf` builds U from the ground state before the kick,
so that W − U is the pure response, O(κ) and oscillating. Built from the
kicked state (what Method does by default at the first step), U leaves a
static O(κ) part in W − U (the exchange change of the kicked orbitals),
which drifts the energy, 170 times faster at Δt = 0.05:

| Δt | plain split | U from the kicked state | U from the ground state (shipped) | run (a.u.) | GMRES iterations per solve with U (max) |
|---|---|---|---|---|---|
| 0.1 | – | 1.2e-6 Ha/a.u. | 2.8e-8 Ha/a.u. (6.5e-9 over the 10 a.u. of the demo) | 0–1 | 9.6 (11) |
| 0.05 | 8.7 Ha/a.u. (stopped at t = 0.3) | 3.2e-7 Ha/a.u. | 1.8e-9 Ha/a.u. | 0–1 (ground state: 0–0.5) | 9.2 (10) |
| 0.02 | – | 3.1e-8 Ha/a.u. | −9e-11 Ha/a.u. | 0–0.4 | 8.2 (9) |
| 0.01 | 6.8e-2 Ha/a.u. | – | – | 0–0.48 | (plain: 7.6) |
| 0.005 | 5.6e-3 Ha/a.u. | 5.4e-9 Ha/a.u. | – | 0–0.04 | 6.8 (8) |

(energy drift of E_el after the kick, κ = 1e-3, so the kick energy is 5e-6
Ha; slope of a linear fit over the run)

- **Ground state without kick:** the plain split at Δt = 0.005 drifts by the
  same 5.6e-3 Ha/a.u., and ⟨z⟩ moves by 1–2e-7 within 0.04 a.u. With the
  reference at Δt = 0.05, |E − E(0)| ≤ 6.6e-11 Ha over 0.5 a.u. and ⟨z⟩
  stays constant to 1e-9.
- **U from the ground state:** the energy jumps once in the first step
  (−5e-9, −1.3e-8, +1.5e-8 Ha at Δt = 0.02, 0.05, 0.1, where W − U is not
  yet small) and then stays flat at Δt ≤ 0.05; at Δt = 0.1 it creeps by
  2–5e-9 Ha per a.u. after the first a.u. (demo below).
- **U from the kicked state:** the drift scales like κ² (κ = 3e-3: 9.0 times
  the increment of κ = 1e-3 at t = 0.5, 3.4 % of the kick energy in both)
  and roughly like Δt² (E(0.4) − E(0+) = 5.3e-7, 1.35e-7, 1.5e-8 Ha at Δt =
  0.1, 0.05, 0.02). A U built from a later state is worse: rebuilt from the
  state at t = 5 (a continued run before `meanFieldReference.dat` existed),
  the energy grew quadratically, 2.7e-6 (t − 5)² Ha.
- **One trajectory:** at Δt = 0.02 the induced dipoles of both references
  agree within 7e-5 (relative) from t = 0.06 on and within 3e-5 from
  t = 0.08 (an offset of 1e-7 a.u., 2e-5 of the largest response, from the
  first step of the ground-state reference dominates at t = 0.02–0.04).
  Against Δt = 0.02 up to t = 0.4, Δt = 0.05 lies within 4.9e-4 (ground
  state) and 5.9e-4 (kicked), Δt = 0.1 within 2.1e-3 and 4.3e-3. The plain
  split at Δt = 0.005, once its spurious ground-state motion (the run
  without kick) is subtracted, agrees with the kicked-state reference at the
  same Δt to 2–3e-5 at t = 0.01–0.04; without the subtraction it is off by
  4–9e-4.
- **Not the RK4** (plain split): with the 8th-order GSL stepper `rk8pd` for
  the mean field the energies at Δt = 0.01 are the same to 1e-8.

**Cost** (dynamics mesh, 2 threads; other users' jobs kept the load at 4–9,
so the times vary by a factor of 2). A step makes 20 Cayley solves (5
orbitals × 4) of about 9–10 GMRES iterations and 4 mean-field evaluations
(15 pair Poisson solves each); the reference adds about 8 ms (10 %: 88
against 80 ms) per GMRES iteration (one multiplication, five inner products
and one pass over the five ACE vectors). Measured: 23–45 s per step at
Δt = 0.05 (460–900 s per a.u.), 24–57 s per step at Δt = 0.1 (240–570 s per
a.u.).

**Demo** (`h2o_kick.json` as shipped: κ = 1e-3 along z, Δt = 0.1, U from
the ground state; run to t = 5 and continued from the checkpoint to t = 10;
load 4–5):
- **Cost:** 100 steps in 2753 s of propagation, i.e. 27.5 s per step and
  275 s per a.u. (46 min in total for the two runs); 2000 Cayley solves of
  9.6 GMRES iterations on average (at most 11); the reference was built in
  0.8 s.
- **Conservation:** E_el − E_el(0+) = 7.1e-8 Ha at t = 5 and 9.0e-8 Ha at
  t = 10 (1.8 % of the kick energy). The increment per a.u. falls from
  3.3e-8 in the first a.u. to 2–5e-9 from t = 4 on, the same before and after
  the restart at t = 5. N − 10 = 2.4e-10 at t = 10; ⟨x⟩, ⟨y⟩ < 3e-12
  (symmetry).
- **Earlier demo for comparison** (U from the kicked state, Δt = 0.05,
  5 a.u., 63 min at load 6–9): E_el − E_el(0+) = 2.0e-6 Ha at t = 5, growing
  linearly. Its induced dipole agrees with the new demo to 1.6e-3 of the
  largest value, the Δt = 0.1 error. Continued from t = 5 with U rebuilt from
  the checkpoint state, its energy grew quadratically (to 3.4e-6 Ha at
  t = 5.7), which is why continued runs now rebuild U from
  `meanFieldReference.dat`.
- **Induced dipole** d_z(t) = ⟨z⟩(t) − ⟨z⟩(0): 3.24e-3 (t = 0.5), 4.75e-3 (1),
  4.94e-3 (1.5), 4.36e-3 (2), 3.36e-3 (2.5), 2.23e-3 (3), 1.17e-3 (3.5),
  2.6e-4 (4), −4.0e-4 (4.5), −8.7e-4 (5), −1.10e-3 (5.5), −1.17e-3 (6),
  −1.10e-3 (6.5), −9.5e-4 (7), −7.7e-4 (7.5), −5.5e-4 (8), −3.8e-4 (8.5),
  −2.7e-4 (9), −1.7e-4 (9.5), −7.8e-5 (10): the valence response, with a
  maximum at t ≈ 1.5 and a shallow minimum at t ≈ 6, damped by the
  dephasing of many lines.
- **Spectrum** (`absorption_spectrum.py dipole.dat --pade`), necessarily
  coarse (Fourier resolution 2π/T = 0.63 Ha):
  - The damped transform (τ = T/4) shows one broad band peaking at 1.12 Ha.
  - The Padé approximant shows one band at 0.87 Ha (0.88 and 0.89 Ha from 90 %
    and 80 % of the samples) and α_zz(0) = 8.5 a.u., in the range of the HF
    polarizabilities of H₂O (8–9 a.u.).
  - Separating the valence lines needs T of 30–50 a.u.
  - TRK: ∫S dω up to the Nyquist frequency 31 Ha is 8.63 of N = 10. The rest
    (O 1s) lies higher, which is also why the initial slope estimated from
    the 0.1-a.u. samples, 8.42, falls short of the grid value 9.9998 at
    t = 0+.

**Estimate for a 50 a.u. spectrum** (per kick direction; the
orientation-averaged spectrum needs three):
- **Δt = 0.1:** 500 steps, about 3.8 h at the load of the demo (27.5 s per
  step) and up to 8 h under heavy load, i.e. 12–24 h for three directions.
  At the demo's rate the energy rises by about 2.5e-7 Ha over 50 a.u. (5 %
  of the kick energy).
- **Δt = 0.05:** twice the cost, with 15 times less drift and induced
  dipoles about 4 times closer to the converged ones (4.9e-4 against 2.1e-3
  relative to Δt = 0.02 over the first 0.4 a.u.).
- **Without the reference:** Δt = 0.0025 (100–200 h per direction) still
  left a drift of about 13 mHa over 50 a.u.

## Mesh convergence

| H₂O mesh | nLocals | hMin O / H | N | E_tot | below the HF limit |
|---|---|---|---|---|---|
| `T_Fem_04` (±8, hMax 3) | 4 | 0.1 / 0.3 | 72 380 | −76.40376 | 336 mHa |
| dynamics (±10, hMax 2) | 5 | 0.1 / 0.3 | 353 367 | −76.11826 | 51 mHa |
| accurate (±10, hMax 2) | 6 | 0.05 / 0.2 | 1 101 464 | −76.06916 | 1.7 mHa |

- **The O center sets the accuracy.** Its hMin and nLocals decide the
  error, which comes from the core. For mHa accuracy at first-row nuclei use
  nLocals 6–8 and hMin ≈ 0.05 or less.
- **H centers are easy (Z = 1).** hMin 0.2–0.3 suffices: H₂ with nLocals 5
  and hMin 0.3 is 0.04 mHa below the HF limit.
- **Errors come from below.** With Galerkin nuclei a lower energy does not
  mean a better mesh. Compare meshes by their convergence toward a limit.
- **Valence properties converge first.** On the dynamics mesh the valence ε
  are within 1e-3 Ha and the dipole within 0.13 % of the accurate mesh.
- **Grading and box.** growth 1.6–1.8 and hMax 2–3 set the resolution of the
  valence and diffuse regions. A ±10 box suffices for neutral ground states.
  The multipole boundary values of the Poisson solve need the density to
  vanish at the faces; a one-time warning reports a box that is too small.
- **Spectra.** Without an absorber the box discretizes the continuum above
  −ε_HOMO. Only states that fit into the box (for H₂ the B state) are
  independent of it; the higher lines are box states.
- **Nuclear spike.** The Galerkin potential at the nucleus node measures
  −40.9 Z/hMin on the accurate mesh (nLocals 6) and −27.6 Z/hMin on the
  dynamics mesh (nLocals 5). It does not limit the Cayley time step, only the
  GMRES iterations (7–8 at Δt/4 = 0.0025, 9–10 at Δt/4 = 0.0125).
