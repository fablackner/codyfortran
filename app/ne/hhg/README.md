# Ne: strong-field MCTDHF with 9 orbitals

Ne (10 electrons) with 9 restricted spatial orbitals 1s 2s 2p 3s 3p (all m of
the p shells) in a z-polarized 800 nm two-cycle pulse, MCTDHF (full CI in the
9 orbitals, 15876 configurations). The programs are the beryllium ones
(`app/be/relax/P_ImagTimeMctdhf`, `app/be/laser/P_LaserMctdhf`); observables
as there (`observables.dat`: E(t), A(t), ⟨z⟩, dipole acceleration −Z⟨z/r³⟩,
electrons in the box, energy; `natocc.dat`).

```bash
source setBuildVars.sh -c gnu -t release
cody build app
cd app/ne/hhg
python3 make_inputs.py                  # inputs/stageA.json, stageB.json, mctdhf_strong2c.json
OMP_NUM_THREADS=4 ./run.sh inputs runs  # or e.g. ./run.sh inputs runs mctdhf_strong2c
```

| step | what it does |
|---|---|
| `relaxA` | imaginary-time MCTDHF ground state at lmax = 1 (exact: the orbitals of the ¹S state keep a pure l) |
| `relaxB` | embedding into the propagation grid, one check of the energy |
| `mctdhf_strong2c` | MCTDHF in the pulse |

## Parameters

| | |
|---|---|
| pulse | F(t) = F₀ cos(ωt) sin²(ωt/4), ω = 0.057, F₀ = 0.1194 (5×10¹⁴ W/cm², below the barrier-suppression intensity 8.6×10¹⁴), velocity gauge; quiver radius 37 a.u., cutoff Ip + 3.17 Up = 4.3 a.u. (harmonic 75) |
| radial grid | FEDVR, 40 elements × 20 points, rmax = 160 a.u. (761 points); cos^¼ absorber from 120 a.u. |
| angular grid | lmax = 40, \|m\| ≤ 1 (`grid.ylm.mmax`; the m of every orbital is conserved) |
| propagation | split step, exact field-free h0 with spectral cut 1000 a.u. and Strang composition, RK4 for the interaction, dt = 0.0025 (the fastest accurate MCTDHF setup of the Be benchmark, `app/be/README.md`) |
| ground state | dτ from 0.02 down to 0.001, non-linear sub-steps ≤ 0.002 |

Radial grid: the exactly diagonalized hydrogenic Ne 1s eigenvalue is −50.0000
with 4 a.u. elements of 20 points, −49.999 with 2 a.u. × 14 (1041 points) and
−49.997 with 1 a.u. × 10 (1441 points). The imaginary-time Strang fixed point
is biased by O(dτ²): ~6 mHa in the HF energy at dτ = 0.005, hence dτ = 0.001.

## Cost

Pentium G4560 (2 cores, 4 threads): 2.1 s per time step, 88 000 steps for the
pulse (≈ 2 days). Four non-linear derivatives per RK4 step take 0.38 s each
(0.28 s in the definite-m node quadrature of the correlation operator,
`src/Method/A_Method.md`; 3.7 s with the generic Gaunt products), four exact
field-free propagations 0.09 s each.

The orbitals keep their m exactly (z-polarized field). Round-off admixtures
of other m (seeded by the Löwdin steps of the relaxation, ~1e-16 of the norm
in the 3p orbitals) are removed by `ProjectOrbitalsOnDefiniteM`
(`app/be/common/BeCommon.inc`) after loading and after every Löwdin step;
otherwise the orbital projector spreads them into all channels and the exact
part cannot skip the empty m channels (≈ 0.5 s more per step). The correlation
operator tolerates admixtures up to 1e-12 of the squared norm; above that it
falls back to the generic operator (≈ 11 s per step).
