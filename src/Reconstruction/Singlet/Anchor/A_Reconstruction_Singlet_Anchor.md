# Anchored Singlet Reconstruction (`reconstruction.singlet.anchor`)

An optional extension of the singlet family (`src/Reconstruction/A_Reconstruction.md`).
`Reconstruction_Singlet_Fabricate` calls `Reconstruction_Singlet_Anchor_Fabricate`
when an `anchor` block exists; the extension then wraps the closure
`Reconstruction_Singlet_FillRdm3aab` of the scheme (adds delta3), with
`nTransitions` > 0 also `Reconstruction_AddCollisionOp` of the family (adds the
spectral part), and binds `Reconstruction_Setup` and
`Reconstruction_SetOrbitals`. Nothing outside this folder refers to it.


A closure R is exact for no correlated state, so a correlated ground state
(MCSCF, imaginary-time MCTDHF) is not stationary under the TD-2RDM equations:
for N2 (10 orbitals, NY with ccQ) the natural occupations drift away from the
ground state within a few a.u. without any field, and D and Q get negative
eigenvalues. The anchor adds the closure error at a reference state,

```
D3 = F(D2) = ccQ( R(D2) + T_S delta3 ),    delta3 = D3_ref - R(D2_ref),
```

with D2_ref, D3_ref the exact RDMs of the reference (`rdm2Bt01_02.in` of
`Coeffs_SaveTwoRdm`, `rdm3Bt01_01_02.in` of
`Reconstruction_Singlet_Anchor_SaveReference`, in
`anchor.dir`). At the reference the reconstructed 3-RDM is exact (ccQ leaves
the exact 3-RDM unchanged), so the reference is a stationary point of the
TD-2RDM equations; R still describes how D3 changes with D2. ccQ acts after
the anchor, so the contractions of D3 are those of D2 and energy conservation
is unaffected.

delta3 is a fixed kernel in the reference orbitals (`Orbs_orbs` at
`Reconstruction_Setup`). T_S transports it to the current ones: ket indices
with S, bra indices with conjg(S), S(i, j) = <phi_i(t)|phi_ref,j>, from the
orbitals that TD-2RDM passes to the generic `Reconstruction_SetOrbitals` before
every `Reconstruction_AddCollisionOp` (also in the split-step gauge where the
orbitals rotate with exp(-i h t); static orbitals: identity). The run must
start from the reference orbitals.

#### Spectral anchor (`anchor.nTransitions` = K > 0, static orbitals)

The anchor does not change the linear response of the closure, and that
response decides the excitation energies: for N2 NY puts the bright 1Σu+
state at 0.809 Ha instead of 0.769 and adds a spurious peak, so a resonant
pulse is detuned (about 0.2 of the 0.494 Ha that TDCI absorbs). The spectral
anchor adds a linear correction that makes chosen excitations exact
eigenmodes of the linearized TD-2RDM equations. From the transition RDMs
T2_k = <ref|..|chi_k>, T3_k of K excited states (`trdm2Bt01_02.in`,
`trdm3Bt01_01_02.in`, ket-first, transition index last) the setup forms the 2K
Hermitian perturbations v_a = t + t^+, t = T2_k and i T2_k, their exact 3-RDM
changes, and the errors e_a of F (the change of F along v_a by a central
difference). The correction

```
D3 += sum_a c_a e_a,    c_a = Re<u_a, D2 - D2_ref>,    u_a = sum_b (G^+)_ab v_b,   G_ab = Re<v_a, v_b>,
```

is the minimal-norm linear map that is exact on span{v_a}: if |ref> and the
chi_k are eigenstates, the exact oscillation D2_ref + eps (T2_k e^{-i w t} +
h.c.) solves the linearized equations. With ccQ the e_a have zero
contractions up to the difference error (F and the exact change both carry
those of v_a) and are projected exactly; ccQ is affine and the collision
operator real-linear (it builds -D3 h as the conjugate of h D3), so the
spectral part enters after ccQ through the collision operators of the e_a,
computed whenever h2 changes. Without ccQ the e_a enter unprojected. Moving
orbitals stop with an error (their response breaks the eigenmode property).
Memory: delta3 and the e_a, (2K + 1) n^6 complex numbers (N2, n = 10,
K = 12: 400 MB; the run peaks at 0.74 GB).

`Reconstruction_Singlet_Anchor_SaveReference(kets)` writes the reference
3-RDM of kets(:, 1) and the transition RDMs of the other kets (the MCSCF
programs with `program.anchor`, `T_Be3d_01`, `T_Be3d_04`). `app/molecules/n2/groundstate/mcscf` takes
the K lowest Ritz vectors of the CI Hamiltonian in the Krylov space of Z|ref>
(`program.anchor`: `nTransitions`, `nLanczos`), the dipole-allowed excitations
for a field along the bond.

```json
"reconstruction": { "singlet": {
    "ny": { "ccQ": true },
    "anchor": { "dir": "../../groundstate/mcscf", "nTransitions": 12 }
} }
```

| Path | Type | Default | Description |
|------|------|---------|-------------|
| `reconstruction.singlet.anchor.dir` | string | "" | Directory of the reference files |
| `reconstruction.singlet.anchor.nTransitions` | int | 0 | Transitions of the spectral part (0: reference only) |

Uses `M_Utils_RdmBasisChange` (transport), `M_Utils_RdmEnforceTraceM3aab`
(projection of the e_a) and `Coeffs_FillTransitionRdm{2,3}Bt` (writer).

Tests:
- `test/simulationTests/Be3d/T_Be3d_03_TimePropagationTd2rdmAnchor`: moving
  orbitals from a rotated start (complex 1s-2s mixing, so a transposed or
  conjugated transport fails); occupations stationary to 1.4e-6 in 1 a.u., the
  split-step error of MCTDHF (unanchored NY 1.8e-4).
- `test/simulationTests/Be3d/T_Be3d_04_Td2rdmSpectralAnchor`: the NY closure at
  the reference (pinned value), stationarity, and the spectral property for
  both quadratures t = T2_k, i T2_k, with NY+ccQ and with Valdemoro without
  ccQ (relative error 2e-8).
- `app/molecules/n2/dynamics/staticOrbs`: N2 against TDCI (`app/README.md`).
