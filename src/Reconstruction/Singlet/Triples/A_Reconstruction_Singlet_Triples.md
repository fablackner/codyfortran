# Connected Triples of the Singlet Closure (`reconstruction.singlet.triples`)

An optional extension of the singlet family (`src/Reconstruction/A_Reconstruction.md`).
`Reconstruction_Singlet_Fabricate` calls `Reconstruction_Singlet_Triples_Fabricate`
when a `triples` block exists; the extension wraps the closure
`Reconstruction_Singlet_FillRdm3aab` bound before it and binds
`Reconstruction_SetHamiltonian`. Nothing outside this folder refers to it.

```json
"reconstruction": { "singlet": {
    "ny": { "ccQ": true },
    "triples": { "moellerPlesset": {} }
} }
```

## What it adds

In the natural orbitals of the 1-RDM, split into the N/2 most occupied
("occupied" o) and the others ("unoccupied" v), the 3-cumulant has one block
that annihilates three occupied and creates three unoccupied orbitals,
Δ3(ooo|vvv), and its adjoint. For a weakly correlated state it is the
connected triple excitation, of second order in the interaction. The
Nakatsuji-Yasuda products of two 2-cumulants are of higher order there. At
the N2 MCSCF state (10 orbitals) the two blocks hold 82% of the squared error
of NY: the exact block has norm 2.86e-2, NY gives 7.5e-4.

The extension builds the block from the 2-cumulant Δ2, the interaction and
the orbital energies, as CCSD(T) builds the triples from the doubles:

```
S(ijk|abc) = P(k/ij) P(a/bc)  sum_e Δ2(ij|ae) <ek||bc>      particle line
           - P(i/jk) P(c/ab)  sum_m Δ2(im|ab) <jk||mc>      hole line
D(ijk|abc) = ε_i + ε_j + ε_k − ε_a − ε_b − ε_c
Δ3(ijk|abc) = S / D                                          moellerPlesset
```

in spin orbitals, ket-first (Δ3(ijk|abc) is the cumulant of
<a+_a a+_b a+_c a_k a_j a_i>), P(x/yz) f = f − f(x↔y) − f(x↔z).

The orbitals are semicanonical: the occupied and the unoccupied natural
orbitals, each set rotated to diagonalize its block of the Fock operator
f = h1 + J − K of the 1-RDM, ε its diagonal. In the natural orbitals
themselves f has large off-diagonal elements (0.69 Ha between the two N 1s
orbitals of N2) and the partition into D and the rest is poor. Rotations
inside a degenerate level do not change the result.

The added term has no onefold contraction (a contraction would pair an
occupied with an unoccupied index), so ccQ leaves it unchanged and the
contractions of D3 stay those of D2.

## Implementation

Only the up-up-down block is computed, in spatial orbitals. For a singlet
the up-down 2-cumulant c obeys c(p,q,r,s) = c(q,p,s,r), and the 18 terms of
S for two up and one down spin on each side become 12 entries each of two
matrix products (`FillSource`):

```
pLine((p,q,r),(x,y,z)) = sum_e c(p,q,r,e) <ex|yz>     e unoccupied
hLine((p,r,s),(x,y,z)) = sum_m c(p,m,r,s) <xy|mz>     m occupied
```

The block is then embedded in all orbitals of the 2-RDM
(`RdmBasisChange_EmbedBlock`, arrays growing from n_occ³ n_virt³ to n⁶) and
its adjoint added. Cost: two transformations of n⁴ tensors (n⁵), the two
products (n_occ³ n_virt⁴ + n_occ⁴ n_virt³) and passes over n⁶; no step of
order n⁷. With 4 threads the extension takes 0.4 to 0.7 of the time of the
NY closure itself (10 orbitals 0.007 s, 12 orbitals 0.034 s, 14 orbitals
0.09 s).

```
Triples/
  M_Reconstruction_Singlet_Triples.f90   contract: Fabricate, FillCumul3, FillAmplitudes pointer
  S_Reconstruction_Singlet_Triples.f90   driver: orbitals, source, embedding, closure wrapper
  MoellerPlesset/                        the variant, binds FillAmplitudes
```

A variant receives the source and the orbital energies and returns the
amplitudes, both as x(i,j,k,a,b,c) with i, j up and k down occupied, a, b up
and c down unoccupied. A new variant is an M_/S_ pair in its own folder
whose fabricate binds `Reconstruction_Singlet_Triples_FillAmplitudes`, and a
branch in `Reconstruction_Singlet_Triples_Fabricate`. A change of the source
(for example a screened interaction at the vertex) belongs in `FillSource`.

With `anchor` in the same input the anchor wraps this closure (delta3 is
then the error of NY plus triples). Its setup evaluates the closure before
the first time derivative, so a program must call
`Reconstruction_SetHamiltonian` before `Reconstruction_Setup`.

## The Hamiltonian

The closure depends on the Hamiltonian. The method passes it through the
generic hook `Reconstruction_SetHamiltonian(h1, h2)` (M_Reconstruction,
default no-op) before every `Reconstruction_AddCollisionOp`; h2 carries every
factor it has in the Hamiltonian. `td2rdmStaticOrbs` does this, with λ h2
during adiabatic switching, so source and orbital energies are those of
H(λ). Methods with moving orbitals do not call the hook yet: the closure
then stops with an error at the first call. To support them, call the hook
with the one-body matrix and h2 in the current orbitals at every place that
calls `Reconstruction_AddCollisionOp` (the split-step parts included).

## What it does to a ground state

Energy error of the adiabatically switched ground state against the exact
energy of the same orbitals, in mHa, switch time 40 / 80 a.u.
(`app/README.md`):

| system | E_corr | NY + ccQ | + triples |
|---|---|---|---|
| C2H4, 10 orbitals | −22.7 | −0.059 / −0.063 | −0.049 / −0.052 |
| C2H4, 12 orbitals | −28.0 | −0.095 / −0.105 | −0.047 / −0.052 |
| Ne, 9 orbitals | −136.0 | −0.009 / −0.009 | +0.016 / +0.016 |
| N2, 10 orbitals | −150.6 | −0.061 / −0.196 | +0.141 / +0.211 |
| Be, 5 orbitals | −44.5 | −0.102 / −0.155 | no triples (parity) |

The term raises the energy in every system. It helps where NY is too low by
more than that (C2H4, where the smallest eigenvalues of D and Q also move a
quarter closer to zero at 40 a.u.) and overshoots in Ne and N2. It is not a
systematic improvement.

N2 is the worst case. At its exact 2-RDM the term reduces |D3 − R(D2)| from
4.5e-2 to 1.9e-2 and the error of the collision term to 45% (block oo|vv
1.34e-2 to 4.8e-3). In a propagation its part on the block ov|vv drives a
coherence between the nearly degenerate pairs σu σ* and π*π* (energy
difference 0.1 Ha). In the exact equations the Nakatsuji-Yasuda term t2*·S
cancels this drive up to the small energy difference; with the 2-cumulant of
N2 the cancellation fails by a factor of ten. The coherence becomes complex
and drains the π*π* pair population at about 4e-4 per a.u., for every switch
time and also from the exact state: min eigenvalue of D −4.5e-3 after a
switch of 40 a.u., −8.4e-3 after 80 a.u.

Tried and not kept (worse or no better in every test, numbers in
`app/README.md`):

- Epstein-Nesbet denominators D − K_diag and the resummation (D − K)⁻¹ S of
  the ladder and ring interactions among the six lines (GMRES). Both made
  the N2 state worse than S / D; they need the amplitudes of all spin
  components.
- Keeping the collision operator of the block only on the doubles blocks
  oo|vv and vv|oo, where it is −¼ S†T, a Hermitian self-energy of the
  2-cumulant that lowers the energy. D stays positive in N2, but the energy
  is no longer conserved (4e-5 Ha in 20 a.u.) and the error grows in C2H4
  and N2 (NY alone is already below the exact energy).

A correction that helps throughout would have to treat the
number-conserving blocks (oo|oo, ov|ov) to third order together with the
triples.

## Tests

- `test/componentTests/Reconstruction/Singlet/Triples/T_Reconstruction_Singlet_Triples`:
  the cumulant (`Reconstruction_Singlet_Triples_FillCumul3`) on a complex
  model with 5 orbitals: antisymmetry, Hermiticity, vanishing contractions,
  covariance under a complex orbital rotation, and its norm, pinned from the
  spin-orbital formula above.
- `app/molecules/n2/groundstate/td2rdm`: N2 ground state by adiabatic
  switching against the CI energy of the same orbitals (`app/README.md`).
