"""Inputs of P_KickTd2rdmFem for an active space of ethylene.

usage: python3 make_inputs.py NAME NVIRT [key=value ...]
  NAME   name of the run (directory runs/NAME/ with input and Hamiltonian)
  NVIRT  number of active virtual orbitals
options (DEFAULTS): select=mp2|cis, ints=FILE, twoRdm=singlet|singletTr0
  (2-RDM storage: ab block with the NY closure and ccQ, or explicit 1-RDM +
  trace-free ab block), switchTime (method.mb.orbBased.adiabaticSwitching),
  kick, polarization, endTime, outputStep, nTimeSteps (program block of
  P_KickTd2rdmFem)

Basis (ints, P_IntegralsFem): the 8 occupied RHF orbitals of
../../groundstate/scf and 39 raw STOs (C 2s 2p 3d, H 1s, 3s 3p 3d on the ghost
center at the C=C midpoint) with their overlap. The STOs are projected on the
complement of the occupied space and canonically orthogonalized. select=mp2
rotates the virtual space to MP2 natural orbitals (C 1s frozen) and takes the
NVIRT most occupied ones; select=cis puts the particle natural transition
orbitals of the bright CIS states (along the kick, state by state from the
lowest) first and fills up with MP2 natural orbitals of the rest of the
virtual space. The C 1s orbitals are folded into the one-body operator.
runs/NAME/hamiltonian.bin: int32 n, real64 core energy, complex128 h1eff(n, n),
h2(n, n, n, n) = <ij|kl>/2, dip(n, n, 3); runs/NAME/hamiltonian.npz (real
arrays): h1 = h1eff, g = <ij|kl> (not halved), dip, ecore, enuc. The input uses a coarse grid with nominal Slater orbitals: they only
give the method its shape (h1 and h2 are replaced).
"""
import json
import os
import sys

import numpy as np

SCF = os.path.abspath('../../groundstate/scf')
N_OCC, CORE = 8, [0, 1]
DEFAULTS = dict(select='mp2', ints='../../analysis/integrals/intsDiffuse.bin', twoRdm='singlet', switchTime='40',
                kick='0.01', polarization='3', endTime='80', outputStep='0.2', nTimeSteps='1')


def read_ints(fn):
    f = open(fn, 'rb')
    n = int(np.fromfile(f, np.int32, 1)[0])
    enuc = float(np.fromfile(f, np.float64, 1)[0])
    rd = lambda m: np.fromfile(f, np.complex128, m)
    h1 = rd(n * n).reshape(n, n, order='F')
    g = 2.0 * rd(n ** 4).reshape((n,) * 4, order='F')        # <ij|kl>
    dip = rd(3 * n * n).reshape(n, n, 3, order='F')
    s = rd(n * n).reshape(n, n, order='F')
    return n, enuc, h1.real, g.real, dip.real, s.real


def transform(h1, g, dip, c):
    h1 = c.T @ h1 @ c
    g = np.einsum('abcd,ai,bj,ck,dl->ijkl', g, c, c, c, c, optimize=True)
    dip = np.einsum('abx,ai,bj->ijx', dip, c, c)
    return h1, g, dip


def orthonormal_basis(s):
    """C with C^T S C = 1: occupied functions as they are, the rest projected
    and canonically orthogonalized (eigenvalues of the overlap above 1e-6)."""
    n = s.shape[0]
    cv = np.eye(n)[:, N_OCC:]
    cv = cv - np.eye(n)[:, :N_OCC] @ s[:N_OCC, N_OCC:]       # project out the occupied orbitals
    sv = cv.T @ s @ cv
    w, u = np.linalg.eigh(sv)
    keep = w > 1e-6 * w.max()
    c = np.hstack([np.eye(n)[:, :N_OCC], cv @ u[:, keep] / np.sqrt(w[keep])])
    return c, int(keep.sum())


def fock(h1, g, occ):
    """Fock matrix of the doubly occupied orbitals occ."""
    return h1 + sum(2.0 * g[:, j, :, j] - g[:, j, j, :] for j in occ)


def mp2_natural_virtuals(h1, g):
    """MP2 (C 1s frozen) in the canonical virtuals: the virtual natural orbitals
    (columns in the virtual functions, descending occupation), their
    occupations, the energy and the virtual density matrix."""
    f = fock(h1, g, range(N_OCC))
    ev, uv = np.linalg.eigh(f[N_OCC:, N_OCC:])
    n = h1.shape[0]
    c = np.eye(n); c[N_OCC:, N_OCC:] = uv
    h1c, gc, _ = transform(h1, g, np.zeros((n, n, 3)), c)
    eps = np.concatenate([np.diag(fock(h1c, gc, range(N_OCC)))[:N_OCC], ev])
    occ = [i for i in range(N_OCC) if i not in CORE]
    vir = list(range(N_OCC, n))
    k = gc[np.ix_(occ, occ, vir, vir)]                        # <ij|ab> = (ia|jb)
    d = eps[occ][:, None, None, None] + eps[occ][None, :, None, None] - eps[vir][None, None, :, None] - eps[vir][None, None, None, :]
    t = k / d
    e2 = np.sum(t * (2.0 * k - k.transpose(0, 1, 3, 2)))
    pvv = 2.0 * np.einsum('ijac,ijbc->ab', t, 2.0 * t - t.transpose(0, 1, 3, 2))
    w, u = np.linalg.eigh(pvv)
    order = np.argsort(-w)
    return uv @ u[:, order], w[order], e2, uv @ pvv @ uv.T


def cis_particle_orbitals(h1, g, nvirt, polarization, dip, emax=0.75, fmin=0.02, wmin=0.1):
    """Particle natural transition orbitals of the bright singlet CIS states
    (C 1s frozen) polarized along `polarization`, below emax: columns in the
    virtual space (orthonormal, at most nvirt)."""
    f = fock(h1, g, range(N_OCC))
    occ = [i for i in range(N_OCC) if i not in CORE]
    vir = list(range(N_OCC, h1.shape[0]))
    no, nv = len(occ), len(vir)
    ei = np.diag(f)[occ]
    fvv = f[np.ix_(vir, vir)]
    ev, uv = np.linalg.eigh(fvv)                                  # canonical virtuals
    c = np.eye(h1.shape[0]); c[N_OCC:, N_OCC:] = uv
    _, gc, dc = transform(h1, g, dip, c)
    a = (2.0 * gc[np.ix_(occ, occ, vir, vir)].transpose(0, 2, 1, 3)
         - gc[np.ix_(occ, vir, occ, vir)]).reshape(no * nv, no * nv)
    a += np.diag((ev[None, :] - ei[:, None]).ravel())
    w, t = np.linalg.eigh(a)
    mu = np.sqrt(2.0) * dc[np.ix_(occ, vir, [polarization - 1])][:, :, 0].ravel() @ t
    osc = 2.0 / 3.0 * w * mu ** 2
    vecs = []
    print('CIS bright states (polarization %d): energy (Ha), oscillator strength' % polarization)
    for k in np.where((w < emax) & (osc > fmin))[0]:
        u, sv, vt = np.linalg.svd(t[:, k].reshape(no, nv))
        print('  %.4f  %.4f  NTO weights %s' % (w[k], osc[k], np.round(sv[:3] ** 2, 3)))
        vecs += [vt[m] for m in range(len(sv)) if sv[m] ** 2 > wmin]
    # Gram-Schmidt in the order of the states (lowest bright state first, its
    # dominant NTO first): a truncated selection keeps the low-lying states
    p = []
    for x in vecs:
        for y in p:
            x = x - (y @ x) * y
        if np.linalg.norm(x) > 0.1:
            p.append(x / np.linalg.norm(x))
    if not p:
        return np.zeros((nv, 0))
    return uv @ np.array(p[:nvirt]).T


def main():
    name, nvirt = sys.argv[1], int(sys.argv[2])
    given = dict(a.split('=', 1) for a in sys.argv[3:])
    if set(given) - set(DEFAULTS):
        sys.exit('unknown options: ' + ', '.join(sorted(set(given) - set(DEFAULTS))))
    opt = {**DEFAULTS, **given}
    if opt['select'] not in ('mp2', 'cis') or opt['twoRdm'] not in ('singlet', 'singletTr0'):
        sys.exit('select must be mp2 or cis, twoRdm singlet or singletTr0')
    n, enuc, h1, g, dip, s = read_ints(opt['ints'])
    c, nv = orthonormal_basis(s)
    h1, g, dip = transform(h1, g, dip, c)
    no, wno, e2, pvv = mp2_natural_virtuals(h1, g)
    nb = N_OCC + nv
    if opt['select'] == 'cis':
        # particle orbitals of the bright CIS states first, the rest MP2
        # natural orbitals of the complement
        p = cis_particle_orbitals(h1, g, nvirt, int(opt['polarization']), dip)
        q = np.eye(nv) - p @ p.T
        w, u = np.linalg.eigh(q @ pvv @ q)
        no = np.hstack([p, u[:, np.argsort(-w)][:, :nv - p.shape[1]]])
        print('%d CIS particle orbitals in the active space' % p.shape[1])
    # rotate the virtual space (selected orbitals first)
    rot = np.eye(nb); rot[N_OCC:, N_OCC:] = no
    h1, g, dip = transform(h1, g, dip, rot)
    print('basis: %d virtual orbitals after orthogonalization, E(MP2, C 1s frozen) = %.6f' % (nv, e2))
    print('MP2 natural occupations of the virtuals:', np.round(wno[:max(nvirt, 12)], 5))

    hc = fock(h1, g, CORE)
    ecore = sum(h1[k, k] + hc[k, k] for k in CORE)
    active = [i for i in range(N_OCC) if i not in CORE] + list(range(N_OCC, N_OCC + nvirt))
    na = len(active)
    ix = np.ix_(active, active)

    run = os.path.join('runs', name)
    os.makedirs(run, exist_ok=True)
    with open(run + '/hamiltonian.bin', 'wb') as f:
        np.array([na], np.int32).tofile(f)
        np.array([ecore], np.float64).tofile(f)
        hc[ix].astype(np.complex128).ravel(order='F').tofile(f)
        (0.5 * g[np.ix_(active, active, active, active)]).astype(np.complex128).ravel(order='F').tofile(f)
        dip[np.ix_(active, active, [0, 1, 2])].astype(np.complex128).ravel(order='F').tofile(f)
    np.savez(run + '/hamiltonian.npz', h1=hc[ix], g=g[np.ix_(active, active, active, active)],
             dip=dip[np.ix_(active, active, [0, 1, 2])], ecore=ecore, enuc=enuc)

    # the method needs orbitals of some grid only for its shape: a coarse grid
    # with na Slater functions (h1 and h2 of the setup are replaced)
    scf = json.load(open(SCF + '/P_GroundScfFem.json'))
    d = {k: scf[k] for k in ['grid', 'sysKinetic', 'sysPotential', 'sysInteraction']}
    gf = d['grid']['fem']
    for a in 'xyz':
        gf[a + 'min'], gf[a + 'max'] = -6.0, 6.0
    gf.update(nLocals=3, hMax=3.0, growth=2.0, multipoleLmax=2)
    gf['centers']['hMin'] = [1.0] * len(gf['centers']['hMin'])
    d['method'] = {'mb': {'nBodyTypes': 2, 'nBodies': [N_OCC - len(CORE)] * 2, 'bodyStatistics': ['f', 'f'],
                          'orbBased': {'nOrbs': [na, na], 'td2rdmStaticOrbs': {},
                                       'adiabaticSwitching': {'switchTime': float(opt['switchTime'])}}}}
    d['orbs'] = {'restrictedQ': True}
    d['orbsInit'] = {'fem': {'sto': {'orbital': list(range(1, na + 1)), 'center': [1] * na, 'n': [1] * na,
                                     'l': [0] * na, 'm': [0] * na,
                                     'zeta': [0.4 + 0.15 * k for k in range(na)], 'coeff': [1.0] * na}}}
    d['twoRdm'] = {'fermionic': {opt['twoRdm']: {}}}
    d['twoRdmInit'] = {'excited': {}}
    if opt['twoRdm'] == 'singlet':
        d['reconstruction'] = {'singlet': {'ny': {'ccQ': True}}}
    else:                                                      # explicit 1-RDM, trace-free 2-RDM
        d['reconstruction'] = {'singletTr0': {'ny': {}}}
    d['integratorList'] = {'rk': {'o4Expl': {}}}
    d['propagator'] = {'single': {}}
    d['program'] = {'hamiltonianFile': 'hamiltonian.bin', 'kick': float(opt['kick']), 'polarization': int(opt['polarization']),
                    'endTime': float(opt['endTime']), 'outputStep': float(opt['outputStep']),
                    'nTimeSteps': int(opt['nTimeSteps']), 'outputName': 'td2rdm'}
    json.dump(d, open(run + '/P_KickTd2rdmFem.json', 'w'), indent=4)
    print('%s: %d active orbitals (6 occupied valence + %d virtual), core energy %.10f' % (name, na, nvirt, ecore))


if __name__ == '__main__':
    main()
