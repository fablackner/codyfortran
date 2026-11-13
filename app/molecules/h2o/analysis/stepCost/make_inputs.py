"""Writes the inputs of P_StepCostFem: H2O on the mesh of ../../groundstate/scf,
MCTDHF (full CI) and TD-2RDM (NY, ccQ) for nOrbs orbitals (STOs on O and H).

usage: python3 make_inputs.py 6 8 10 12 14 16 (writes the inputs to runs/)
"""
import copy
import json
import os
import sys

# (center, n, l, m, zeta) of one STO per orbital, in the order they are added
STOS = [(1, 1, 0, 0, 7.66), (1, 2, 0, 0, 2.25), (1, 2, 1, 1, 2.23), (1, 2, 1, 0, 2.23), (1, 2, 1, -1, 2.23),
        (2, 1, 0, 0, 1.0), (3, 1, 0, 0, 1.0), (1, 3, 0, 0, 1.0), (1, 3, 1, 1, 1.0), (1, 3, 1, 0, 1.0),
        (1, 3, 1, -1, 1.0), (1, 3, 2, 0, 1.2), (1, 3, 2, 1, 1.2), (1, 3, 2, -1, 1.2), (1, 3, 2, 2, 1.2),
        (1, 3, 2, -2, 1.2)]

scf = json.load(open('../../groundstate/scf/P_GroundScfFem.json'))
base = {k: scf[k] for k in ['grid', 'sysKinetic', 'sysPotential', 'sysInteraction']}
base['orbs'] = {'restrictedQ': True}
base['program'] = {'nRepeat': 2}
os.makedirs('runs', exist_ok=True)

for n in map(int, sys.argv[1:]):
    sto = STOS[:n]
    init = {'fem': {'sto': {'orbital': list(range(1, n + 1)), 'center': [s[0] for s in sto],
                            'n': [s[1] for s in sto], 'l': [s[2] for s in sto], 'm': [s[3] for s in sto],
                            'zeta': [s[4] for s in sto], 'coeff': [1.0] * n}}}
    mb = {'nBodyTypes': 2, 'nBodies': [5, 5], 'bodyStatistics': ['f', 'f']}

    d = copy.deepcopy(base)
    d['orbsInit'] = init
    d['method'] = {'mb': dict(mb, orbBased={'nOrbs': [n, n], 'mctdhx': {}})}
    d['configList'] = {'allActive1': {'fermionic': {'bodyTarget': 1}},
                       'allActive2': {'fermionic': {'bodyTarget': 2}}}
    d['coeffs'] = {'generic': {}}
    d['coeffsInit'] = {'unary': {}}
    json.dump(d, open('runs/mctdhf_%02d.json' % n, 'w'), indent=1)

    d = copy.deepcopy(base)
    d['orbsInit'] = init
    d['method'] = {'mb': dict(mb, orbBased={'nOrbs': [n, n], 'td2rdm': {}})}
    d['twoRdm'] = {'fermionic': {'singlet': {}}}
    d['twoRdmInit'] = {'excited': {}}
    d['reconstruction'] = {'singlet': {'singletNy': {'ccQ': True}}}
    json.dump(d, open('runs/td2rdm_%02d.json' % n, 'w'), indent=1)
