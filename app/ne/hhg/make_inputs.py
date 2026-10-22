#!/usr/bin/env python3
"""Generate the JSON inputs of the neon strong-field MCTDHF runs.

Ne (10 electrons) with 9 restricted spatial orbitals 1s 2s 2p 3s 3p (all m of
the p shells, grid.ylm.mmax = 1) in a z-polarized 800 nm few-cycle pulse. The
programs are the beryllium ones (app/be): relax/P_ImagTimeMctdhf for the
ground state, laser/P_LaserMctdhf for the pulse.

    python3 make_inputs.py              # production inputs (inputs/)
    python3 make_inputs.py --test DIR   # small, fast pipeline for checking
"""
import argparse
import copy
import json
import math
import os

OMEGA = 0.057  # 800 nm

# F(t) = F0 cos(wt) sin^2(wt / 2Nc), F0 = sqrt(I / 3.51e16 W/cm^2)
PULSES = {
    # I = 5e14 W/cm^2 (below the barrier-suppression intensity 8.6e14 of Ne),
    # quiver radius F0/w^2 = 37 a.u., cutoff Ip + 3.17 Up = 4.3 a.u. (harmonic 75)
    "strong2c": {"fieldStrength": 0.1194, "nCycles": 2.0},
}

PRODUCTION = {
    "charge": 10.0,
    # radial grid: FEDVR, 4 a.u. elements with 20 points (hydrogenic Ne 1s
    # eigenvalue exact to 5e-5; 2 a.u. x 14 points: 1e-3 with more points),
    # absorber from 120 a.u. (3.3 quiver radii)
    "rmax": 160.0,
    "nElements": 40,
    "nLocals": 20,
    "absorberOnset": 120.0,
    "absorberOrder": 4,
    # angular grid: l <= lmax, |m| <= 1 (m of every orbital is conserved)
    "lmax": 40,
    # the ground-state orbitals keep a pure l, so lmax = 1 is exact for stage A
    "lmaxStageA": 1,
    # orbitals (n, l, m)
    "orbitals": [(1, 0, 0), (2, 0, 0), (2, 1, -1), (2, 1, 0), (2, 1, 1),
                 (3, 0, 0), (3, 1, -1), (3, 1, 0), (3, 1, 1)],
    "orbitalCharge": 8.0,
    # time stepping: dt = outputStep / nTimeSteps = 0.0025
    "outputStep": 0.2,
    "nTimeSteps": 80,
    # checkpoint (orbitals, CI vector) every 800 steps
    "saveStep": 2.0,
    "postPulseTime": 0.0,
    # imaginary-time relaxation (the Strang fixed point is biased by O(dtau^2):
    # ~6 mHa in the HF energy at dtau = 0.005 for Ne)
    "dtauMin": 0.001,
    "relaxConvThresh": 1e-11,
    # exact field-free h0 with a spectral cut and the Strang composition
    # (fastest accurate MCTDHF setup of the Be benchmark, app/be/README.md)
    "spectralCut": 1000.0,
    "pulses": PULSES,
}


def grid(p, lmax):
    return {"ylm": {"lmax": lmax, "mmax": 1, "rmin": 0.0, "rmax": p["rmax"],
                    "fedvr": {"nElements": p["nElements"], "nLocals": p["nLocals"]}}}


def hamiltonian(p, lmax):
    return {
        "grid": grid(p, lmax),
        "sysKinetic": {"ylm": {"laplacian": {"bodyMass": [1.0, 1.0], "fedvr": {}}}},
        "sysPotential": {"ylm": {"coulomb": {"charge": p["charge"], "stdImpl": {}}}},
        "sysInteraction": {"ylm": {"coulomb": {"strength": 1.0, "blockEq": {}}}},
    }


def method(p):
    n = len(p["orbitals"])
    nElectrons = int(p["charge"])
    # declared m of every spin orbital (both spins): binds the definite-m node
    # quadrature and the m-selective exact part, and the orbitals are
    # projected on it once when they are initialized (src/Method/A_Method.md)
    m_of_orbs = [o[2] for o in p["orbitals"]] * 2
    orb_based = {"nOrbs": [n, n], "mOfOrbs": m_of_orbs, "mctdhx": {}, "exactFieldFreeQ": True}
    if p["spectralCut"] > 0.0:
        orb_based["exactFieldFreeCut"] = p["spectralCut"]
    return {"mb": {"nBodyTypes": 2, "nBodies": [nElectrons // 2, nElectrons // 2],
                   "bodyStatistics": ["f", "f"], "orbBased": orb_based}}


def propagator(p):
    scheme = "strang" if p["spectralCut"] > 0.0 else "blanesMoan6Swap"
    return {"splitStep": {"order2": {}, "exactPartScheme": scheme}}


def ci():
    # CI restricted to the total-M sector of the closed-shell ground state
    # (M = 0; requires mOfOrbs, src/Coeffs/A_Coeffs.md)
    return {
        "configList": {"allActive1": {"fermionic": {"bodyTarget": 1}},
                       "allActive2": {"fermionic": {"bodyTarget": 2}}},
        "coeffs": {"generic": {"totalM": 0}},
    }


def split_step(imag=False):
    sil = {"krylovDim": 40, "tolerance": 1e-10, "maxRestarts": 10}
    if imag:
        sil["imagTimeQ"] = True
        sil["tolerance"] = 1e-12
    return {"integratorList": {"sil1": sil, "rk2": {"o4Expl": {}}}}


def stage_a(p):
    d = hamiltonian(p, p["lmaxStageA"])
    d["method"] = method(p)
    d["orbs"] = {"restrictedQ": True}
    d["orbsInit"] = {"ylm": {"hydrogenLike": {
        "charge": p["orbitalCharge"],
        "n": [o[0] for o in p["orbitals"]],
        "l": [o[1] for o in p["orbitals"]],
        "m": [o[2] for o in p["orbitals"]]}}}
    d.update(ci())
    d["coeffsInit"] = {"unary": {}}
    d.update(split_step(imag=True))
    d["program"] = {"lmaxIn": -1, "dtau": 0.02, "dtauMin": p["dtauMin"], "dtauFactor": 0.5,
                    "dtauSwitch": 1e-8, "convThresh": p["relaxConvThresh"], "nStepsCheck": 50,
                    "maxChecks": 4000, "dtauNonLinMax": 0.002}
    return d


def stage_b(p):
    # embeds the stage A state into the propagation grid; the orbitals keep a
    # pure l, so one check only confirms the energy
    d = stage_a(p)
    d["grid"] = grid(p, p["lmax"])
    d["coeffsInit"] = {"load": {}}
    d["program"].update({"lmaxIn": p["lmaxStageA"], "dtau": p["dtauMin"], "maxChecks": 1,
                         "nStepsCheck": 10})
    return d


def laser(p, pulse):
    d = hamiltonian(p, p["lmax"])
    d["sysGauge"] = {"ylm": {"velocity": {
        "fieldStrength": pulse["fieldStrength"], "omega": OMEGA,
        "nCycles": pulse["nCycles"], "cep": 0.0, "tStart": 0.0,
        "envelope": "field", "bodyMass": [1.0, 1.0], "fedvr": {}}}}
    d["method"] = method(p)
    d["orbs"] = {"restrictedQ": True}
    d["orbsInit"] = {"load": {}}
    d.update(ci())
    d["coeffsInit"] = {"load": {}}
    d.update(split_step())
    d["propagator"] = propagator(p)
    d["absorber"] = {"ylm": {"cosinus": {"onset": p["absorberOnset"],
                                         "order": p["absorberOrder"]}}}
    duration = 2.0 * math.pi * pulse["nCycles"] / OMEGA
    end = math.ceil((duration + p["postPulseTime"]) / p["outputStep"]) * p["outputStep"]
    d["program"] = {"startTime": 0.0, "endTime": round(end, 6),
                    "outputStep": p["outputStep"], "nTimeSteps": p["nTimeSteps"],
                    "saveStep": p["saveStep"]}
    return d


def write(path, d):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        json.dump(d, f, indent=4)
        f.write("\n")
    print("wrote", path)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--test", metavar="DIR",
                    help="write a small, fast configuration into DIR instead")
    ap.add_argument("--out", metavar="DIR", help="output directory (default: inputs/)")
    ap.add_argument("--lmax", type=int, default=PRODUCTION["lmax"])
    ap.add_argument("--norbs", type=int, default=None,
                    help="keep only the first NORBS orbitals (5: Hartree-Fock)")
    args = ap.parse_args()

    p = copy.deepcopy(PRODUCTION)
    p["lmax"] = args.lmax
    if args.norbs:
        p["orbitals"] = p["orbitals"][:args.norbs]
    out = args.out or os.path.join(os.path.dirname(os.path.abspath(__file__)), "inputs")
    if args.test:
        out = args.test
        p.update({"rmax": 40.0, "nElements": 10, "nLocals": 20, "lmax": 6,
                  "absorberOnset": 30.0, "dtauMin": 0.005, "relaxConvThresh": 1e-8})
        for pulse in p["pulses"].values():
            pulse["nCycles"] = 0.5

    write(os.path.join(out, "stageA.json"), stage_a(p))
    write(os.path.join(out, "stageB.json"), stage_b(p))
    for name, pulse in p["pulses"].items():
        write(os.path.join(out, f"mctdhf_{name}.json"), laser(p, pulse))


if __name__ == "__main__":
    main()
