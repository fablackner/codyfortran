#!/usr/bin/env python3
"""Generate the JSON inputs of the beryllium HHG project from one parameter set.

All programs of a pipeline must share the radial grid, the orbital space and
(for the laser runs) the pulse, absorber and propagation scheme; writing every
input from the same dictionary keeps them consistent.

    python3 make_inputs.py              # production inputs (inputs/)
    python3 make_inputs.py --test DIR   # small, fast pipeline for checking

Paper: F. Lackner, I. Brezinova, T. Sato, K. L. Ishikawa, J. Burgdoerfer,
"High-harmonic spectra from time-dependent two-particle reduced density matrix
theory", Phys. Rev. A 95, 033414 (2017), arXiv:1611.00888.
"""
import argparse
import copy
import json
import math
import os

OMEGA = 0.057  # 800 nm (paper)

# Laser pulses of Fig. 4 of the paper: F(t) = F0 cos(wt) sin^2(wt / 2Nc)
PULSES = {
    # (a,b): I = 0.5e14 W/cm^2, four cycles, no purification
    "weak4c": {"fieldStrength": 0.038, "nCycles": 4.0, "purification": False},
    # (c,d): I = 4e14 W/cm^2, two cycles, 40 purification iterations per step
    "strong2c": {"fieldStrength": 0.107, "nCycles": 2.0, "purification": True},
}

PRODUCTION = {
    # radial grid of the paper: 200 a.u., 50 FEDVR elements x 20 points
    "rmax": 200.0,
    "nElements": 50,
    "nLocals": 20,
    # angular grid (paper: Lmax = 47, all m); m = 0 orbitals only here
    "lmax": 32,
    "lmaxStageA": 2,
    # absorber: cos^(1/4) mask from r = 160 a.u. to the box edge (paper)
    "absorberOnset": 160.0,
    "absorberOrder": 4,
    # orbitals (m = 0 only): 1s 2s 2p 3s 3p 3d
    "orbitals": [(1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (3, 2)],
    "orbitalCharge": 3.0,
    # time stepping (output every outputStep, nTimeSteps split steps in between)
    "outputStep": 0.2,
    # dt = 0.0025; dt = 0.01 (20 steps) is unstable (TD-2RDM) or drifts
    # (MCTDHF) with 6 or more orbitals
    "nTimeSteps": 80,
    "postPulseTime": 0.0,
    "averagingTime": 40.0,
    # imaginary-time relaxation: final step and convergence (energy change per check)
    "dtauMin": 0.005,
    "relaxConvThresh": 1e-12,
    # TD-2RDM: number of weak spatial natural orbitals (active/weak split,
    # src/Method/A_Method.md; 0 = off)
    "nWeakOrbs": 0,
    # exact field-free one-body propagation: exp(-i h0 dt) per l, SIL only for
    # the gauge coupling, fourth-order Blanes-Moan composition inside every
    # linear split step (method.mb.orbBased.exactFieldFreeQ,
    # src/Method/A_Method.md); same accuracy as SIL at this dt, 2.5-5x faster
    # (README); not available with weak orbitals
    "exactFieldFree": True,
    # run the exact field-free transformations with cuBLAS (requires a build
    # with cmake -DCODY_USE_CUDA=ON)
    "gpu": False,
    # purify the time-averaged initial state once (P_Td2rdmTimeAverage)
    "purifyAverageQ": False,
    "pulses": PULSES,
}


def grid(p, lmax):
    return {"ylm": {"lmax": lmax, "mmax": 0, "rmin": 0.0, "rmax": p["rmax"],
                    "fedvr": {"nElements": p["nElements"], "nLocals": p["nLocals"]}}}


def hamiltonian(p, lmax):
    return {
        "grid": grid(p, lmax),
        "sysKinetic": {"ylm": {"laplacian": {"bodyMass": [1.0, 1.0], "fedvr": {}}}},
        "sysPotential": {"ylm": {"coulomb": {"charge": 4.0, "stdImpl": {}}}},
        "sysInteraction": {"ylm": {"coulomb": {"strength": 1.0, "blockEq": {}}}},
    }


def method(p, kind):
    n = len(p["orbitals"])
    opts = {}
    if kind == "td2rdm" and p["nWeakOrbs"] > 0:
        opts["nWeakOrbs"] = p["nWeakOrbs"]
    orb_based = {"nOrbs": [n, n], kind: opts}
    if p["exactFieldFree"]:
        orb_based["exactFieldFreeQ"] = True
        if p["gpu"]:
            orb_based["exactFieldFreeGpuQ"] = True
    return {"mb": {"nBodyTypes": 2, "nBodies": [2, 2], "bodyStatistics": ["f", "f"],
                   "orbBased": orb_based}}


def purification():
    # D- and Q-condition purification of the paper (40 iterations)
    return {"singlet": {"reduce": {"nSteps": 40, "threshold": -1e-5,
                                   "conserveEnergyQ": False}}}


def propagator(p):
    # second-order split step; with the exact field-free propagator the linear
    # step is the fourth-order Blanes-Moan composition of exp(-i h0 dt) and
    # the gauge SIL (src/Propagator/A_Propagator.md)
    split = {"order2": {}}
    if p["exactFieldFree"]:
        split["exactPartScheme"] = "blanesMoan6"
    return {"splitStep": split}


def ci():
    return {
        "configList": {"allActive1": {"fermionic": {"bodyTarget": 1}},
                       "allActive2": {"fermionic": {"bodyTarget": 2}}},
        "coeffs": {"generic": {}},
    }


def two_rdm():
    return {
        "twoRdm": {"fermionic": {"singlet": {}}},
        "twoRdmInit": {"load": {}},
        # contraction-consistent Nakatsuji-Yasuda reconstruction (paper)
        "reconstruction": {"singlet": {"singletNy": {"ccQ": True}}},
    }


def split_step(imag=False):
    sil = {"krylovDim": 40, "tolerance": 1e-12, "maxRestarts": 10}
    if imag:
        sil["imagTimeQ"] = True
    return {"integratorList": {"sil1": sil, "rk2": {"o4Expl": {}}}}


def stage_a(p):
    d = hamiltonian(p, p["lmaxStageA"])
    d["method"] = method(p, "mctdhx")
    d["orbs"] = {"restrictedQ": True}
    d["orbsInit"] = {"ylm": {"hydrogenLike": {
        "charge": p["orbitalCharge"],
        "n": [o[0] for o in p["orbitals"]],
        "l": [o[1] for o in p["orbitals"]],
        "m": [0] * len(p["orbitals"])}}}
    d.update(ci())
    d["coeffsInit"] = {"unary": {}}
    d.update(split_step(imag=True))
    d["program"] = {"lmaxIn": -1, "dtau": 0.05, "dtauMin": p["dtauMin"], "dtauFactor": 0.5,
                    "dtauSwitch": 1e-9, "convThresh": p["relaxConvThresh"], "nStepsCheck": 20,
                    "maxChecks": 4000, "nSubStepsNonLin": 10}
    return d


def stage_b(p):
    d = stage_a(p)
    d["grid"] = grid(p, p["lmax"])
    # orbsInit only provides placeholders: the program embeds the stage A
    # orbitals (program.lmaxIn) over them
    d["coeffsInit"] = {"load": {}}
    # stage A already ends at dtauMin; the l-mixing on the larger grid is small
    d["program"].update({"lmaxIn": p["lmaxStageA"], "dtau": p["dtauMin"], "maxChecks": 4000,
                         "nSubStepsNonLin": 2})
    return d


def time_average(p):
    d = hamiltonian(p, p["lmax"])
    d["method"] = method(p, "td2rdm")
    d["orbs"] = {"restrictedQ": True}
    d["orbsInit"] = {"load": {}}
    d.update(two_rdm())
    d.update(split_step())
    d["propagator"] = propagator(p)
    if p["purifyAverageQ"]:
        d["purification"] = purification()
    d["program"] = {"averagingTime": p["averagingTime"], "outputStep": p["outputStep"],
                    "nTimeSteps": p["nTimeSteps"]}
    return d


def laser(p, pulse, kind):
    d = hamiltonian(p, p["lmax"])
    d["sysGauge"] = {"ylm": {"velocity": {
        "fieldStrength": pulse["fieldStrength"], "omega": OMEGA,
        "nCycles": pulse["nCycles"], "cep": 0.0, "tStart": 0.0,
        "envelope": "field", "bodyMass": [1.0, 1.0], "fedvr": {}}}}
    d["method"] = method(p, kind)
    d["orbs"] = {"restrictedQ": True}
    d["orbsInit"] = {"load": {}}
    if kind == "mctdhx":
        d.update(ci())
        d["coeffsInit"] = {"load": {}}
    else:
        d.update(two_rdm())
        if pulse["purification"]:
            d["purification"] = purification()
    d.update(split_step())
    d["propagator"] = propagator(p)
    d["absorber"] = {"ylm": {"cosinus": {"onset": p["absorberOnset"],
                                         "order": p["absorberOrder"]}}}
    duration = 2.0 * math.pi * pulse["nCycles"] / OMEGA
    end = duration + p["postPulseTime"]
    end = math.ceil(end / p["outputStep"]) * p["outputStep"]
    d["program"] = {"startTime": 0.0, "endTime": round(end, 6),
                    "outputStep": p["outputStep"], "nTimeSteps": p["nTimeSteps"]}
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
    ap.add_argument("--norbs", type=int, default=None,
                    help="keep only the first NORBS orbitals (default: all)")
    ap.add_argument("--nweak-orbs", type=int, default=PRODUCTION["nWeakOrbs"],
                    help="TD-2RDM weak natural orbitals, time average and laser "
                         "inputs (default: 0, off)")
    ap.add_argument("--no-exact-field-free", action="store_true",
                    help="propagate the whole one-body part with SIL instead of "
                         "the exact field-free propagator")
    ap.add_argument("--gpu", action="store_true",
                    help="exact field-free transformations on the GPU (cuBLAS; "
                         "build with cmake -DCODY_USE_CUDA=ON)")
    args = ap.parse_args()

    p = copy.deepcopy(PRODUCTION)
    if args.norbs:
        p["orbitals"] = p["orbitals"][:args.norbs]
    p["nWeakOrbs"] = args.nweak_orbs
    if args.no_exact_field_free:
        p["exactFieldFree"] = False
    if args.gpu:
        p["gpu"] = True
    if p["nWeakOrbs"] > 0 and p["exactFieldFree"]:
        # TD-2RDM weak orbitals replace h by H_L; keep both methods on one scheme
        print("note: --nweak-orbs switches off the exact field-free propagator")
        p["exactFieldFree"] = False
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "inputs")
    if args.test:
        out = args.test
        p.update({"rmax": 40.0, "nElements": 10, "nLocals": 12, "lmax": 6,
                  "absorberOnset": 30.0, "dtauMin": 0.025, "relaxConvThresh": 1e-9,
                  "averagingTime": 4.0})
        for pulse in p["pulses"].values():
            pulse["nCycles"] = 0.5

    write(os.path.join(out, "stageA.json"), stage_a(p))
    write(os.path.join(out, "stageB.json"), stage_b(p))
    write(os.path.join(out, "timeAverage.json"), time_average(p))
    for name, pulse in p["pulses"].items():
        write(os.path.join(out, f"mctdhf_{name}.json"), laser(p, pulse, "mctdhx"))
        write(os.path.join(out, f"td2rdm_{name}.json"), laser(p, pulse, "td2rdm"))

    # TDHF = MCTDHF with as many orbitals as electrons per spin (1s, 2s)
    hf = copy.deepcopy(p)
    hf["orbitals"] = p["orbitals"][:2]
    write(os.path.join(out, "stageA_hf.json"), stage_a(hf))
    write(os.path.join(out, "stageB_hf.json"), stage_b(hf))
    for name, pulse in p["pulses"].items():
        write(os.path.join(out, f"tdhf_{name}.json"), laser(hf, pulse, "mctdhx"))


if __name__ == "__main__":
    main()
