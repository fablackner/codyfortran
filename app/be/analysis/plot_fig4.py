#!/usr/bin/env python3
"""Compare MCTDHF and TD-2RDM for HHG in Be (Fig. 4 and 6 of Lackner et al.,
Phys. Rev. A 95, 033414 (2017)).

Reads <runDir>/{mctdhf,td2rdm,tdhf}_{weak4c,strong2c}/observables.dat written by the
laser programs (columns: time, E(t), A(t), <z>, -Z<z/r^3>, N_box, energy) and
plots for each pulse
  - the induced dipole moment d(t) = -<z> versus the optical cycle tau = t w/2pi,
  - the HHG spectrum from the dipole acceleration (Larmor formula, Eq. 37),
        I(w) = 2/(3c) |int a(t) e^{iwt} dt|^2 ,  a(t) = -Z <z/r^3>,
versus the harmonic order w/w0,
and optionally the time-frequency spectrum of the strong pulse (Fig. 6) from
a short-time Fourier transform with a Blackman window (alpha = 0.16).

usage: plot_fig4.py <runDir> [--out fig4.png] [--window none|hann|blackman]
                    [--spectrogram fig6.png]
Requires numpy and matplotlib.
"""
import argparse
import os

import numpy as np

OMEGA = 0.057
C_LIGHT = 137.035999
PULSES = [("weak4c", r"$I=0.5\times10^{14}$ W/cm$^2$, 4 cycles"),
          ("strong2c", r"$I=4\times10^{14}$ W/cm$^2$, 2 cycles")]
METHODS = [("mctdhf", "MCTDHF", "k", "-"), ("td2rdm", "TD-2RDM", "C3", "--"),
           ("tdhf", "TDHF", "C0", ":")]


def load(run_dir, method, pulse):
    path = os.path.join(run_dir, f"{method}_{pulse}", "observables.dat")
    if not os.path.exists(path):
        return None
    d = np.loadtxt(path, comments="#")
    return {"t": d[:, 0], "field": d[:, 1], "a_pot": d[:, 2], "z": d[:, 3],
            "acc": d[:, 4], "n": d[:, 5], "energy": d[:, 6]}


def window(n, kind):
    if kind == "hann":
        return np.hanning(n)
    if kind == "blackman":
        return np.blackman(n)
    return np.ones(n)


def hhg_spectrum(t, acc, kind="none", pad=8):
    """Larmor spectrum 2/(3c)|FT[a]|^2 on a harmonic-order axis."""
    dt = t[1] - t[0]
    a = (acc - acc[0]) * window(len(acc), kind)
    n = pad * len(a)
    spec = np.fft.rfft(a, n) * dt
    omega = 2.0 * np.pi * np.fft.rfftfreq(n, dt)
    return omega / OMEGA, 2.0 / (3.0 * C_LIGHT) * np.abs(spec) ** 2


def spectrogram(t, acc, width_cycles=0.12, alpha=0.16, max_order=60):
    """Short-time Fourier transform of a(t) with a Blackman window."""
    dt = t[1] - t[0]
    period = 2.0 * np.pi / OMEGA
    half = int(0.5 * width_cycles * period / dt)
    x = np.arange(-half, half + 1) / half  # in [-1, 1]
    blackman = (1 - alpha) / 2 + 0.5 * np.cos(np.pi * x) + alpha / 2 * np.cos(2 * np.pi * x)
    centers = np.arange(half, len(t) - half, max(1, half // 4))
    n = 8 * len(x)
    omega = 2.0 * np.pi * np.fft.rfftfreq(n, dt)
    keep = omega / OMEGA <= max_order
    rows = []
    for c in centers:
        seg = acc[c - half:c + half + 1] * blackman
        rows.append(np.abs(np.fft.rfft(seg, n) * dt)[keep] ** 2)
    return t[centers] / period, omega[keep] / OMEGA, np.array(rows).T


def main():
    ap = argparse.ArgumentParser(description="Fig. 4 / Fig. 6 of PRA 95, 033414 (2017)")
    ap.add_argument("run_dir")
    ap.add_argument("--out", default="fig4.png")
    ap.add_argument("--window", default="none", choices=["none", "hann", "blackman"])
    ap.add_argument("--max-order", type=float, default=70.0)
    ap.add_argument("--spectrogram", default=None, help="also write the Fig. 6 analogue")
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5))
    period = 2.0 * np.pi / OMEGA

    for row, (pulse, label) in enumerate(PULSES):
        ax_d, ax_s = axes[row]
        for method, name, color, style in METHODS:
            data = load(args.run_dir, method, pulse)
            if data is None:
                continue
            tau = data["t"] / period
            ax_d.plot(tau, -data["z"], color=color, ls=style, lw=1.2, label=name)
            order, spec = hhg_spectrum(data["t"], data["acc"], args.window)
            sel = (order > 0) & (order <= args.max_order)
            ax_s.semilogy(order[sel], spec[sel], color=color, ls=style, lw=1.0, label=name)
            np.savetxt(os.path.join(args.run_dir, f"{method}_{pulse}", "hhg.dat"),
                       np.column_stack([order[sel], spec[sel]]),
                       header="harmonic order   2/(3c)|a(w)|^2")
            final = data["n"][-1]
            print(f"{name:8s} {pulse:9s}: electrons left in the box {final:.6f}, "
                  f"energy drift {data['energy'][-1] - data['energy'][0]:+.3e}")
            if method == "mctdhf":
                ax_f = ax_d.twinx()
                ax_f.plot(tau, data["field"], color="0.7", lw=0.8, zorder=0)
                ax_f.set_ylabel("E(t) [a.u.]", color="0.5")
                ax_f.tick_params(axis="y", colors="0.5")
        ax_d.set_title(f"({'ac'[row]}) dipole, {label}", fontsize=10)
        ax_d.set_xlabel(r"$\tau = t\,\omega/2\pi$")
        ax_d.set_ylabel(r"$d(t) = -\langle z \rangle$ [a.u.]")
        ax_d.legend(loc="upper left", fontsize=9)
        ax_s.set_title(f"({'bd'[row]}) HHG spectrum, {label}", fontsize=10)
        ax_s.set_xlabel(r"harmonic order $\omega/\omega_0$")
        ax_s.set_ylabel(r"$I^{\rm HHG}(\omega)$ [arb. u.]")
        ax_s.legend(loc="upper right", fontsize=9)
        ax_s.set_xlim(0, args.max_order)

    fig.suptitle("Be, 800 nm: MCTDHF vs TD-2RDM (contraction-consistent Nakatsuji-Yasuda)")
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print("wrote", args.out)

    if args.spectrogram:
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
        for ax, (method, name, _, _) in zip(axes, METHODS[:2]):
            data = load(args.run_dir, method, "strong2c")
            if data is None:
                continue
            tau, order, s = spectrogram(data["t"], data["acc"])
            im = ax.pcolormesh(tau, order, np.log10(s + 1e-30), shading="auto",
                               vmin=np.log10(s.max()) - 9, vmax=np.log10(s.max()), cmap="viridis")
            ax.set_title(f"{name}, strong 2-cycle pulse", fontsize=10)
            ax.set_xlabel(r"$\tau$")
        axes[0].set_ylabel(r"harmonic order")
        fig.colorbar(im, ax=axes, label=r"$\log_{10}$ intensity")
        fig.savefig(args.spectrogram, dpi=150)
        print("wrote", args.spectrogram)


if __name__ == "__main__":
    main()
