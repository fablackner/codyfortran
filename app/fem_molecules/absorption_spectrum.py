"""Linear-response absorption spectrum of a delta-kick run (P_FemKickTdhf).

The kick phi -> exp(i kappa r_a) phi of every orbital at t = 0 is the
response to the field E(t) = -kappa e_a delta(t). The induced dipole
d(t) = <r_a>(t) - <r_a>(0) (position expectation summed over the electrons,
column a of dipole.dat) gives the dynamic polarizability

    alpha(omega) = (1/kappa) int_0^T d(t) w(t) e^{i omega t} dt,   w(t) = e^{-t/tau}

(trapezoidal rule on the output times; the damping turns every line into a
Lorentzian of half width 1/tau), the photo-absorption cross section (atomic
units; 1 a.u. of area = 28.0029 Mb)

    sigma(omega) = 4 pi omega Im alpha(omega) / c,   c = 137.036,

and the dipole strength function S(omega) = (2 omega / pi) Im alpha(omega).
For d(t) = kappa sum_n (f_n / omega_n) sin(omega_n t), alpha is
sum_n f_n / (omega_n^2 - (omega + i/tau)^2) and the integral of S over all
omega is sum_n f_n = N, the number of electrons, for any tau (Thomas-Reiche-
Kuhn sum rule). The script checks it twice:
  - initial slope: d(t) = N kappa t + O(t^3), from the first two samples
    (Richardson extrapolation of d(t)/t to t = 0);
  - spectral sum: int_0^Omega S(omega) d omega in closed form from the samples,
    (2/pi) sum_k c_k [sin(Omega t_k)/t_k^2 - Omega cos(Omega t_k)/t_k] (c_k: the
    weighted damped samples), for several cutoffs Omega up to the Nyquist
    frequency pi/dt of the output step (aliasing above). The cumulative sums
    show where the oscillator strength lies; at the Nyquist frequency the sum
    equals the initial slope up to the sampling.
Peaks: the local maxima of S(omega) (positions interpolated; strength = the
integral of S between the neighbouring minima). Lines closer than about
2 pi/T cannot be told apart by the Fourier transform of a run of length T;
--pade adds the Pade approximant of the damped series (Bruner, LaMaster,
Lopata, J. Chem. Theory Comput. 12, 3741 (2016)), which extrapolates the
signal and resolves such lines from the same samples (damping time
--pade-tau, default T), and lists its peaks too. Pade poles of a dense
spectrum (e.g. the discretized continuum of the box) can split or move when
the series changes slightly, so the lowest peaks are repeated from the first
90% and 80% of the samples: trust only lines that stay.

usage: python3 absorption_spectrum.py dipole.dat [dipole2.dat ...]
           [--tau TAU] [--omega-max W] [--d-omega DW] [--kappa K]
           [--direction A] [--electrons N] [--min-peak F] [--output FILE]
           [--pade] [--pade-tau TAU_P]
kappa, the kick direction (1, 2, 3 = x, y, z) and N are read from the header
of dipole.dat (P_FemKickTdhf); the options override them. tau defaults to
T/4. With several files (e.g. kicks along x, y and z) alpha is their average,
(alpha_xx + alpha_yy + alpha_zz)/3 for the orientation-averaged spectrum.
Writes FILE (default spectrum.dat): omega [a.u.], omega [eV], Re alpha,
Im alpha, sigma [a.u.], sigma [Mb], S(omega), up to W (default 2.0 a.u.);
with --pade also FILE with the suffix _pade (same columns, step DW/4)."""
import argparse, cmath, math, re, sys

HARTREE_EV = 27.211386245988
C_LIGHT = 137.035999084
AU_AREA_MB = 28.00285205  # bohr^2 in megabarn


def read_dipole(path):
    """Header values (kappa, direction, electrons) and the rows of dipole.dat."""
    header = {}
    times, rows = [], []
    number = r"([-+]?[0-9]*\.?[0-9]+(?:[EeDd][-+]?[0-9]+)?)"
    for line in open(path):
        if line.startswith("#"):
            for key, pattern in (("kappa", r"kappa\s*=\s*" + number),
                                 ("direction", r"direction\s*=\s*([123])"),
                                 ("electrons", r"electrons\s*=\s*([0-9]+)")):
                match = re.search(pattern, line)
                if match and key not in header:
                    header[key] = float(match.group(1).replace("D", "E").replace("d", "e"))
            continue
        if not line.strip():
            continue
        values = [float(x.replace("D", "E")) for x in line.split()]
        times.append(values[0])
        rows.append(values[1:4])
    return header, times, rows


def trapezoid_weights(times):
    n = len(times)
    return [0.5 * ((times[k] - times[k - 1] if k > 0 else 0.0) + (times[k + 1] - times[k] if k < n - 1 else 0.0))
            for k in range(n)]


def trk_sum(samples, cutoff):
    """(2/pi) int_0^cutoff omega Im alpha(omega) d omega for Im alpha = sum_k c_k sin(omega t_k)."""
    total = 0.0
    for t, c in samples:
        if t > 0.0:
            total += c * (math.sin(cutoff * t) / t ** 2 - cutoff * math.cos(cutoff * t) / t)
    return 2.0 * total / math.pi


def spectrum_row(omega, alpha):
    sigma = 4.0 * math.pi * omega * alpha.imag / C_LIGHT
    return (omega, omega * HARTREE_EV, alpha.real, alpha.imag, sigma, sigma * AU_AREA_MB, 2.0 * omega * alpha.imag / math.pi)


def write_spectrum(path, comment, rows):
    with open(path, "w") as f:
        f.write("# %s\n" % comment)
        f.write("# omega[a.u.] omega[eV] Re_alpha Im_alpha sigma[a.u.] sigma[Mb] S(omega)\n")
        for row in rows:
            f.write("%12.6f %12.6f %16.8e %16.8e %16.8e %16.8e %16.8e\n" % row)


def integrate(omegas, values):
    """Trapezoidal integral from omega = 0 (value 0 there) to the last omega."""
    return 0.5 * values[0] * omegas[0] + sum(0.5 * (values[k] + values[k + 1]) * (omegas[k + 1] - omegas[k])
                                             for k in range(len(values) - 1))


def find_peaks(omegas, S, minPeak):
    """Local maxima of S above minPeak * max(S): (position (parabolic interpolation), height,
    integral of S between the neighbouring minima, omega of those minima)."""
    sMax = max(S)
    peaks = []
    for i in range(1, len(S) - 1):
        if S[i] > S[i - 1] and S[i] >= S[i + 1] and S[i] > minPeak * sMax:
            curvature = S[i - 1] - 2.0 * S[i] + S[i + 1]
            shift = 0.5 * (S[i - 1] - S[i + 1]) / curvature if curvature else 0.0
            position = omegas[i] + shift * (omegas[i + 1] - omegas[i])
            lo, hi = i, i
            while lo > 0 and S[lo - 1] < S[lo]:
                lo -= 1
            while hi < len(S) - 1 and S[hi + 1] < S[hi]:
                hi += 1
            strength = sum(0.5 * (S[k] + S[k + 1]) * (omegas[k + 1] - omegas[k]) for k in range(lo, hi))
            peaks.append((position, S[i], strength, omegas[lo], omegas[hi]))
    return peaks


def print_peaks(omegas, S, minPeak, title):
    print("\n%s (local maxima above %.0f%% of the maximum %.4g):" % (title, 100 * minPeak, max(S)))
    print("   omega [a.u.]   omega [eV]       S(omega)   strength (integral of S between the neighbouring minima)")
    peaks = find_peaks(omegas, S, minPeak)
    for position, height, strength, lo, hi in peaks:
        print("   %12.5f %12.4f %14.6e   %.4f  (omega %.4f .. %.4f)"
              % (position, position * HARTREE_EV, height, strength, lo, hi))
    if not peaks:
        print("   none")


def pade_coefficients(g):
    """Pade approximant a(z)/b(z), b_0 = 1, both of degree N, of sum_k g_k z^k (k = 0 .. 2N):
    sum_{m=1..N} b_m g_{k-m} = -g_k for k = N+1 .. 2N (Gaussian elimination with partial pivoting),
    a_k = sum_{m=0..k} b_m g_{k-m}."""
    if len(g) % 2 == 0:
        g = g[:-1]
    n = (len(g) - 1) // 2
    rows = [[g[n + i - j] for j in range(1, n + 1)] + [-g[n + i]] for i in range(1, n + 1)]
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(rows[r][col]))
        rows[col], rows[pivot] = rows[pivot], rows[col]
        rowPivot = rows[col]
        if rowPivot[col] == 0.0:
            raise ZeroDivisionError("singular Pade system")
        for r in range(col + 1, n):
            row = rows[r]
            factor = row[col] / rowPivot[col]
            if factor:
                for c in range(col, n + 1):
                    row[c] -= factor * rowPivot[c]
    b = [0.0] * n
    for i in range(n - 1, -1, -1):
        b[i] = (rows[i][n] - sum(rows[i][j] * b[j] for j in range(i + 1, n))) / rows[i][i]
    b = [1.0] + b
    a = [sum(b[m] * g[k - m] for m in range(k + 1)) for k in range(n + 1)]
    return a, b


def horner(coefficients, z):
    value = 0j
    for c in reversed(coefficients):
        value = value * z + c
    return value


def pade_spectrum(series, dt, tauPade, omegas):
    """Spectrum rows of the Pade approximant of alpha(omega) = dt sum_k g_k z^k, z = e^{i omega dt},
    g_k = series_k e^{-k dt/tauPade} (series: d(t_k)/kappa from t = 0 on), and its coefficients."""
    g = [s * math.exp(-k * dt / tauPade) for k, s in enumerate(series)]
    a, b = pade_coefficients(g)
    rows = []
    for omega in omegas:
        z = cmath.exp(1j * omega * dt)
        rows.append(spectrum_row(omega, dt * horner(a, z) / horner(b, z)))
    return rows, a, b


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("files", nargs="+", help="dipole.dat of P_FemKickTdhf (several: average)")
    parser.add_argument("--tau", type=float, help="damping time of w(t) = exp(-t/tau) (default T/4)")
    parser.add_argument("--omega-max", type=float, default=2.0, help="largest omega of the output (a.u.)")
    parser.add_argument("--d-omega", type=float, help="omega step (default min(0.002, pi/(4T)))")
    parser.add_argument("--kappa", type=float, help="kick strength (default: header)")
    parser.add_argument("--direction", type=int, choices=(1, 2, 3), help="kick direction (default: header)")
    parser.add_argument("--electrons", type=float, help="number of electrons N (default: header)")
    parser.add_argument("--min-peak", type=float, default=0.01, help="peaks below this fraction of max S are not listed")
    parser.add_argument("--output", default="spectrum.dat", help="output file (default spectrum.dat)")
    parser.add_argument("--pade", action="store_true", help="also the Pade spectrum and its peaks")
    parser.add_argument("--pade-tau", type=float, help="damping time of the Pade series (default T)")
    args = parser.parse_args()

    runs = []
    for path in args.files:
        header, times, rows = read_dipole(path)
        kappa = args.kappa if args.kappa is not None else header.get("kappa")
        direction = args.direction if args.direction is not None else header.get("direction")
        electrons = args.electrons if args.electrons is not None else header.get("electrons")
        if kappa is None or direction is None:
            sys.exit("%s: kappa and the direction are not in the header, give --kappa and --direction" % path)
        if not kappa:
            sys.exit("%s: kappa = 0 (no kick)" % path)
        if len(times) < 3:
            sys.exit("%s: fewer than 3 samples" % path)
        if abs(times[0]) > 1e-9:
            print("note: %s starts at t = %g, not at the kick (t = 0)" % (path, times[0]))
        direction = int(direction)
        d = [r[direction - 1] - rows[0][direction - 1] for r in rows]
        runs.append(dict(path=path, kappa=kappa, direction=direction, electrons=electrons, times=times, d=d))

    T = min(r["times"][-1] for r in runs)
    dtMax = max(max(b - a for a, b in zip(r["times"], r["times"][1:])) for r in runs)
    omegaNyquist = math.pi / dtMax
    tau = args.tau if args.tau else T / 4.0
    dOmega = args.d_omega if args.d_omega else min(0.002, math.pi / T / 4.0)
    electrons = runs[0]["electrons"]
    if any(r["electrons"] != electrons for r in runs):
        print("note: the files have different electron numbers %s" % [r["electrons"] for r in runs])

    print("%d file(s): %s" % (len(runs), ", ".join("%s (kick along %s, kappa %.3e)" % (r["path"], "xyz"[r["direction"] - 1],
                                                                                         r["kappa"]) for r in runs)))
    print("T = %.3f a.u. (%d samples, dt = %.4g), tau = %.4g a.u.: Lorentzian half width 1/tau = %.4f a.u. = %.3f eV"
          % (T, len(runs[0]["times"]), dtMax, tau, 1.0 / tau, HARTREE_EV / tau))
    print("Fourier resolution 2 pi/T = %.4f a.u.; Nyquist frequency pi/dt = %.3f a.u. = %.1f eV"
          % (2.0 * math.pi / T, omegaNyquist, omegaNyquist * HARTREE_EV))

    # Damped, weighted samples of the average: alpha(omega) = sum_k c_k e^{i omega t_k}
    samples = []
    for r in runs:
        weights = trapezoid_weights(r["times"])
        for t, w, dk in zip(r["times"], weights, r["d"]):
            if t <= T + 1e-9:
                samples.append((t, w * dk * math.exp(-t / tau) / (r["kappa"] * len(runs))))

    omegas = []
    k = 1
    while k * dOmega <= args.omega_max + 1e-12 and k * dOmega < omegaNyquist:
        omegas.append(k * dOmega)
        k += 1
    rows = [spectrum_row(omega, sum(c * cmath.exp(1j * omega * t) for t, c in samples)) for omega in omegas]
    comment = "absorption_spectrum.py %s: tau = %.6g, kappa = %s, directions = %s" % (
        " ".join(r["path"] for r in runs), tau, [r["kappa"] for r in runs], [r["direction"] for r in runs])
    write_spectrum(args.output, comment, rows)
    print("wrote %s (%d frequencies, 0 < omega <= %.3f a.u., d omega = %.4g)" % (args.output, len(rows), omegas[-1], dOmega))

    print("\nstatic polarizability of the damped signal (lines count with f/(omega^2 + 1/tau^2)): alpha(0) = %.6f a.u."
          % sum(c for t, c in samples))
    S = [row[6] for row in rows]
    print_peaks(omegas, S, args.min_peak, "peaks of S(omega)")

    if args.pade:
        sameGridQ = all(len(r["times"]) == len(runs[0]["times"]) and
                        max(abs(a - b) for a, b in zip(r["times"], runs[0]["times"])) < 1e-9 for r in runs)
        steps = [b - a for a, b in zip(runs[0]["times"], runs[0]["times"][1:])]
        if not sameGridQ or max(steps) - min(steps) > 1e-6 * max(steps):
            print("\nPade: needs equidistant samples on the same times in every file; skipped")
        else:
            dt = steps[0]
            times = runs[0]["times"]
            series = [sum(r["d"][k] / r["kappa"] for r in runs) / len(runs) for k in range(len(times))]
            omegasPade = [k * dOmega / 4.0 for k in range(1, 4 * len(omegas) + 1)]
            tauPade = args.pade_tau if args.pade_tau else T
            rowsPade, a, b = pade_spectrum(series, dt, tauPade, omegasPade)
            outputPade = re.sub(r"(\.dat)?$", "_pade.dat", args.output, count=1)
            commentPade = "absorption_spectrum.py %s: Pade approximant, tau = %.6g, kappa = %s, directions = %s" % (
                " ".join(r["path"] for r in runs), tauPade, [r["kappa"] for r in runs], [r["direction"] for r in runs])
            write_spectrum(outputPade, commentPade, rowsPade)
            print("\nPade approximant of degree %d (damping time %.4g a.u., half width %.4f a.u.): wrote %s"
                  % (len(b) - 1, tauPade, 1.0 / tauPade, outputPade))
            print("static polarizability of the Pade approximant: alpha(0) = %.6f a.u." % (dt * sum(a) / sum(b)))
            print_peaks(omegasPade, [row[6] for row in rowsPade], args.min_peak, "peaks of the Pade S(omega)")
            # The Pade poles can split or move when the series changes slightly (dense spectra, e.g. the
            # discretized continuum of a box): repeat with shortened series and compare
            print("   stability (lowest Pade peaks from the first part of the samples; trust only stable lines):")
            for fraction in (0.9, 0.8):
                n = int(round(fraction * (len(series) - 1))) + 1
                tauPart = args.pade_tau if args.pade_tau else times[n - 1] - times[0]
                rowsPart, a, b = pade_spectrum(series[:n], dt, tauPart, omegasPade)
                peaks = find_peaks(omegasPade, [row[6] for row in rowsPart], args.min_peak)
                print("      T = %8.3f a.u.: %s" % (times[n - 1] - times[0], "  ".join("%.5f" % p[0] for p in peaks[:4])))

    # TRK sum rule
    print("\nTRK sum rule (integral of S over all omega = N%s):" % (" = %g" % electrons if electrons else ""))
    for r in runs:
        t, d = r["times"], r["d"]
        slope = (4.0 * d[1] / t[1] - d[2] / t[2]) / 3.0
        print("   initial slope d<r_%s>/dt(0+)/kappa = %.8f (%s; first two samples, Richardson)"
              % ("xyz"[r["direction"] - 1], slope / r["kappa"], r["path"]))
    cutoffs = [c for c in (0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0) if c < omegaNyquist] + [omegaNyquist]
    print("   int_0^Omega S(omega) d omega (closed form from the samples):")
    for cutoff in cutoffs:
        value = trk_sum(samples, cutoff)
        note = " (Nyquist)" if cutoff is cutoffs[-1] else ""
        print("      Omega = %8.3f a.u.%-10s %10.6f%s" % (cutoff, note, value,
                                                         "   (%.4f of N)" % (value / electrons) if electrons else ""))
    print("   check: trapezoidal integral of the written S up to %.3f a.u. = %.6f, closed form %.6f"
          % (omegas[-1], integrate(omegas, S), trk_sum(samples, omegas[-1])))


if __name__ == "__main__":
    main()
