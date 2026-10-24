"""HHG spectrum of a run: |FT[a(t) w(t)]|^2 (dipole acceleration, Hann
window) and omega^4 |FT[z(t) w(t)]|^2 (length form) versus harmonic order of
omega0 = 0.057, in steps of 0.1 up to maxOrder (default 300).

usage: python3 hhg_spectrum.py <observables.dat> [maxOrder]
Writes hhg.dat (order, acceleration form, length form) in the working
directory."""
import cmath, math, sys
omega0 = 0.057
rows = [list(map(float, l.split())) for l in open(sys.argv[1]) if not l.startswith("#")]
t = [r[0] for r in rows]; z = [r[3] for r in rows]; a = [r[4] for r in rows]
n = len(t); T = t[-1] - t[0]; dt = T / (n - 1)
w = [math.sin(math.pi * (ti - t[0]) / T) ** 2 for ti in t]
zm = sum(z) / n
out = []
dOrder = 0.1
order = dOrder
maxOrder = float(sys.argv[2]) if len(sys.argv) > 2 else 300.0
while order <= maxOrder and order * omega0 < math.pi / dt:
    om = order * omega0
    fa = sum(a[k] * w[k] * cmath.exp(-1j * om * t[k]) for k in range(n)) * dt
    fz = sum((z[k] - zm) * w[k] * cmath.exp(-1j * om * t[k]) for k in range(n)) * dt
    out.append((order, abs(fa) ** 2, om ** 4 * abs(fz) ** 2))
    order += dOrder
with open("hhg.dat", "w") as f:
    f.write("# harmonic order, |a(w)|^2, w^4 |z(w)|^2\n")
    for o, sa, sz in out:
        f.write("%8.2f %14.6e %14.6e\n" % (o, sa, sz))
