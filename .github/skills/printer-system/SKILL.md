---
name: printer-system
description: "Use when working on CodyFortran printer modules in src/Printer, including adding new printer routines, extending existing dump formats, wiring printer calls in apps/tests, and generating plots from printer output with gnuplot. Keywords: printer, dump, observable, density, rdm, output file, gnuplot, .gnup."
argument-hint: "[target printer module] [new observable or format] [plot goal]"
---

# Extend Printer System And Generate Gnuplot Plots

Use this skill when you need to:

1. Add a new printer routine in the Printer modules.
2. Extend or adjust existing output files and columns.
3. Wire printer output into simulation programs or tests.
4. Generate publication-ready plots from printer data using gnuplot.

## Required Context Before Editing

Load these files first:

- [Framework overview](../../../AGENTS.md)
- [Printer module doc](../../../src/Printer/A_Printer.md)
- [Printer density base](../../../src/Printer/Density/M_Printer_Density.f90)
- [Printer density linear](../../../src/Printer/Density/Linear/M_Printer_Density_Linear.f90)
- [Printer observable base](../../../src/Printer/Observable/M_Printer_Observable.f90)
- [Printer observable linear](../../../src/Printer/Observable/Linear/M_Printer_Observable_Linear.f90)
- [Printer observable lattice](../../../src/Printer/Observable/Lattice/M_Printer_Observable_Lattice.f90)
- [Printer observable fem](../../../src/Printer/Observable/Fem/M_Printer_Observable_Fem.f90)
- [Printer RDM base](../../../src/Printer/Rdm/M_Printer_Rdm.f90)
- [Printer RDM linear](../../../src/Printer/Rdm/Linear/M_Printer_Rdm_Linear.f90)
- [Printer natural orbital](../../../src/Printer/NatOrb/M_Printer_NatOrb.f90)
- [Printer spectrum](../../../src/Printer/Spectrum/M_Printer_Spectrum.f90)

For plotting conventions, inspect:

- [HHG spectrum plot](../../../app/atoms/be/analysis/plot_hhg.py)

## Printer Module Map

- `M_Printer_Density[_Linear]`: one-body/two-body densities.
- `M_Printer_Observable[_<Family>]`: energy, norm, dipole, potential, lattice-site observables.
- `M_Printer_Rdm[_Linear]`: full RDM dumps and RDM symmetry diagnostics.
- `M_Printer_NatOrb`: natural-orbital occupations from 1-RDM.
- `M_Printer_Spectrum`: diagonalizer eigenvalue export.

Use base modules for representation-independent output and specialized modules (`Linear`, `Lattice`, `Ylm`, `Fem`) for geometry-dependent projections.

## Workflow For Adding New Printer Features

1. Pick the right module level.
- Generic quantity independent of geometry: add to base module.
- Quantity tied to `Grid_Linear_*` or `Grid_Lattice_*`: add to corresponding specialized module.

2. Copy a neighboring routine and keep interface conventions.
- Keep argument order consistent: data arrays, `filename`, `toScreenQ`, optional `time`.
- Preserve `complex(R64)` and `integer(I32)` typing style.
- Keep `contiguous` on array arguments where existing routines use it.

3. Compute once and guard early when no output is requested.
- Use pattern: `needOutput = toScreenQ .or. (filename .ne. "")`.
- Return early when no output is needed.

4. Keep screen and file sections separated.
- Use clear `Print to Screen` and `Print to file` blocks.
- Preserve append semantics: `status="unknown", position="append"`.
- Keep open error handling pattern with `iostat` + `iomsg`, then `error stop`.

5. Keep data layout stable.
- Do not reorder existing columns for established files.
- If adding columns, append them and document meaning.

6. Wire into caller program/test.
- Add `use M_Printer_...` in the caller.
- Call routine after state/RDM refresh so printed values match current time step.
- Prefer deterministic output filenames.

## Data-Format Cheatsheet For Plotting

- `Printer_Observable_DumpEnergy`: column 1 = time, columns 2-3 = complex energy (real/imag parts in scientific format).
- `Printer_Observable_Linear_DumpDipole`: column 1 = time, columns 2-3 = complex dipole.
- `Printer_Observable_Fem_DumpDipole`: column 1 = time, columns 2-4 = <x>, <y>, <z> (grid.fem).
- `Printer_Observable_Fem_DumpNorm`: column 1 = time, column 2 = N = Tr(rho S) (grid.fem).
- `Printer_Observable_Fem_DumpDipoleVelocity`: column 1 = time, columns 2-4 = v_x, v_y, v_z = i Tr(rho [T, r_a]), the velocity form of the dipole (grid.fem).
- `Printer_Observable_Lattice_DumpDensityOnSite`: column 1 = time, columns 2..N+1 = site densities.
- `Printer_Spectrum_DumpSpectrum`: column 1 = eigenvalue index, column 2 = eigenvalue.
- `Printer_Density_Linear_DumpOneBodyDensityOnGrid`: column 1 = x, columns 2-3 = complex density.

Always confirm exact column count by checking one output line before plotting.

## Gnuplot Workflow

1. Start from templates in this skill:
- [Line plot template](./templates/gnuplot-lines.gnup)
- [Lattice heatmap template](./templates/gnuplot-lattice-heatmap.gnup)

2. Set paths and columns.
- Keep time on x-axis (`using 1:<column>`).
- For lattice files, use the AWK flatten command from the heatmap template.

3. Render plots.

```bash
gnuplot my_plot.gnup
```

4. Verify artifacts.

```bash
ls -lh *.png
```

## Optional Resource For New Routines

Use [printer routine template](./templates/printer-routine-template.f90) as a starting point, then adapt to your target observable.

## Validation Commands

```bash
source setCodyVars.sh -c gnu -t release
cmake -B build && cmake --build build && cmake --install build
cd build && ctest -R T_
```

If you changed one simulation/test program, run a focused regex with `ctest -R <name_fragment>`.

## Done Checklist

- New routine added in the correct printer module level.
- Caller program/test imports and invokes the routine at the right point.
- Output file format documented and stable.
- At least one `.gnup` script added or updated for the new output.
- Build and relevant tests pass.
