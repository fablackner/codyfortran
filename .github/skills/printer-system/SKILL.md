---
name: printer-system
description: "Use when working on CodyFortran printer modules in src/Utils, including adding new printer routines, extending existing dump formats, wiring printer calls in apps/tests, and generating plots from printer output with gnuplot. Keywords: printer, dump, observable, density, rdm, output file, gnuplot, .gnup."
argument-hint: "[target printer module] [new observable or format] [plot goal]"
---

# Extend Printer System And Generate Gnuplot Plots

Use this skill when you need to:

1. Add a new printer routine in the Utils printer modules.
2. Extend or adjust existing output files and columns.
3. Wire printer output into simulation programs or tests.
4. Generate publication-ready plots from printer data using gnuplot.

## Required Context Before Editing

Load these files first:

- [Framework overview](../../../AGENTS.md)
- [Printer density base](../../../src/Utils/M_Utils_PrinterDensity.f90)
- [Printer density linear](../../../src/Utils/M_Utils_PrinterDensityLinear.f90)
- [Printer observable base](../../../src/Utils/M_Utils_PrinterObservable.f90)
- [Printer observable linear](../../../src/Utils/M_Utils_PrinterObservableLinear.f90)
- [Printer observable lattice](../../../src/Utils/M_Utils_PrinterObservableLattice.f90)
- [Printer RDM base](../../../src/Utils/M_Utils_PrinterRdm.f90)
- [Printer RDM linear](../../../src/Utils/M_Utils_PrinterRdmLinear.f90)
- [Printer natural orbital](../../../src/Utils/M_Utils_PrinterNatOrb.f90)
- [Printer spectrum](../../../src/Utils/M_Utils_PrinterSpectrum.f90)

For plotting conventions, inspect:

- [Hubbard line plot example](../../../app/cnnHubbard/checkModel/density_site3.gnup)
- [Hubbard heatmap example](../../../app/cnnHubbard/checkModel/colorDensity.gnup)

## Printer Module Map

- `M_Utils_PrinterDensity*`: one-body/two-body densities.
- `M_Utils_PrinterObservable*`: energy, norm, dipole, potential, lattice-site observables.
- `M_Utils_PrinterRdm*`: full RDM dumps and RDM symmetry diagnostics.
- `M_Utils_PrinterNatOrb`: natural-orbital occupations from 1-RDM.
- `M_Utils_PrinterSpectrum`: diagonalizer eigenvalue export.

Use base modules for representation-independent output and specialized modules (`Linear`, `Lattice`) for geometry-dependent projections.

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
- Add `use M_Utils_Printer...` in the caller.
- Call routine after state/RDM refresh so printed values match current time step.
- Prefer deterministic output filenames.

## Data-Format Cheatsheet For Plotting

- `PrinterObservable_DumpEnergy`: column 1 = time, columns 2-3 = complex energy (real/imag parts in scientific format).
- `PrinterObservableLinear_DumpDipole`: column 1 = time, columns 2-3 = complex dipole.
- `PrinterObservableLattice_DumpDensityOnSite`: column 1 = time, columns 2..N+1 = site densities.
- `PrinterSpectrum_DumpSpectrum`: column 1 = eigenvalue index, column 2 = eigenvalue.
- `PrinterDensityLinear_DumpOneBodyDensityOnGrid`: column 1 = x, columns 2-3 = complex density.

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
source setBuildVars.sh -c gnu -t release
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
