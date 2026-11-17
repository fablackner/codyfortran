# Printer Module — Developer Onboarding Guide

## Overview

The **Printer** modules write observables, densities, reduced density matrices
and spectra to the screen and to files. They are plain procedure libraries like
Utils: no `Fabricate`, no `Setup`, no JSON block. Unlike Utils they read the
state of the simulation (grid, orbitals, method), so they sit above the domain
modules and are called by programs only. No module outside `src/Printer/`
imports a printer.

---

## Directory Structure

```
src/Printer/
├── A_Printer.md                               # This file
├── Density/
│   ├── M_Printer_Density.f90                  # One- and two-body densities in the orbital basis
│   └── Linear/
│       └── M_Printer_Density_Linear.f90       # Densities on a 1D grid
├── NatOrb/
│   └── M_Printer_NatOrb.f90                   # Natural-orbital occupations
├── Observable/
│   ├── M_Printer_Observable.f90               # Norm, orthonormality, orbital energy matrix, energy
│   ├── Fem/
│   │   └── M_Printer_Observable_Fem.f90       # Dipole, norm, dipole velocity on grid.fem
│   ├── Lattice/
│   │   └── M_Printer_Observable_Lattice.f90   # Site densities
│   ├── Linear/
│   │   └── M_Printer_Observable_Linear.f90    # Dipole and potential on 1D grids
│   └── Ylm/
│       └── M_Printer_Observable_Ylm.f90       # Dipole and dipole acceleration on Ylm grids
├── Rdm/
│   ├── M_Printer_Rdm.f90                      # RDM dumps and symmetry diagnostics
│   └── Linear/
│       └── M_Printer_Rdm_Linear.f90           # 1-RDM on a 1D grid
└── Spectrum/
    └── M_Printer_Spectrum.f90                 # Eigenvalues of a diagonalizer
```

The first level is the quantity. The second level is the grid family, for
output that needs the coordinates of that family. A module is a single file
without a submodule: a printer binds nothing.

---

## Modules

| Module | Public routines | Notes |
|--------|-----------------|-------|
| `M_Printer_Density` | `Printer_Density_DumpOneBodyDensity`, `_DumpTwoBodyDensity` | Orbital basis |
| `M_Printer_Density_Linear` | `Printer_Density_Linear_DumpOneBodyDensityOnGrid`, `_DumpTwoBodyDensityOnGrid` | Linear grids |
| `M_Printer_NatOrb` | `Printer_NatOrb_DumpNatOrbOccupation` | Occupations from the 1-RDM |
| `M_Printer_Observable` | `Printer_Observable_DumpNorm`, `_DumpOrthonormality`, `_DumpOrbitalEnergyMatrix`, `_DumpEnergy` | Method-aware: uses `M_Method`/`M_Method_Mb_OrbBased` (Apply*Op, `Method_Mb_OrbBased_FillInverseRdm1` for the orbital energy matrix), so not a stateless utility; one inner product per orbital (norms) and per pair (overlap matrix) |
| `M_Printer_Observable_Linear` | `Printer_Observable_Linear_DumpDipole` (line "time Re⟨x⟩ Im⟨x⟩", optional `dipole_` out), `_DumpPotentialFullExpansion` | Linear grids; ⟨x⟩ = Σ_ij ρ_s(i,j)⟨φ_j\|x\|φ_i⟩ with ρ_s = `Method_Mb_OrbBased_FillSpatialRdm1` of rdm1 |
| `M_Printer_Observable_Lattice` | `Printer_Observable_Lattice_DumpDensityCenter`, `_DumpDensityOnSite` | Lattice grids |
| `M_Printer_Observable_Ylm` | `Printer_Observable_Ylm_DumpDipole` (optional `dipole_` out), `Printer_Observable_Ylm_DumpDipoleAcceleration` (optional `acceleration_` out) | Ylm grids; z = r√(4π/3)Y₁₀ applied as spatial product, with ρ_s = `Method_Mb_OrbBased_FillSpatialRdm1` of rdm1; acceleration form −Z⟨cosθ/r²⟩ |
| `M_Printer_Observable_Fem` | `Printer_Observable_Fem_DumpDipole(rdm1, orbs, filename, toScreenQ, time, dipole_)` (line "time ⟨x⟩ ⟨y⟩ ⟨z⟩", optional `dipole_(3)` out), `Printer_Observable_Fem_DumpNorm(..., norm_)` (line "time N", N = Tr(ρ S)), `Printer_Observable_Fem_DumpDipoleVelocity(..., velocity_)` (line "time v_x v_y v_z", the velocity form v_a = i Tr(ρ_s [T, r_a]) of the dipole, d⟨r_a⟩/dt wherever the equations of motion obey Ehrenfest's theorem for r_a, e.g. TDHF; optional `velocity_(3)` out) | grid.fem; Method-aware like `M_Printer_Observable`: ρ_s = `Method_Mb_OrbBased_FillSpatialRdm1` of rdm1 (rdm1 itself, or the sum of the two spin blocks for restricted storage), body type of stored orbital i = `Method_Mb_OrbBased_bodyTypeOfOrb(i)`. Dipole and N: one pass over the nodes with the density ρ = Σ_ij ρ_s(i,j) φ_j* φ_i and the nodal weights (no orthonormality assumed), node coordinates from `Grid_Fem_FillNodeCoordinate`. Velocity: Re(i Σ_ij ρ_s(i,j)(⟨Tφ_j\|r_aφ_i⟩ − ⟨r_aφ_j\|Tφ_i⟩)) with `SysKinetic_MultiplyWithKineticOp` (one T application per orbital, 2 nOS² inner products per axis); the potential, the length gauge and the interaction commute with r_a, so after a kick e^{iκr_a} it gives the TRK value Nκ/m on the grid (P_KickTdhfFem and FemKick.inc take their initial slope from it); early return without output request |
| `M_Printer_Rdm` | `Printer_Rdm_DumpRdm1`, `_DumpRdm1AbsMatrixForm`, `_DumpRdm2`, `_DumpRdm3`, `_DumpRdm2Symmetries` | `_DumpRdm2Symmetries` reads the particle and orbital numbers of `M_Method_Mb` and `M_Method_Mb_OrbBased` |
| `M_Printer_Rdm_Linear` | `Printer_Rdm_Linear_DumpRdm1FullExpansion` | Linear grids |
| `M_Printer_Spectrum` | `Printer_Spectrum_DumpSpectrum` | Eigenvalues of `DiagonalizerList` |

---

## Conventions

- Names: module `M_Printer_<Quantity>[_<Family>]`, routines
  `Printer_<Quantity>[_<Family>]_Dump<What>`.
- Arguments: the data arrays, then `filename`, `toScreenQ` and an optional
  `time`. Nothing is computed if neither the screen nor a file is requested.
- Files are appended to (`position="append"`), one line per call for a time
  series.
- Established columns keep their position. New columns are appended.

Adding a routine and plotting its output: `.github/skills/printer-system/SKILL.md`.
