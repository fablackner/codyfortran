---
name: setup-simulation-app
description: "Use when creating a new CodyFortran simulation app scaffold under app/, including folder creation and a matched Fortran program (.f90) + JSON input (.json). Provides a complete starter workflow (fabricate, setup, propagate, output) and reusable templates. Keywords: app scaffold, simulation setup, create app folder, generate f90 json, new simulation program."
argument-hint: "[category] [system] [groundstate|dynamics] [method] [program base name]"
---

# Setup Complete Simulation App Scaffold

Use this skill when you want to create a new runnable simulation app under `app/` with:

1. A folder `app/<category>/<system>/<stage>/<method>/` (layout: `app/README.md`).
2. A program source file `P_*.f90`.
3. A matching JSON input `P_*.json`.

A folder holds exactly one program and one JSON input with the same base
name. A variant (another pulse, another orbital space) gets its own folder,
e.g. `dynamics/mctdhf/weak4c/`. Programs run in their folder; a dynamics
input reads the ground state from its folder with `orbsInit.load.dir`
(and `coeffsInit.load.dir`, `twoRdmInit.load.dir`).

Analysis scripts (`*.py`) of a system go to `app/<category>/<system>/analysis/`.

## Required References

- [Project architecture and module index](../../../AGENTS.md)
- [Build and executable wiring](../../../CMakeLists.txt)
- [Fortran template](./templates/P_Simulation_Scaffold.f90)
- [JSON template](./templates/P_Simulation_Scaffold.json)
- [Optional scaffold script](./scripts/create_simulation_app.sh)

## Input Contract

Provide five values:

1. `category`: `atoms`, `molecules`, or a new one.
2. `system`: the physical system (example: `be`, `h2o`, `chain6`).
3. `stage`: `groundstate` or `dynamics`.
4. `method`: method folder (example: `scf`, `mcscf`, `mctdhf`, `tdci`, `td2rdm`).
5. `programBase`: base program file name without extension (example: `P_TdciPropagate`).

Resulting files:

- `app/<category>/<system>/<stage>/<method>/<programBase>.f90`
- `app/<category>/<system>/<stage>/<method>/<programBase>.json`

## Standard Workflow

1. Validate naming.
- `programBase` should start with `P_`.
- Use ASCII letters, digits, `_` for folder and file names.

2. Create target directory.
- Make `app/<category>/<system>/<stage>/<method>/` if missing.

3. Create `.f90` from template.
- Copy [Fortran template](./templates/P_Simulation_Scaffold.f90).
- Replace `__PROGRAM_NAME__` with `<programBase>`.
- Replace `__JSON_REL_PATH__` with `app/<category>/<system>/<stage>/<method>/<programBase>.json`.

4. Create `.json` from template.
- Copy [JSON template](./templates/P_Simulation_Scaffold.json).
- Keep schema complete for a first runnable lattice TD-CI propagation.

5. Optional script-based creation.
- Use [create_simulation_app.sh](./scripts/create_simulation_app.sh) for one-command scaffolding.

6. Build and run.

```bash
source setBuildVars.sh -c gnu -t release
cmake -B build && cmake --build build && cmake --install build
./app/<category>/<system>/<stage>/<method>/<programBase>.exe
```

## Notes

- `CMakeLists.txt` already discovers `app/**/*.f90`, so no CMake edits are needed for new app programs.
- The scaffold program explicitly loads JSON using project-relative path, so execution does not depend on current working directory.
- After scaffold creation, customize JSON blocks (`grid`, `sysPotential`, `method`, etc.) for the target physics case.

## Example Invocation

`/setup-simulation-app hubbard chain6 dynamics tdci P_TdciPropagate`

Expected output:

- `app/hubbard/chain6/dynamics/tdci/P_TdciPropagate.f90`
- `app/hubbard/chain6/dynamics/tdci/P_TdciPropagate.json`

## Done Checklist

- Target app directory exists.
- `.f90` and `.json` files exist with matching base names.
- Program name and JSON path placeholders were replaced.
- Project builds with the new app executable.