---
name: setup-simulation-app
description: "Use when creating a new CodyFortran simulation app scaffold under app/, including folder creation and a matched Fortran program (.f90) + JSON input (.json). Provides a complete starter workflow (fabricate, setup, propagate, output) and reusable templates. Keywords: app scaffold, simulation setup, create app folder, generate f90 json, new simulation program."
argument-hint: "[app group] [scenario folder] [program base name]"
---

# Setup Complete Simulation App Scaffold

Use this skill when you want to create a new runnable simulation app under `app/` with:

1. A dedicated folder path.
2. A program source file `P_*.f90`.
3. A matching JSON input `P_*.json`.

## Required References

- [Project architecture and module index](../../../AGENTS.md)
- [Build and executable wiring](../../../CMakeLists.txt)
- [Fortran template](./templates/P_Simulation_Scaffold.f90)
- [JSON template](./templates/P_Simulation_Scaffold.json)
- [Optional scaffold script](./scripts/create_simulation_app.sh)

## Input Contract

Provide three values:

1. `appGroup`: top-level app bucket under `app/` (example: `myProject`).
2. `scenario`: scenario folder under that group (example: `baseline`).
3. `programBase`: base program file name without extension (example: `P_MySystem_TdciPropagate`).

Resulting files:

- `app/<appGroup>/<scenario>/<programBase>.f90`
- `app/<appGroup>/<scenario>/<programBase>.json`

## Standard Workflow

1. Validate naming.
- `programBase` should start with `P_`.
- Use ASCII letters, digits, `_` for folder and file names.

2. Create target directory.
- Make `app/<appGroup>/<scenario>/` if missing.

3. Create `.f90` from template.
- Copy [Fortran template](./templates/P_Simulation_Scaffold.f90).
- Replace `__PROGRAM_NAME__` with `<programBase>`.
- Replace `__JSON_REL_PATH__` with `app/<appGroup>/<scenario>/<programBase>.json`.

4. Create `.json` from template.
- Copy [JSON template](./templates/P_Simulation_Scaffold.json).
- Keep schema complete for a first runnable lattice TD-CI propagation.

5. Optional script-based creation.
- Use [create_simulation_app.sh](./scripts/create_simulation_app.sh) for one-command scaffolding.

6. Build and run.

```bash
source setBuildVars.sh -c gnu -t release
cmake -B build && cmake --build build && cmake --install build
./app/<appGroup>/<scenario>/<programBase>.exe
```

## Notes

- `CMakeLists.txt` already discovers `app/**/*.f90`, so no CMake edits are needed for new app programs.
- The scaffold program explicitly loads JSON using project-relative path, so execution does not depend on current working directory.
- After scaffold creation, customize JSON blocks (`grid`, `sysPotential`, `method`, etc.) for the target physics case.

## Example Invocation

`/setup-simulation-app hubbard1d sweepA P_Hubbard1d_TdciPropagate`

Expected output:

- `app/hubbard1d/sweepA/P_Hubbard1d_TdciPropagate.f90`
- `app/hubbard1d/sweepA/P_Hubbard1d_TdciPropagate.json`

## Done Checklist

- Target app directory exists.
- `.f90` and `.json` files exist with matching base names.
- Program name and JSON path placeholders were replaced.
- Project builds with the new app executable.