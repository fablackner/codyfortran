---
name: extend-module
description: "Use when adding or extending CodyFortran modules/backends (for example: a new Grid type, new SysPotential model, or new implementation variant). Follow this workflow to add files, wire Fabricate dispatchers, bind procedure pointers, configure JSON keys, and validate with tests. Keywords: extend module, new grid, new potential, Fabricate chain, procedure pointer wiring, JSON dispatch."
argument-hint: "[target module] [new type] [json path]"
---

# Extend CodyFortran Module

Use this skill when implementing new functionality in CodyFortran that follows the existing Interface Module + Implementation Submodule architecture.

## Core Pattern

CodyFortran extensions usually follow this chain:

1. Create a new folder and `M_*.f90` + `S_*.f90` pair.
2. Parse JSON keys in `*_Fabricate` and store parameters in module data.
3. Bind exported procedure pointers to concrete routines.
4. Register the new type in the parent dispatcher `S_*.f90` (`if`/`else if` branch).
5. Add/update tests and JSON input.

This pattern applies to Grid, SysPotential, SysInteraction, SysKinetic, Reconstruction, and many other modules.

## Required Context Before Editing

Load the relevant module docs before implementing:

- [Framework overview](../../../AGENTS.md)
- [Grid onboarding](../../../src/Grid/A_Grid.md)
- [SysPotential onboarding](../../../src/SysPotential/A_SysPotential.md)
- [Build system](../../../CMakeLists.txt)

Then inspect sibling implementations in the same folder and copy their structure.

## Extension Workflow (Canonical)

1. Identify insertion level in the hierarchy.
- New family (for example `grid.cylindrical`) means wiring in top-level dispatcher (for example `src/Grid/S_Grid.f90`).
- New variant inside a family (for example `sysPotential.linear.myPotential`) means wiring in family dispatcher (for example `src/SysPotential/Linear/S_SysPotential_Linear.f90`).
- Optional backend implementation level (for example `stdImpl`) may require one more nested dispatcher.

2. Create new files by copying the nearest sibling.
- Keep naming convention aligned with path.
- Example for new linear potential `myPotential`:
  - `src/SysPotential/Linear/MyPotential/M_SysPotential_Linear_MyPotential.f90`
  - `src/SysPotential/Linear/MyPotential/S_SysPotential_Linear_MyPotential.f90`
  - Optional concrete backend:
    - `src/SysPotential/Linear/MyPotential/StdImpl/M_SysPotential_Linear_MyPotential_StdImpl.f90`
    - `src/SysPotential/Linear/MyPotential/StdImpl/S_SysPotential_Linear_MyPotential_StdImpl.f90`

3. Implement `*_Fabricate` in the new submodule.
- Call `Say_Fabricate("...")` with the exact JSON path.
- Read parameters with `Json_Get("path", default)`.
- Validate shapes/ranges and `error stop` with clear messages.
- Set module flags if needed (for `SysPotential`: `SysPotential_timeIndependentQ`, `SysPotential_bodyTypeIndependentQ`).

4. Bind procedure pointers at the correct level.
- Family-level submodule may bind generic operations (example: `SysPotential_MultiplyWithExternalPotential`).
- Concrete implementation submodule binds actual compute routine (example: `SysPotential_FillExternalPotential => FillExternalPotential`).
- For Grid variants, bind pointers such as `Grid_Setup => Setup` and family-level `Grid_InnerProduct => InnerProduct` where appropriate.

5. Register new type in parent dispatcher.
- Add `use M_..._YourType` in dispatcher submodule.
- Add branch:

```fortran
if (Json_GetExistence("<json.path.to.newType>")) then
  call <YourType>_Fabricate
else if (...)
  ...
else
  error stop "...missing one of: ..."
end if
```

- Keep branch ordering and error message list consistent.

6. If introducing a new top-level family, wire top-level module/factory too.
- Update top-level dispatcher (for example `src/Grid/S_Grid.f90` or `src/SysPotential/S_SysPotential.f90`).
- Ensure required `use M_<Family>` is present in that dispatcher.

7. Add JSON usage example.
- New potential example:

```json
"sysPotential": {
  "linear": {
    "myPotential": {
      "strength": 1.0,
      "stdImpl": {}
    }
  }
}
```

- New grid variant example:

```json
"grid": {
  "linear": {
    "myGrid": {
      "nPoints": 256
    }
  }
}
```

8. Add or update tests.
- Copy a nearby `test/T_*.f90` and matching `test/*.json` input, then adapt minimal knobs.
- Prefer small deterministic tests that verify fabrication wiring and one numerical sanity check.

9. Build and run focused tests.
- Standard project commands:

```bash
source setBuildVars.sh -c gnu -t release
cmake -B build && cmake --build build && cmake --install build
cd build && ctest -R T_
```

- If you added a specific feature test, run it directly with `ctest -R <test_name_fragment>`.

## Grid-Specific Notes

When adding a new Grid family or variant:

1. Maintain metric consistency: any `Grid_InnerProduct` implementation must apply correct weights/Jacobians.
2. Ensure `Grid_nPoints` is set during fabricate before setup-dependent allocations.
3. Ensure `Grid_Setup` allocates all coordinates/weights consumed by downstream operators.
4. Preserve shared top-level bindings in `src/Grid/S_Grid.f90`:
- `Grid_Orthonormalize => Orthonormalize`
- `Grid_ProjectOnSubspace => ProjectOnSubspace`

## SysPotential-Specific Notes

When adding a new potential model:

1. Keep JSON path shape consistent: `sysPotential.<gridFamily>.<model>[.<impl>]`.
2. Bind multiplication at family level when diagonal in basis.
3. Bind fill routine at concrete implementation level.
4. Set independence flags correctly for caching/optimization.
- Time-dependent models: `SysPotential_timeIndependentQ = .false.`
- Body-type dependent models: `SysPotential_bodyTypeIndependentQ = .false.`

## Build-System Note

This repository uses `file(GLOB_RECURSE ... "src/*.f90")` in [CMakeLists.txt](../../../CMakeLists.txt), so new `.f90` files are discovered automatically after reconfigure. You generally do not need to edit per-module CMake files.

## Done Checklist

- New files follow naming and folder conventions.
- Parent dispatcher has `use` + JSON branch + updated error message.
- Procedure pointers are bound exactly once at the correct abstraction level.
- JSON example for the new feature exists.
- At least one test/json pair added or updated.
- Build + relevant tests pass.
