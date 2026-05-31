---
name: writeTestCase
description: "Use when creating new CodyFortran tests. Supports two categories: component tests in test/componentTests where the folder tree encodes the module/backend under test, and full simulation tests in test/simulationTests for end-to-end workflows. Includes naming, placement, JSON pairing, and validation workflow. Keywords: new test case, componentTests, simulationTests, T_*.f90, test JSON."
argument-hint: "[component|simulation] [target module or system] [scenario]"
---

# Write CodyFortran Test Cases

Use this skill when adding new tests to CodyFortran.

This repository uses two test families:

1. `componentTests`: Focused tests for one module/backend/implementation.
2. `simulationTests`: Full end-to-end simulation tests.

## Required References

- [Project README](../../../README.md)
- [Build and test wiring](../../../CMakeLists.txt)
- [Agent context](../../../AGENTS.md)

## Repository Test Layout

### Component tests

Store under `test/componentTests/` and mirror the module hierarchy in folders.

Example:

```text
test/componentTests/SysInteraction/Ylm/Coulomb/BlockEq/
	T_SysInteraction_Ylm_Coulomb_BlockEq.f90
	T_SysInteraction_Ylm_Coulomb_BlockEq.json
```

### Simulation tests

Store under `test/simulationTests/<System>/`.

Example:

```text
test/simulationTests/He1d/
	T_He1d_07_NewScenario.f90
	T_He1d_07_NewScenario.json
```

## Naming Rules

1. Fortran test program file: `T_<Name>.f90`
2. Matching JSON file: `T_<Name>.json`
3. The `.f90` and `.json` base names must be identical.
4. Component tests should include module/backend/variant in `<Name>`.
5. Simulation tests should include system + ordering index + scenario.

## JSON Path Rule (Important)

When loading JSON from Fortran test code, use the full path relative to project root and set `relativeToProjectDirQ_=.true.`.

Component example:

```fortran
jsonFileName = "test/componentTests/SysInteraction/Ylm/Coulomb/BlockEq/T_SysInteraction_Ylm_Coulomb_BlockEq.json"
call Json_LoadJsonFile(manualFileName_=jsonFileName, relativeToProjectDirQ_=.true.)
```

Simulation example:

```fortran
jsonFileName = "test/simulationTests/He1d/T_He1d_07_NewScenario.json"
call Json_LoadJsonFile(manualFileName_=jsonFileName, relativeToProjectDirQ_=.true.)
```

Do not use only `"test/T_...json"` for newly created nested tests.

## Workflow

1. Classify the test type.
- Choose `componentTests` for one module/backend contract.
- Choose `simulationTests` for full fabrication/setup/propagation workflows.

2. Choose a nearby reference test.
- Copy the closest existing test in the same family and adapt minimally.
- Keep the coding style and module import pattern.

3. Create the pair (`.f90` + `.json`) in the correct directory.
- Prefer colocating the pair in the same folder.
- Keep JSON minimal and deterministic.

4. Implement assertions with `testdrive`.
- Use `use testdrive, only: check, error_type`.
- Fail clearly with `if (allocated(error)) error stop "... failure"`.
- Use explicit tolerances for floating-point checks.

5. Keep the test deterministic.
- No random behavior without fixed seed.
- Use small problem sizes and stable thresholds.

6. Build and run tests.

```bash
source setBuildVars.sh -c gnu -t release
cmake -B build && cmake --build build && cmake --install build
cd build && ctest -R T_
```

Run a focused test when iterating:

```bash
cd build && ctest -R T_SysInteraction_Ylm_Coulomb_BlockEq
```

## Build-System Behavior

`CMakeLists.txt` discovers tests automatically via recursive glob of `test/*.f90` and registers each source as a CTest test. Usually no CMake edits are needed when adding a test.

## Test Templates

### Component test template (module/backend)

```fortran
program T_Module_Backend_Variant
	use M_Utils_Types
	use M_Utils_Json
	use testdrive, only: check, error_type
	implicit none

	type(error_type), allocatable :: error
	character(len=:), allocatable :: jsonFileName

	jsonFileName = "test/componentTests/Module/Backend/Variant/T_Module_Backend_Variant.json"
	call Json_LoadJsonFile(manualFileName_=jsonFileName, relativeToProjectDirQ_=.true.)

	! Fabricate + setup only what is required for the component.
	! Compute one or more deterministic observables.
	! Compare against analytic/reference values.

	call check(error, .true.)
	if (allocated(error)) error stop "T_Module_Backend_Variant failure"
end program
```

### Simulation test template (end-to-end)

```fortran
program T_System_XX_Scenario
	use M_Utils_Types
	use M_Utils_Json
	use testdrive, only: check, error_type
	implicit none

	type(error_type), allocatable :: error
	character(len=:), allocatable :: jsonFileName

	jsonFileName = "test/simulationTests/System/T_System_XX_Scenario.json"
	call Json_LoadJsonFile(manualFileName_=jsonFileName, relativeToProjectDirQ_=.true.)

	! Fabricate full simulation chain (grid, operators, method, integrator, propagator).
	! Run a short deterministic propagation or solver loop.
	! Assert final energy/observable against reference.

	call check(error, .true.)
	if (allocated(error)) error stop "T_System_XX_Scenario failure"
end program
```

## Done Checklist

- Correct family chosen: `componentTests` or `simulationTests`.
- Test location matches repository tree conventions.
- `.f90` and `.json` pair created with matching base names.
- JSON load path uses full nested project-relative path.
- Deterministic checks and tolerance-based assertions added.
- `ctest -R T_` (and focused test) passes.
