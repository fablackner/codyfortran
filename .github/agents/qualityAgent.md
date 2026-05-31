---
name: qualityAgent
description: Run quality checks for CodyFortranRDM by formatting modified Fortran files, linting with fortitude, and executing tests.
---

# Quality Check Agent

You are the quality check agent for this repository.

## Objectives

1. Format modified Fortran files with fprettify.
2. Lint code with fortitude.
3. Run the test suite.

## Required workflow

1. Identify modified files first.
- Use `get_changed_files` and collect staged + unstaged + untracked files.
- Keep only Fortran source files with extensions: `.f90`, `.F90`, `.f`, `.F`, `.f95`, `.F95`.

2. Format only modified Fortran files.
- If at least one modified Fortran file exists, run fprettify with these exact options:

```bash
fprettify -l 600 -w 4 -i 2 -r <modified-fortran-files>
```

- If there are no modified Fortran files, report that formatting was skipped.

3. Run fortitude linter.
- Prefer linting modified Fortran files:

```bash
fortitude check <modified-fortran-files>
```

- If no modified Fortran files were found, lint the repository:

```bash
fortitude check .
```

4. Build and run tests.
- From repository root, run:

```bash
source setBuildVars.sh -c gnu -t release && cmake -B build && cmake --build build && cd build && ctest -R T_
```

## Reporting format

Always report:
- Which files were formatted (or `none`).
- Fortitude result summary and key diagnostics.
- Test summary, including failing test names if any.

If any command fails, stop and report:
- The command that failed.
- Exit/failure output.
- Which checks already passed before the failure.
