# CodyFortran: Correlated Dynamics with Fortran

![CodyFortran Logo](assets/img_850x248.png)

![Language](https://img.shields.io/badge/fortran-2008-brightgreen.svg?style=flat-square)
![License](https://img.shields.io/badge/license-BSD--3--Clause-blue.svg?style=flat-square)

## Abstract

**CodyFortran** is a Modern Fortran framework for simulating interacting quantum many-body systems. Designed to be flexible and versatile, it enables advanced simulations of fermionic, bosonic, and mixed quantum systems with customizable potential and interaction functions. The framework supports time propagation and ground-state calculations, suitable for lattice and continuous systems across one, two, and three dimensions.

## Table of Contents
- [Installation](#installation)
  - [Prerequisites](#prerequisites)
  - [Building Third-Party Libraries Locally](#building-third-party-libraries-locally)
  - [Building CodyFortran with CMake](#building-codyfortran-with-cmake)
- [Usage](#usage)
- [Example Programs](#example-programs)
- [Capabilities](#capabilities)
  - [Simulation Methods](#simulation-methods)
- [The CodyFortran Architecture: A Deep Dive into Modularity](#the-codyfortran-architecture-a-deep-dive-into-modularity)
  - [1. Mirrored Abstraction Hierarchies](#1-mirrored-abstraction-hierarchies)
  - [2. Separation of Interface (M_) and Implementation (S_)](#2-separation-of-interface-m_-and-implementation-s_)
  - [3. Fabrication and Procedure Pointers](#3-fabrication-and-procedure-pointers)
    - [The Special Case: `*List` Modules](#the-special-case--list--modules)
- [Setting Up a New Simulation: A Step-by-Step Guide](#setting-up-a-new-simulation-a-step-by-step-guide)
  - [Step 1: Create the Program File and Load Configuration](#step-1-create-the-program-file-and-load-configuration)
  - [Step 2: Import Required Module Interfaces](#step-2-import-required-module-interfaces)
  - [Step 3: Initialize Modules with `Fabricate` Calls](#step-3-initialize-modules-with-fabricate-calls)
  - [Step 4: Allocate Resources with `Setup` Calls](#step-4-allocate-resources-with-setup-calls)
  - [Step 5 & 6: Implement Logic and Output](#step-5--6-implement-logic-and-output)
- [Extending CodyFortran](#extending-codyfortran)
- [Contributing](#contributing)
- [Coding Conventions](#coding-conventions)
  - [Naming Conventions](#naming-conventions)
  - [Function Signatures](#function-signatures)
  - [Boolean Operators](#boolean-operators)
- [Abbreviations](#abbreviations)
- [License](#license)

## Installation

### Prerequisites

CodyFortran depends on several specialized libraries:

- **EXPOKIT**: Matrix exponential routines
- **FFTW3**: Fast Fourier Transform library
- **SHTNS**: Spherical Harmonic Transform library
- **Fortran stdlib**: Fortran standard library
- **Fortran test-drive**: Testing framework for Fortran
- **BLAS**: Basic Linear Algebra Subprograms
- **LAPACK**: Linear Algebra Package
- **JSON-Fortran**: JSON parsing capabilities for Fortran
- **gsl**: GNU Scientific Library (optional, for additional numerical methods)
- **arpack**: ARPACK library for large-scale eigenvalue problems

These dependencies can be provided in two ways:

1. **System-provided libraries**: Using your system's package manager (if available)
2. **Locally compiled libraries**: Built in the project's `thirdParty` folder

> **Important**: Fortran libraries that use modules must be compiled with the same Fortran compiler you plan to use for CodyFortran. Different compilers may produce incompatible binary interfaces, which can cause linking errors.

### Building Third-Party Libraries Locally

To build the required libraries locally:

1. Configure which libraries to download by editing `thirdParty/getSources.sh`:

```bash
# Set to "true" to download the library, "false" to skip
DOWNLOAD_EXPOKIT="true"
DOWNLOAD_SHTNS="true"
DOWNLOAD_STDLIB="true"
DOWNLOAD_TEST_DRIVE="true"
DOWNLOAD_JSON_FORTRAN="true"
DOWNLOAD_GSL="true"
DOWNLOAD_ARPACK="true"
# these are not needed when building with intel compilers and MKL
DOWNLOAD_FFTW3="true"
DOWNLOAD_OPENBLAS="true"
```

2. Build the libraries:

```bash
# Clone repository
git clone https://github.com/flackner/codyfortran
cd CodyFortran
cd thirdParty

# source the setBuildVars.sh file (the one in thirdParty not the top-level one)
# flags:
# -h help: display help message
# -c compiler family: gnu or intel
# -t build type: debug or release
# For third-party libraries, release build type is recommended unless
# you need to debug those libraries
source setBuildVars.sh -c gnu -t release

# Download the selected libraries (delete the old folders fetching updates)
bash getSources.sh

# Build and install the libraries
cmake -B build
cmake --build build
cmake --install build
```

The libraries will be installed in the `thirdParty/installation/` directory.  When building CodyFortran, it will first check for libraries in the `thirdParty/installation/lib/` and `thirdParty/installation/lib64/` directories. If a required library isn't found there, it will search in the standard system paths.

### Building CodyFortran with CMake

After building the third-party libraries, it is recommended to start a new terminal session before sourcing `setBuildVars.sh` in the main project directory.

```bash
# Clone repository (can be skipped if already done)
git clone https://github.com/flackner/codyfortran
cd CodyFortran

# source the setBuildVars.sh file
# flags:
# -h help: display help message
# -c compiler family: gnu or intel
# -t build type: debug or release
source setBuildVars.sh -c gnu -t release

# Build and install CodyFortran
cmake -B build
cmake --build build
cmake --install build

# Run the test suite
cd build
ctest -R T_
```

CodyFortran's console output uses ANSI color escape codes and a few non-ASCII characters. If you prefer plain ASCII (for example when redirecting output to a file) or you don't like the styled output, disable the fancy printing by setting
```bash
export CODY_PRETTY_PRINT="OFF"
```
Per default, `CODY_PRETTY_PRINT` is set to `"ON"` in `setBuildVars.sh`. You can change this default behavior by editing that script.

## Usage

This section shows how to run provided executables and example programs, and how their JSON inputs are paired.

All `.f90` files in the `app` folder are automatically compiled into executable files. When running an executable, you can provide the JSON configuration input file in two ways:

1. Explicit JSON file specification:

```bash
# Set the dynamic library paths (needs to be done once per terminal session)
source setBuildVars.sh -c gnu -t release

# Run the executable with the specified input file
./myProg.exe myInput.json

# Example: Run the executable with input file myProg.json in current directory
./myProg.exe
```
## Example Programs

The \`test/\` folder contains several runnable programs together with a matching JSON configuration that shares the same base name. For example:

- `T_He1d_03_ImagTimePropagationTdhx.f90` ↔ `T_He1d_03_ImagTimePropagationTdhx.json`
- `T_FermiHubbard_04_ImagTimePropagationTdci.f90` ↔ `T_FermiHubbard_04_ImagTimePropagationTdci.json`

These examples serve as excellent starting points for understanding how to set up different types of simulations using CodyFortran.

## Capabilities

-   **Grids**
    -   Lattice grids in 1D/2D/3D with periodic/hard boundaries.
    -   Continuous grids in 1D: equidistant and FEDVR.
    -   Continuous grids in 2D: polar and cartesian.
    -   Continuous grids in 3D: spherical (Ylm) and cartesian with radial FEDVR.
-   **Physics Operators**
    -   Kinetic: different hopping schemes for lattices; laplacian for continuous.
    -   External potentials: Many presets (Coulomb, harmonic) and easy to extend.
    -   Interactions: Many presets (lattice on-site, Coulomb) and easy to extend.
-   **Time Evolution and Solvers**
    -   Real-time and imaginary-time propagation.
    -   Various Runge–Kutta/Crank-Nicholson integrators and interface to gsl odeiv2.
    -   Split-step propagation schemes using several integrators for different Hamiltonian parts.
    -   Ground-state search via ED, imaginary time propagation or iterative scf methods.
-   **I/O, Analysis, and Utilities**
    -   JSON-driven configuration mirroring the code architecture.
    -   Data output for densities, observables, spectra, and RDMs
    -   Many test programs under \`test/\` serve as runnable examples with matching JSON.

### Simulation Methods

CodyFortran is built to support a wide array of simulation methods, each has a different definition of the state representation of the physical system. The implementation of each method is located as a subdirectory inside the \`src/Method/\` folder. The choice of method is configured in your JSON input file. The table below details the primary methods available.

| Category | Method | State Representation | JSON Configuration Key |
| :--- | :--- | :--- | :--- |
| **Single-Body** | Single Particle | Single-particle wavefunction on the grid. | `"sb": {}` |
| **Grid-Based Many-Body** | Full Expansion on Grid | Full wavefunction on a direct-product grid basis. | `"gridBased": {"full": {}}` |
| **Geminal-Based Many-Body** | Geminal Expansion | Wavefunction expanded in geminals (two-particle functions). | `"gemBased": {}` |
| **Orbital-Based Many-Body** | **TDCI** (Time-Dependent CI) | `CI Coefficients (Dynamic)` + `Orbitals (Fixed)` | `"orbBased": {"tdci": {}}` |
| | **MCTDHX** | `CI Coefficients (Dynamic)` + `Orbitals (Dynamic)` | `"orbBased": {"mctdhx": {}}` |
| | **TDHX** (Time-Dependent Hartree-Fock) | Single Slater determinant (\`Orbitals are Dynamic\`). | \`"orbBased": {"tdhx": {}}\` |


## The CodyFortran Architecture: A Deep Dive into Modularity

CodyFortran is a modular, extensible framework for performing quantum dynamics simulations. Its architecture is designed for maximum flexibility and clarity, allowing users to easily swap physical models, numerical methods, and computational backends by configuring a simple JSON file.

The power and extensibility of CodyFortran stem from three core architectural principles:

1.  **Mirrored Abstraction Hierarchies**: The abstraction hierarchy of the code is perfectly mirrored in both the source folder structure and the JSON configuration file structure. This creates an intuitive and predictable system where the configuration is a direct map of the code's design.
2.  **Strict Separation of Interface and Implementation**: Each module's public API (the interface) is defined in a single `M_*.f90` file, completely separate from its various concrete implementations, which reside in `S_*.f90` files.
3.  **API via Procedure Pointers and Fabrication**: A module's API functions are defined as `deferred` procedure pointers. These pointers are assigned to concrete implementations during a "Fabrication" step, where a chain of `*_Fabricate` calls traverses the abstraction hierarchy to wire the final routines.

### 1. Mirrored Abstraction Hierarchies

CodyFortran is organized around a clear abstraction hierarchy. The most general concepts form the "base modules" at the top level of the `src/` directory (e.g., `Method`), which are made concrete by "child modules" in subdirectories. This hierarchy is strictly mirrored in two places:

*   **The Folder Structure**: The nesting of folders within `src/` represents specialization.
*   **The JSON Configuration**: The nesting of objects within the JSON file represents the same specialization.

This creates a one-to-one mapping that makes the system easy to navigate. The relationship between a base module and a child module is always an **"is-a"** relationship (like inheritance), not a "has-a" relationship (like composition).

| Abstraction Hierarchy Path (`src/`) | Mirrored JSON Structure | Interpretation |
| :--- | :--- | :--- |
| `Method/Mb/OrbBased/Mctdhx/` | `{"method": {"mb": {"orbBased": {"mctdhx": {}}}}}`| The method **is a** Mb method, which **is an** OrbBased method, which **is an** Mctdhx method. |

Because of this "is-a" logic, each base module object in the JSON file must contain **exactly one** child module sub-object to define its concrete type.

### 2. Separation of Interface (`M_`) and Implementation (`S_`)

Every module consists of two distinct parts:

*   **The Interface (`M_*.f90`)**: This is the stable, public face of the module. It defines the "contract" for what the module can do by declaring public data, types, and a set of `deferred` procedure pointers for its API functions (e.g., `Apply`, `Setup`). User code should **only ever** interact with this `M_` interface.
*   **The Implementation (`S_*.f90`)**: These files contain the concrete algorithms that fulfill the interface's contract. There can be many different implementations for a single interface, and the choice of which one to use is made at runtime based on the JSON configuration.

This separation ensures that you can change or add new implementations without ever modifying the code that uses the module.

### 3. Fabrication and Procedure Pointers

CodyFortran's API is not built on static function calls but on dynamically assigned procedure pointers. This dynamic binding happens during the **Fabrication** phase.

The process works as a chain of command down the abstraction hierarchy:

1.  Your program calls the `*_Fabricate` routine of a top-level base module (e.g., `Method_Fabricate`).
2.  This routine reads the JSON configuration to identify which child module was selected (e.g., `mb`).
3.  It then calls the `*_Fabricate` routine *of that child module*.
4.  This process repeats, traversing deeper into the nested folders and corresponding JSON objects until the most concrete implementation is reached.
5.  The final, most concrete `*_Fabricate` routine assigns its specific functions to the procedure pointers defined in the top-level `M_` interface file.

This cascading fabrication process is how the framework wires the entire simulation together, connecting the abstract API calls in your program to the concrete scientific code chosen in your JSON file.

### The Special Case: `*List` Modules

Modules whose names end in `List` (e.g., `M_IntegratorList`, `M_DiagonalizerList`) are containers for multiple objects. They are used for orchestration. For example, a split-step propagator needs two distinct integrators.

In the JSON, you configure multiple child objects by appending a unique key to each.

```json
{
  "integratorList": {
    "rk1": {
      "o2Impl": { }
    },
    "rk2": {
      "o4Expl": { }
    }
  }
}
```

The program then uses these named instances as required, for instance, by wiring `rk1` to the first propagator sub-step and `rk2` to the second.


## Setting Up a New Simulation: A Step-by-Step Guide

This guide walks through the creation of a new simulation program, keeping the architectural principles in mind.

### Step 1: Create the Program File and Load Configuration

Start with a new Fortran program (`.f90`) in the `app/` directory. The first step is to load the JSON configuration file, which will dictate the entire simulation's structure.

```fortran
program P_MySimulation
  use M_Utils_Json
  ! ... other modules
  call Json_LoadJsonFile
end program
```

### Step 2: Import Required Module Interfaces

Import the `M_*` interface modules needed for your simulation. You will only ever interact with these public interfaces.

```fortran
use M_Grid
use M_SysKinetic
use M_Method
use M_IntegratorList
use M_Propagator
```

### Step 3: Initialize Modules with `Fabricate` Calls

Call the `*_Fabricate` routine for each module. This phase reads the JSON config file, navigates the configuration hierarchy, and wires the procedure pointers to the selected `S_*` implementations.

```fortran
call Grid_Fabricate
call SysKinetic_Fabricate
call Method_Fabricate

! For a list module, the fabricate call might take programmatic input
! to specify how the list items are used.
type(T_IntegratorList_FabricateInput) :: IntegratorListInput(2)
IntegratorListInput(1) % TimeDerivative => Method_Mb_OrbBased_TimeDerivativeOrbsLin
IntegratorListInput(2) % TimeDerivative => Method_Mb_OrbBased_TimeDerivativeCoeffsPlusOrbsNonLin
call IntegratorList_Fabricate(IntegratorListInput)

call Propagator_Fabricate()
```

### Step 4: Allocate Resources with `Setup` Calls

Call the `*_Setup` routines. This is the phase where memory is allocated and operators are constructed.

```fortran
call Grid_Setup
call SysKinetic_Setup
call Method_Setup
call Propagator_Setup
```

### Step 5 & 6: Implement Logic and Output

With all modules fabricated and set up, you can now call the interface procedures (e.g., `Propagator_Propagate`) to run the simulation and generate output.

```fortran
! The call to Propagator_Propagate is simple, but behind the scenes it uses the
! integrators from the IntegratorList, which in turn call the time derivatives
! provided by the Method module—all wired up during Fabricate.
call Propagator_Propagate(Method_state, t, t + timeStep)
```

## Extending CodyFortran

This architecture makes adding new functionality straightforward. To add a new 1D potential:
1.  **Create a new Subfolder**: Create a new folder `src/SysPotential/Linear/YourPotential` with two files inside: `M_SysPotential_Linear_YourPotential.f90` and `S_SysPotential_Linear_YourPotential.f90`. Most of the time, you can copy an existing potential implementation and modify it.
2.  **Create an Implementation Submodule**: Inside `S_SysPotential_Linear_YourPotential.f90`, implement your logic and a `*_Fabricate` routine that reads JSON keys (e.g., `sysPotential.linear.yourPotential.strength`) and assigns the `SysPotential_Apply` procedure pointer.
3.  **Register the New Type**: In `src/SysPotential/Linear/S_SysPotential_Linear.f90`, add an `if` block to its `*_Fabricate` routine to detect your new JSON object key and call your submodule's `*_Fabricate`.
4.  **Use it in JSON**: You can now select your new potential with the structure:
    ```json
    "sysPotential": {
      "linear": {
        "yourPotential": { 
          "strength": 1.0 
          }
      }
    }
    ```

This pattern of extending CodyFortran by creating a new folder, creating the implementation in the S_*.f90 file, registering it in the parent dispatcher, and configuring it via JSON applies universally across the framework.

## Contributing
This is a small project, and there is no established community. Contributions are welcome, whether they are bug fixes, improvements, or extensions. If you are interested in contributing, please feel free to open an issue or a pull request. All contributions will be considered.

## Coding Conventions

To maintain consistency and readability, CodyFortran follows specific coding conventions.

### Naming Conventions

Since Fortran lacks namespaces, the underscore `_` is used as a namespace separator.

*   **Functions, Subroutines, and Modules**: Use **UpperCamelCase**.
    *   Example: `Method_Fabricate`, `Grid_InnerProduct`.
    *   Namespaced: `Namespace_SubNamespace_FunctionName`.
*   **Variables**: Use **lowerCamelCase**.
    *   Example: `timeStep`, `methodState`, `gridPoints`.
    *   Namespaced: `Namespace_SubNamespace_variableName` (e.g., `Grid_Linear_xmin`).

### Function Signatures

Argument order in procedures follows a strict pattern:
1.  **Output** arguments (`intent(out)`)
2.  **Input/Output** arguments (`intent(inout)`)
3.  **Input** arguments (`intent(in)`)
4.  **Optional Output** arguments
5.  **Optional Input/Output** arguments
6.  **Optional Input** arguments

**Optional Arguments**: All optional arguments must have a trailing underscore `_` and must be called using keyword arguments.

```fortran
subroutine ExampleRoutine(fOut, fInOut, fIn, optOut_, optInOut_, optIn_)
  real(R64), intent(out) :: fOut
  real(R64), intent(inout) :: fInOut
  real(R64), intent(in)  :: fIn
  real(R64), intent(out), optional :: optOut_
  real(R64), intent(inout), optional :: optInOut_
  real(R64), intent(in),  optional :: optIn_
  ! ...
end subroutine
```

### Boolean Operators

Prefer the use of symbolic comparison operators for inequalities and dot-operators for equality:
*   Use `<`, `<=`, `>`, `>=`
*   Use `.eq.`, `.ne.`

## Abbreviations
The codebase and JSON configuration use a number of recurring abbreviations. This list summarizes the most common ones for quick reference.

| Abbreviation | Explanation |
| :--- | :--- |
| `sb` | **Single-Body**. Method for a single particle where the wavefunction is represented directly on a grid. |
| `mb` | **Many-Body**. General category for methods simulating interacting particles (includes orbital-, geminal-, and grid-based approaches). |
| `gridBased` | **Grid-Based Expansion**. Many-body method where the full wavefunction is expanded on a direct-product grid basis. |
| `gemBased` | **Geminal-Based Expansion**. Many-body method where the wavefunction is expanded in geminals (two-particle functions). |
| `orbBased` | **Orbital-Based Expansion**. Many-body method where the wavefunction is expanded using single-particle orbitals. |
| `tdci` | **Time-Dependent Configuration Interaction**. Orbital-based method where CI coefficients evolve in time while orbitals remain fixed. |
| `mctdhx` | **Multi-Configuration Time-Dependent Hartree-Fock for Bosons, Fermions, and mixtures.**. Orbital-based method where both CI coefficients and orbitals evolve in time. |
| `tdhx` | **Time-Dependent Hartree-Fock for Bosons, Fermions, and mixtures**. Mean-field approximation where the wavefunction is a single Slater determinant (or permanent) with time-evolving orbitals. |
| `ci` | **Configuration Interaction**. Expansion of the many-body wavefunction in a basis of Slater determinants or permanents. |
| `rdm` | **Reduced Density Matrix**. Matrices representing the quantum state of a subsystem (1-body, 2-body, etc.) by tracing out other degrees of freedom. |
| `fedvr` | **Finite-Element Discrete Variable Representation**. A high-order grid discretization method combining finite elements with local DVR bases. |
| `ecs` | **Exterior Complex Scaling**. A method using a complex-valued coordinate contour to absorb outgoing wavepackets at boundaries. |
| `ylm` | **Spherical Harmonics**. Angular basis functions $Y_l^m$ used for representing wavefunctions on spherical grids. |
| `rk` | **Runge–Kutta**. A family of iterative methods for the approximate numerical solution of ordinary differential equations (ODEs). |
| `cn` | **Crank–Nicolson**. An implicit, unitary time-stepping scheme often used for the Schrödinger equation. |
| `sil` | **Short Iterative Lanczos**. A Krylov subspace time-propagation method based on the Lanczos algorithm. |
| `json` | **JavaScript Object Notation**. The human-readable data format used for CodyFortran configuration files. |
| `o4Expl` | **Explicit 4th order**. Explicit integrator of order 4. |
| `o2Impl` | **Implicit 2nd order**. Implicit integrator of order 2. |
| `fabricate` | **Fabrication**. The initialization phase where abstract procedure pointers are bound to concrete implementations based on JSON input. |
| `setup` | **Setup**. The initialization phase following fabrication where memory is allocated and operators are constructed. |
| `M_` | **Module Interface**. Prefix for files defining public interfaces and deferred procedure pointers (e.g., `M_Grid.f90`). |
| `S_` | **Submodule Implementation**. Prefix for files containing concrete implementations of interfaces (e.g., `S_Grid_Linear_Fedvr.f90`). |
| `ctx` | **Context**. A derived type holding the configuration and state for a specific subsystem or module. |
| `e` | **Element**. Generic name for an element within a container or list structure (e.g., in `*List` modules). |
| `evecs` / `evals` | **Eigenvectors / Eigenvalues**. Common variable names for the results of a diagonalization procedure. |
| `bt` | **Body Type**. Index identifying a distinguishable particle species. |
| `*Q` | **Query / Question**. Suffix convention for logical (boolean) variables (e.g., `toScreenQ` means "write to screen?"). |

## License

Distributed under the BSD 3-Clause License. See `LICENSE` for more information.
