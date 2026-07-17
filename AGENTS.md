# CodyFortranRDM Agent Context

## WHAT
CodyFortranRDM is a modular Fortran simulation tool for quantum many-body systems. It can compute ground-state properties and time-dependent properties using multiple many-body methods. Methods differ in the representation of the many-body quantum state and/or in the equation of motion. The JSON configuration mirrors the code hierarchy, and a chain of `*_Fabricate` calls binds deferred procedure pointers to concrete implementations at runtime.

## WHY
The goal is to combine different methods in one modular framework while reusing shared components such as Grid, Coeffs, Diagonalizer, Propagator, SysInteraction, SysKinetic, SysPotential, and Utils. Modularization is enabled by small module interfaces with clear interface contracts.

## Coding Guidelines

- Use `_` as namespace separator because Fortran lacks namespaces.
- Use **UpperCamelCase** for functions, subroutines, and modules (for example: `Method_Fabricate`, `Grid_InnerProduct`).
- Use **lowerCamelCase** for variables; keep namespaced variables as `Namespace_SubNamespace_variableName`.
- Keep procedure argument order strict: `intent(out)`, `intent(inout)`, `intent(in)`, then optional variants in the same order.
- Suffix all optional arguments with `_` and call them with keyword arguments.
- Prefer symbolic inequality operators: `<`, `<=`, `>`, `>=`.
- Prefer `.eq.` and `.ne.` for equality and inequality checks.

## Architecture Guidelines
- **Interface Modules** (`M_*.f90`): Define contracts via deferred procedure pointers. Rule: Only import `M_*` files.
- **Implementation Submodules** (`S_*.f90`): Implement concrete algorithms.
- **Fabrication**: `*_Fabricate()` reads JSON and binds procedure pointers at runtime.

## Module Documentation Index (A_*.md)

### Required Delegation Workflow

The following modules are currently available. If programming work requires understanding one or more modules, load the dedicated module doc file using the `read_file` tool.

When work targets a specific module, follow this flow:

1. Identify the target module in the index below.
2. Load the corresponding `A_*.md` file before implementing changes.
3. Prefer a dedicated module specialist subagent when the change is non-trivial.

### Unified Module Index and Reference

**Grid**
Domain: Spatial & Grid Foundation.
Purpose: Spatial discretization, metric-aware inner products, orthonormalization, subspace projection.
Key interfaces/notes: Backends include Linear, Square, Polar, Spherical, Ylm, and Lattice. Key exports include `Grid_nPoints`, `Grid_InnerProduct`, `Grid_Orthonormalize`, and `Grid_ProjectOnSubspace`.
Module doc: `src/Grid/A_Grid.md`.

**SysKinetic**
Domain: Hamiltonian Operators.
Purpose: Applies kinetic operator.
Key interfaces/notes: `SysKinetic_MultiplyWithKineticOp(dOrb, orb, t, bt)`; supports Linear finite-difference/FFT, Lattice hopping, and Ylm radial+centrifugal forms.
Module doc: `src/SysKinetic/A_SysKinetic.md`.

**SysPotential**
Domain: Hamiltonian Operators.
Purpose: Builds/applies external one-body potential.
Key interfaces/notes: `SysPotential_FillExternalPotential` and `MultiplyWithExternalPotential`; supports time/body-type flags and multiple backend models.
Module doc: `src/SysPotential/A_SysPotential.md`.

**SysInteraction**
Domain: Hamiltonian Operators.
Purpose: Builds/applies two-body interaction potential.
Key interfaces/notes: `FillInteractionSrc -> FillInteractionPotential -> MultiplyWithInteractionPotential`; includes Linear SoftYukawa, Lattice Hubbard, and Ylm Coulomb-Poisson variants.
Module doc: `src/SysInteraction/A_SysInteraction.md`.

**SysGauge**
Domain: Hamiltonian Operators.
Purpose: Applies the velocity-gauge laser coupling A(t)·p/m + A²/(2m) (non-diagonal, hence not a SysPotential).
Key interfaces/notes: `SysGauge_MultiplyWithGaugeOp(dOrb, orb, time, bt_)` (SysKinetic-style operator contract); Ylm backend with sin²-envelope vector potential and Δl = ±1 channel coupling on Fedvr/FedvrEcs radial grids; consumed by Method via `Method_Mb_OrbBased_ApplyGaugeOp`.
Module doc: `src/SysGauge/A_SysGauge.md`.

**Orbs**
Domain: Quantum State Containers.
Purpose: Orbital coefficient container and orbital-space operations.
Key interfaces/notes: `Orbs_orbs(nBasis, nOrbs)`, `Orbs_Orthonormalize`, `Orbs_ProjectOnSubspace`, `Orbs_SaveOrbs`; supports restricted/unrestricted modes.
Module doc: `src/Orbs/A_Orbs.md`.

**OrbsInit**
Domain: Quantum State Containers.
Purpose: Initializes orbitals for selected grid/physics model.
Key interfaces/notes: `OrbsInit_Initialize(orbs)`; supports HO, hydrogen-like, on-site, load-from-file, and grid-point initialization.
Module doc: `src/OrbsInit/A_OrbsInit.md`.

**Coeffs**
Domain: Quantum State Containers.
Purpose: CI coefficients, Hamiltonian application, and RDM filling.
Key interfaces/notes: `Coeffs_ApplyH1FillRdm1`, `Coeffs_ApplyH2FillRdm2`, `Coeffs_ApplyExcitation`, and `Coeffs_FillRdm{1,2,3}Bt`; includes Generic and Hubbard variants.
Module doc: `src/Coeffs/A_Coeffs.md`.

**CoeffsInit**
Domain: Quantum State Containers.
Purpose: Initializes CI coefficient vectors.
Key interfaces/notes: `CoeffsInit_Initialize(coeffs)`; supports unary ground state, excited-state operator chains, and load from `coeffs.in`.
Module doc: `src/CoeffsInit/A_CoeffsInit.md`.

**TwoRdm**
Domain: Quantum State Containers.
Purpose: Two-body RDM state representation and algebraic operations.
Key interfaces/notes: `TwoRdm_ApplyCommutatorH1/H2`, `FillRdm1/2FromTwoRdmState`, and `FillTwoRdmStateFromRdm2`; includes Generic, Singlet, and SingletTr0 forms.
Module doc: `src/TwoRdm/A_TwoRdm.md`.

**TwoRdmInit**
Domain: Quantum State Containers.
Purpose: Initializes two-body RDM states.
Key interfaces/notes: Supports excited-state and load-from-file initialization using `rdm2Bt{bt1}_{bt2}.in`.
Module doc: `src/TwoRdmInit/A_TwoRdmInit.md`.

**ConfigList**
Domain: Many-Body Basis.
Purpose: Fock-basis enumeration and excitation connectivity.
Key interfaces/notes: `ExciteConfiguration(iCNew, factor, creates, destroys, iC)`; includes precomputed singles/doubles and supports fermionic and bosonic cases.
Module doc: `src/ConfigList/A_ConfigList.md`.

**Method**
Domain: Dynamics & Propagation.
Purpose: Unified method time-derivative interface for propagators.
Key interfaces/notes: `Method_TimeDerivative(dState, state, t)`; includes Sb, TDHX, TDCI, MCTDHX, and TD-2RDM variants.
Module doc: `src/Method/A_Method.md`.

**Propagator**
Domain: Dynamics & Propagation.
Purpose: Time-evolution facade for state propagation.
Key interfaces/notes: `Propagator_Propagate(state, t0, t1)`; supports single integrator, split-step (order2/order4), and eigen-expansion modes.
Module doc: `src/Propagator/A_Propagator.md`.

**IntegratorList**
Domain: Dynamics & Propagation.
Purpose: Pluggable ODE integrators.
Key interfaces/notes: `Integrate(state, t0, t1)` with callback; includes RK variants, CrankNicolson, Expokit, SIL, and GSL ODE backends.
Module doc: `src/IntegratorList/A_IntegratorList.md`.

**GroundSolver**
Domain: Dynamics & Propagation.
Purpose: SCF ground-state solver facade.
Key interfaces/notes: `GroundSolver_Approach(state, alpha, t)`; supports Hartree-Fock assembly, diagonalization, gauge alignment (`Orbs_AlignOnReference`), mixing (via Mixing), and orthonormalization.
Module doc: `src/GroundSolver/A_GroundSolver.md`.

**Mixing**
Domain: Numerical Infrastructure.
Purpose: Self-consistency mixers (damping/acceleration) for fixed-point iterations.
Key interfaces/notes: `Mixing_SetupMixer(mixer, dim)`, `Mixing_Mix(mixer, x, xNew, alpha)` on generic complex vectors; Linear and DIIS (Pulay) backends; consumed by GroundSolver (orbital or Hartree-potential/density targets).
Module doc: `src/Mixing/A_Mixing.md`.

**Reconstruction**
Domain: RDM Reconstruction & Purification.
Purpose: BBGKY closure for approximating 3-RDM effects in 2-RDM evolution.
Key interfaces/notes: `Reconstruction_AddCollisionOp(dD2, D2, h2)`; includes Valdemoro, Nakatsuji-Yasuda, and NN-based closures.
Module doc: `src/Reconstruction/A_Reconstruction.md`.

**Purification**
Domain: RDM Reconstruction & Purification.
Purpose: Enforces N-representability constraints on 2-RDM states.
Key interfaces/notes: `Purification_PurifyTwoRdmState(D2, h2)`; iterative D/Q positivity correction with optional energy projection.
Module doc: `src/Purification/A_Purification.md`.

**DiagonalizerList**
Domain: Numerical Infrastructure.
Purpose: Unified dense/iterative eigensolver facade.
Key interfaces/notes: `Diagonalize(t, evecsQ)` with `ApplyMatOnVec`; supports LAPACK and ARPACK modes and outputs evals/evecs/nFound.
Module doc: `src/DiagonalizerList/A_DiagonalizerList.md`.

**Absorber**
Domain: Numerical Infrastructure.
Purpose: Boundary damping masks to suppress edge reflections.
Key interfaces/notes: `Absorber_ApplyAbsorber(orbs)`; uses linear cosinus masks near boundaries.
Module doc: `src/Absorber/A_Absorber.md`.

**Utils**
Domain: Numerical Infrastructure.
Purpose: Stateless numerical/algebra utility modules.
Key interfaces/notes: Includes types, BLAS/LAPACK/ARPACK/FFT/GSL wrappers, RDM algebra utilities, derivatives/convolutions, JSON/data I/O, and logging.
Module doc: `src/Utils/A_Utils.md`.

---

## Documentation Pointers
- **Build:** `README.md`, `CMakeLists.txt`
- **Tests:** `test/componentTests/<Module>/T_*.f90` (unit-level) and `test/simulationTests/<System>/T_*.f90` paired with `T_*.json` configs (end-to-end)
- **Add component:** Mirror directory structure, e.g., `src/SysPotential/.../MyPotential/`
