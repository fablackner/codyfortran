! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Data and operator interfaces for orbital-based many-body methods.
!>
!> This module defines shared data structures (orbital partitioning, reduced
!> density matrices, and Hamiltonian tensors) and declares abstract interfaces
!> for building and applying one- and two-body operators in an orbital basis.
!> Concrete implementations are bound at runtime by
!> `Method_Mb_OrbBased_Fabricate`.
!>
!> S_Method_Mb_OrbBased_FieldFreeMeanField (the frozen mean-field reference
!> U) is an add-on of the Cayley variant of the exact field-free part, not a
!> third exact-part variant next to FieldFreeEigen and FieldFreeCayley.
module M_Method_Mb_OrbBased
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Binds orbital-based procedures at runtime and initializes related data
    !> according to the selected configuration.
    module subroutine Method_Mb_OrbBased_Fabricate
    end subroutine

    !> Setup shared by the orbital-based methods, called first by their
    !> Method_Setup (after the Grid and Sys* setups and Orbs_Fabricate): builds
    !> the skeleton of the definite-m node path if it is bound and the
    !> exact field-free part (Method_Mb_OrbBased_SetupFieldFree) if the method
    !> binds Method_PropagateExactPart (error stop if that setup is not bound,
    !> e.g. Method_PropagateExactPart bound by other code).
    module subroutine Method_Mb_OrbBased_Setup
    end subroutine

    !> Whether the one-body Hamiltonian of Method_Mb_OrbBased_ApplySingleBodyOp
    !> (and of Method_Mb_OrbBased_FillH1) is time-independent: the kinetic
    !> operator, the external potential and, if bound (a sysGauge block
    !> exists), the gauge coupling all are (their *_timeIndependentQ flags,
    !> read at call time). The methods with static orbitals (TDCI,
    !> TD-2RDM-StaticOrbs) refill h1 at every time derivative otherwise.
    module function Method_Mb_OrbBased_H1TimeIndependentQ() result(res)
      logical :: res
    end function

    !> Eigen variant of Method_Mb_OrbBased_SetupFieldFree (Ylm grids): builds
    !> and diagonalizes the blocks of the field-free one-body Hamiltonian
    !> h0 = T + V_ext per l and checks the requirements of the exact part
    !> (S_Method_Mb_OrbBased_FieldFreeEigen).
    module subroutine Method_Mb_OrbBased_SetupFieldFreeEigen
    end subroutine

    !> Cayley variant of Method_Mb_OrbBased_SetupFieldFree (grids whose
    !> kinetic backend binds SysKinetic_SolveShifted, e.g. grid.fem): checks
    !> that h0 = T + V_static is static, caches V_static, reads the GMRES
    !> controls and ends with Method_Mb_OrbBased_SetupFieldFreeMeanField
    !> (S_Method_Mb_OrbBased_FieldFreeCayley).
    module subroutine Method_Mb_OrbBased_SetupFieldFreeCayley
    end subroutine

    !> Fills the spatial 1-RDM in the orbitals of the state: rdm1 itself for
    !> unrestricted or single-species (nOrbsSum = nOrbsInState), the sum of
    !> the two spin blocks for restricted (nOrbsSum = 2 nOrbsInState).
    module subroutine Method_Mb_OrbBased_FillSpatialRdm1(rdm1Spatial, rdm1)
      !> Output: spatial 1-RDM (nOrbsInState x nOrbsInState).
      complex(R64), intent(out), allocatable :: rdm1Spatial(:, :)
      !> One-body reduced density matrix (nOrbsSum x nOrbsSum).
      complex(R64), intent(in), contiguous   :: rdm1(:, :)
    end subroutine

    !> Fills the matrix h1(i, j) = <orbs_i|hOrbs_j> of a one-body operator
    !> from its action hOrbs on the stored orbitals, for orbitals of equal
    !> body type (zero otherwise), in the spin-orbital layout of
    !> Method_Mb_OrbBased_FillH1: nOrbsSum x nOrbsSum, restricted storage:
    !> the spatial block in both spin blocks (the spin-orbital counterpart of
    !> Method_Mb_OrbBased_FillSpatialRdm1).
    module subroutine Method_Mb_OrbBased_FillOneBodyMatrix(h1, orbs, hOrbs)
      !> Output: one-body matrix (nOrbsSum x nOrbsSum).
      complex(R64), intent(out), allocatable :: h1(:, :)
      !> Input: stored orbitals (nPoints x nOrbsInState).
      complex(R64), intent(in), contiguous   :: orbs(:, :)
      !> Input: the operator applied to them (nPoints x nOrbsInState).
      complex(R64), intent(in), contiguous   :: hOrbs(:, :)
    end subroutine

    !> Fills the regularized inverse of the spatial 1-RDM,
    !> M = U n/(n^2 + eps) U^dagger (natural orbitals U, occupations n,
    !> eps = Method_Mb_OrbBased_regularizationParameter), spin-summed in the
    !> restricted case; the weights contracted with the 2-RDM in
    !> Method_Mb_OrbBased_ApplyCorrelationOp. The natural orbitals are
    !> computed within the blocks Orbs_blockOfOrb, so M is exactly block
    !> diagonal.
    module subroutine Method_Mb_OrbBased_FillInverseRdm1(invRdm1, rdm1)
      !> Output: inverse spatial 1-RDM (nOrbsInState x nOrbsInState).
      complex(R64), intent(out), allocatable :: invRdm1(:, :)
      !> One-body reduced density matrix (nOrbsSum x nOrbsSum).
      complex(R64), intent(in), contiguous   :: rdm1(:, :)
    end subroutine

    !> Writes the body-type blocks rdm2BtXX_YY.in of the 2-RDM (the format
    !> read by twoRdmInit.load).
    module subroutine Method_Mb_OrbBased_SaveRdm2Blocks(rdm2)
      !> Two-body reduced density matrix (nOrbsSum^4).
      complex(R64), intent(in), contiguous :: rdm2(:, :, :, :)
    end subroutine

    !> Eigen variant of Method_Mb_OrbBased_PropagateFieldFree (Ylm grids):
    !> orbs <- exp(-i h0 dt) orbs, or exp(-h0 dt) orbs if imagTimeQ, in the
    !> eigenbasis of the blocks of Method_Mb_OrbBased_SetupFieldFreeEigen.
    module subroutine Method_Mb_OrbBased_PropagateFieldFreeEigen(orbs, dt, imagTimeQ)
      !> Orbitals (nPoints x nOrbsInState), propagated in place.
      complex(R64), intent(inout), contiguous :: orbs(:, :)
      !> Time step (may be negative).
      real(R64), intent(in) :: dt
      !> Imaginary-time propagation exp(-h0 dt).
      logical, intent(in) :: imagTimeQ
    end subroutine

    !> Cayley variant of Method_Mb_OrbBased_PropagateFieldFree (grids with
    !> SysKinetic_SolveShifted): orbs <- (1 + i dt h0/2)^-1 (1 - i dt h0/2) orbs
    !> per orbital, solved by right-preconditioned GMRES in the grid metric
    !> (M_Utils_Gmres); h0 + U instead of h0 once the frozen mean-field
    !> reference U is built (exactFieldFreeMeanFieldQ). Real time only (error
    !> stop if imagTimeQ). Not reentrant (shared GMRES context, work vectors
    !> and step data; error stop inside an active OpenMP parallel region).
    module subroutine Method_Mb_OrbBased_PropagateFieldFreeCayley(orbs, dt, imagTimeQ)
      !> Orbitals (nPoints x nOrbsInState), propagated in place.
      complex(R64), intent(inout), contiguous :: orbs(:, :)
      !> Time step (may be negative).
      real(R64), intent(in) :: dt
      !> Imaginary-time propagation (not supported: error stop).
      logical, intent(in) :: imagTimeQ
    end subroutine

    !> Internal plumbing of the Cayley exact part (called only by
    !> S_Method_Mb_OrbBased_FieldFreeCayley, not for programs): the setup of
    !> the frozen mean-field reference, called at the end of
    !> Method_Mb_OrbBased_SetupFieldFreeCayley, so every setup of the Cayley
    !> part releases the reference. Releases a reference built before and,
    !> with exactFieldFreeMeanFieldQ, checks that the interaction is
    !> body-type independent and prints the mean-field line of the setup.
    !> The reference and its procedures live in
    !> S_Method_Mb_OrbBased_FieldFreeMeanField.
    module subroutine Method_Mb_OrbBased_SetupFieldFreeMeanField
    end subroutine

    !> Builds the frozen mean-field reference U = J0 - K0 of the Cayley exact
    !> part (method.mb.orbBased.exactFieldFreeMeanFieldQ; replaces one built
    !> earlier): J0 the Hartree potential of the density of rdm1 (one
    !> potential solve of the total source, SysInteraction_FillDirectPotential),
    !> K0 the exchange operator of rdm1 within every body type (restricted
    !> storage: the spin-up block for all stored orbitals), compressed by
    !> SysInteraction_FillExchangeAce on the stored orbitals of that body
    !> type (n(n+1)/2 pair potentials). For a single determinant h0 + U is
    !> its Fock operator. From then on the Cayley step propagates h0 + U
    !> (S_Method_Mb_OrbBased_FieldFreeCayley, through
    !> Method_Mb_OrbBased_ApplyFieldFreeMeanFieldOnOrb).
    !> MCTDHX builds it from the state at its first exact-part (or nonlinear)
    !> call after Method_Setup unless a program called this before (for a
    !> kick: from the ground state before the kick, which leaves no static
    !> part in the remainder; continued runs: from the same state again;
    !> A_Method.md, "Which state").
    module subroutine Method_Mb_OrbBased_BuildFieldFreeMeanField(orbs, rdm1, time)
      !> Stored orbitals of the state (nPoints x nOrbsInState).
      complex(R64), intent(in), contiguous :: orbs(:, :)
      !> One-body reduced density matrix of the state (nOrbsSum x nOrbsSum,
      !> rdm1(i, j) = <a+_j a_i>).
      complex(R64), intent(in), contiguous :: rdm1(:, :)
      !> Time of the state (interaction potentials).
      real(R64), intent(in) :: time
    end subroutine

    !> dOrbs(:, i) += U orbs(:, i) with the frozen mean-field reference U of
    !> the body type of orbital i (U = J0 - K0, K0 by its ACE vectors); with
    !> h1_ also its matrix h1_(i, j) = <orbs_i|U orbs_j> (equal body types;
    !> nOrbsSum x nOrbsSum, restricted: both spin blocks), the convention of
    !> Method_Mb_OrbBased_FillH1 (Method_Mb_OrbBased_FillOneBodyMatrix); the
    !> orbital form used by the MCTDHX nonlinear part. Error stop if the
    !> reference is not built.
    module subroutine Method_Mb_OrbBased_ApplyFieldFreeMeanField(dOrbs, orbs, h1_)
      !> In/out: orbitals to accumulate into (nPoints x nOrbsInState).
      complex(R64), intent(inout), contiguous :: dOrbs(:, :)
      !> Input: orbitals the reference acts on (nPoints x nOrbsInState).
      complex(R64), intent(in), contiguous :: orbs(:, :)
      !> Optional output: matrix of U in the orbitals.
      complex(R64), intent(out), allocatable, optional :: h1_(:, :)
    end subroutine

    !> Internal plumbing of the Cayley exact part (called only by
    !> S_Method_Mb_OrbBased_FieldFreeCayley, not for programs): dOrb += U orb
    !> with the frozen mean-field reference U of body type bt once it is built
    !> and exactFieldFreeMeanFieldQ is set, a no-op otherwise; the per-vector
    !> form used inside the Cayley step (its Arnoldi operator and true
    !> residual). SysInteraction_AddMeanField: J0 orb is formed in scratch,
    !> the ACE terms follow in one fused pass; reentrant.
    module subroutine Method_Mb_OrbBased_ApplyFieldFreeMeanFieldOnOrb(dOrb, scratch, orb, bt)
      !> In/out: vector to accumulate into (nPoints).
      complex(R64), intent(inout), contiguous :: dOrb(:)
      !> Work vector (nPoints), overwritten when U is applied.
      complex(R64), intent(inout), contiguous :: scratch(:)
      !> Input: vector the reference acts on (nPoints).
      complex(R64), intent(in), contiguous :: orb(:)
      !> Body type whose ACE vectors act on orb.
      integer(I32), intent(in) :: bt
    end subroutine

    !> Whether the frozen mean-field reference is built (since the last
    !> Method_Mb_OrbBased_SetupFieldFreeMeanField, which releases it).
    module function Method_Mb_OrbBased_FieldFreeMeanFieldBuiltQ() result(res)
      logical :: res
    end function

    !> Eigen variant of Method_Mb_OrbBased_ApplyFieldFreePhi (Ylm grids):
    !> orbs = sum_m phi_m(-i h0 dt) orbsList(:, :, m), evaluated in the
    !> h0 eigenbasis (kept eigenvectors only) (S_Method_Mb_OrbBased_FieldFreeEigen).
    module subroutine Method_Mb_OrbBased_ApplyFieldFreePhiEigen(orbs, orbsList, dt)
      !> Output: orbitals (nPoints x nOrbsInState).
      complex(R64), intent(out), contiguous :: orbs(:, :)
      !> Input: orbitals (nPoints x nOrbsInState x 0:K), slice m multiplies phi_m.
      complex(R64), intent(in), contiguous :: orbsList(:, :, 0:)
      !> Time step.
      real(R64), intent(in) :: dt
    end subroutine

    !> Eigen variant of Method_Mb_OrbBased_ProjectBelowFieldFreeCut (Ylm
    !> grids): projects orbitals (or orbital derivatives) onto the h0
    !> eigenvectors below `Method_Mb_OrbBased_exactFieldFreeCut` (spectral
    !> cut), per (l,m) block (S_Method_Mb_OrbBased_FieldFreeEigen).
    module subroutine Method_Mb_OrbBased_ProjectBelowFieldFreeCutEigen(orbs)
      !> Orbitals (nPoints x nOrbsInState), projected in place.
      complex(R64), intent(inout), contiguous :: orbs(:, :)
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Classification of basis states by body type to enable type-dependent
  !> potentials and interactions (e.g., spin-dependent terms).
  !> Method_Mb_OrbBased_bodyTypeOfOrb(i) gives the body type of basis state i.
  integer(I32), allocatable :: Method_Mb_OrbBased_bodyTypeOfOrb(:)

  !> Declared m of every spin orbital (nOrbsSum, layout of nOrbs), from
  !> method.mb.orbBased.mOfOrbs; allocated only if declared. Every orbital
  !> lives in the Ylm channels of its m (A_Method.md, "Declared m").
  integer(I32), allocatable :: Method_Mb_OrbBased_mOfOrb(:)

  !> True if method.mb.orbBased.mOfOrbs is declared.
  logical :: Method_Mb_OrbBased_definiteMQ = .false.

  !> Number of one-body orbitals per body type.
  !> Method_Mb_OrbBased_nOrbs(i) gives the number of basis states for body type i.
  integer(I32), allocatable :: Method_Mb_OrbBased_nOrbs(:)

  !> Start index of one-body orbitals per body type.
  !> Method_Mb_OrbBased_nOrbsStart(i) gives the start index for basis states of type i.
  integer(I32), allocatable :: Method_Mb_OrbBased_nOrbsStart(:)

  !> End index of one-body orbitals per body type.
  !> Method_Mb_OrbBased_nOrbsEnd(i) gives the end index for basis states of type i.
  integer(I32), allocatable :: Method_Mb_OrbBased_nOrbsEnd(:)

  !> Total number of one-body basis states across all body types.
  !> May differ from Orbs_nOrbsInState in restricted calculations where
  !> the same orbital is used for different body types.
  integer(I32) :: Method_Mb_OrbBased_nOrbsSum

  !> One-body reduced density matrix in the orbital basis.
  !> Method_Mb_OrbBased_rdm1(i,j) = \(\langle\Psi| a^\dagger_j a_i |\Psi\rangle\)
  !> where i is the lower index and j is the upper index. Hermitian by
  !> construction: rdm1(i,j) = conjg(rdm1(j,i)).
  complex(R64), allocatable :: Method_Mb_OrbBased_rdm1(:, :)

  !> Two-body reduced density matrix in the orbital basis.
  !> Method_Mb_OrbBased_rdm2(i1,i2,j1,j2) = \(\langle\Psi|
  !> a^\dagger_{j1} a^\dagger_{j2} a_{i2} a_{i1} |\Psi\rangle\) with
  !> (i1,i2) lower and (j1,j2) upper indices. Symmetry properties follow
  !> particle statistics (antisymmetric for fermions, symmetric for bosons).
  complex(R64), allocatable :: Method_Mb_OrbBased_rdm2(:, :, :, :)

  !> One-body Hamiltonian matrix in the orbital basis.
  !> Method_Mb_OrbBased_h1(i,j) = \(\langle i| H_1 |j\rangle\) representing
  !> kinetic energy and external potentials.
  complex(R64), allocatable :: Method_Mb_OrbBased_h1(:, :)

  !> Regularization parameter for numerical stability.
  !> Used to avoid small denominators or ill-conditioned inversions.
  real(R64) :: Method_Mb_OrbBased_regularizationParameter

  !> Two-body Hamiltonian in the orbital basis.
  !> Method_Mb_OrbBased_h2(i1,i2,j1,j2) = \(\langle i1,i2| H_2 |j1,j2\rangle\)
  !> representing interaction energy between particles.
  complex(R64), allocatable :: Method_Mb_OrbBased_h2(:, :, :, :)

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> Pointer to the procedure for calculating the one-body Hamiltonian matrix.
  procedure(I_Method_Mb_OrbBased_FillH1), pointer :: Method_Mb_OrbBased_FillH1
  abstract interface
    !> Allocates and fills the one-body Hamiltonian matrix.
    !> The matrix elements include kinetic and external potential terms.
    !> \[
    !> H^1_{ij} = \langle \phi_i | \hat{T} + \hat{V} | \phi_j \rangle
    !> \]
    subroutine I_Method_Mb_OrbBased_FillH1(h1, orbs, time)
      import :: R64
      !> Output: one-body Hamiltonian matrix, allocated within the subroutine.
      complex(R64), intent(out), allocatable :: h1(:, :)
      !> Input: orbitals defining the basis for the Hamiltonian (column-major).
      complex(R64), intent(in), contiguous :: orbs(:, :)
      !> Current time, allowing for time-dependent potentials.
      real(R64), intent(in) :: time
    end subroutine
  end interface

  !> Pointer to the procedure for calculating the two-body Hamiltonian matrix.
  procedure(I_Method_Mb_OrbBased_FillH2), pointer :: Method_Mb_OrbBased_FillH2
  abstract interface
    !> Allocates and fills the two-body Hamiltonian matrix from orbital integrals.
    !> \[
    !> H^2_{i_1 i_2 j_1 j_2} = \int \phi_{i_1}^*(r_1) \phi_{i_2}^*(r_2) V(r_1,r_2) \phi_{j_1}(r_1) \phi_{j_2}(r_2) dr_1 dr_2
    !> \]
    subroutine I_Method_Mb_OrbBased_FillH2(h2, orbs, time)
      import :: R64
      !> Output: two-body Hamiltonian matrix, allocated within the subroutine.
      complex(R64), intent(out), allocatable :: h2(:, :, :, :)
      !> Input: orbitals defining the basis for the Hamiltonian (column-major).
      complex(R64), intent(in), contiguous :: orbs(:, :)
      !> Current time, allowing for time-dependent interactions.
      real(R64), intent(in) :: time
    end subroutine
  end interface

  !> Pointer to the procedure for calculating the one-body reduced density matrix.
  procedure(I_Method_Mb_OrbBased_FillRdm1), pointer :: Method_Mb_OrbBased_FillRdm1
  abstract interface
    subroutine I_Method_Mb_OrbBased_FillRdm1(rdm1, state)
      import :: R64
      !> Output: one-body reduced density matrix, allocated within the subroutine.
      complex(R64), intent(out), allocatable :: rdm1(:, :)
      !> Input: packed quantum state from which rdm1 is computed.
      complex(R64), intent(in), contiguous, target  :: state(:)
    end subroutine
  end interface

  !> Pointer to the procedure for calculating the two-body reduced density matrix.
  procedure(I_Method_Mb_OrbBased_FillRdm2), pointer :: Method_Mb_OrbBased_FillRdm2
  abstract interface
    subroutine I_Method_Mb_OrbBased_FillRdm2(rdm2, state)
      import :: R64
      !> Output: two-body reduced density matrix, allocated within the subroutine.
      complex(R64), intent(out), allocatable :: rdm2(:, :, :, :)
      !> Input: packed quantum state from which rdm2 is computed.
      complex(R64), intent(in), contiguous, target  :: state(:)
    end subroutine
  end interface

  !> Pointer to the procedure for applying the kinetic operator to orbitals.
  procedure(I_Method_Mb_OrbBased_ApplyKineticOp), pointer :: Method_Mb_OrbBased_ApplyKineticOp
  abstract interface
    !> Applies the kinetic operator to the orbitals (e.g., second-derivative or
    !> tight-binding hopping contribution) and accumulates in `dOrbs`.
    subroutine I_Method_Mb_OrbBased_ApplyKineticOp(dOrbs, orbs, time)
      import :: R64
      !> In/out: derivative orbitals after applying the kinetic operator (accumulate).
      complex(R64), intent(inout), contiguous :: dOrbs(:, :)
      !> Input: orbitals to which the kinetic operator is applied.
      complex(R64), intent(in), contiguous :: orbs(:, :)
      !> Current time, allowing for time-dependence.
      real(R64), intent(in) :: time
    end subroutine
  end interface

  !> Pointer to the procedure for applying the potential energy operator.
  procedure(I_Method_Mb_OrbBased_ApplyPotentialOp), pointer :: Method_Mb_OrbBased_ApplyPotentialOp
  abstract interface
    !> Applies the potential energy operator to a set of orbitals and accumulates
    !> the result in `dOrbs` (dOrbs += V(orbs)).
    subroutine I_Method_Mb_OrbBased_ApplyPotentialOp(dOrbs, orbs, time)
      import :: I32, R64
      !> In/out: where the result is added. Initialize before calling.
      complex(R64), intent(inout), contiguous :: dOrbs(:, :)
      !> Input: orbitals to which the potential operator is applied.
      complex(R64), intent(in), contiguous :: orbs(:, :)
      !> Current time, allowing for time-dependence.
      real(R64), intent(in) :: time
    end subroutine
  end interface

  !> Pointer to the procedure for applying the velocity-gauge coupling operator.
  procedure(I_Method_Mb_OrbBased_ApplyGaugeOp), pointer :: Method_Mb_OrbBased_ApplyGaugeOp
  abstract interface
    !> Applies the gauge coupling operator (A(t)·p/m + A²/(2m), velocity gauge)
    !> to a set of orbitals and accumulates the result in `dOrbs`.
    subroutine I_Method_Mb_OrbBased_ApplyGaugeOp(dOrbs, orbs, time)
      import :: R64
      !> In/out: where the result is added. Initialize before calling.
      complex(R64), intent(inout), contiguous :: dOrbs(:, :)
      !> Input: orbitals to which the gauge operator is applied.
      complex(R64), intent(in), contiguous :: orbs(:, :)
      !> Current time at which the vector potential is evaluated.
      real(R64), intent(in) :: time
    end subroutine
  end interface

  !> Pointer to the procedure for applying the two-body interaction operator.
  procedure(I_Method_Mb_OrbBased_ApplyInteractionOp), pointer :: Method_Mb_OrbBased_ApplyInteractionOp
  abstract interface
    !> Applies the two-body interaction operator to the orbitals (e.g., Coulomb
    !> or Hubbard-like terms) and accumulates the result in `dOrbs`.
    subroutine I_Method_Mb_OrbBased_ApplyInteractionOp(dOrbs, orbs, i2, j2, time)
      import :: I32, R64
      complex(R64), intent(inout), contiguous :: dOrbs(:, :)
      complex(R64), intent(in), contiguous   :: orbs(:, :)
      integer(I32), intent(in)               :: i2
      integer(I32), intent(in)               :: j2
      !> Current time, allowing for time-dependence.
      real(R64), intent(in) :: time
    end subroutine
  end interface

  !> Pointer to the procedure for applying a correlation operator (beyond mean-field).
  procedure(I_Method_Mb_OrbBased_ApplyCorrelationOp), pointer :: Method_Mb_OrbBased_ApplyCorrelationOp
  abstract interface
    !> Applies a correlation operator that depends on the reduced density matrices
    !> (typical of multiconfigurational methods) and accumulates the result in `dOrbs`.
    subroutine I_Method_Mb_OrbBased_ApplyCorrelationOp(dOrbs, orbs, invRdm1, rdm2, time, h2_)
      import :: R64
      !> Input/output derivative orbitals. On exit, contains the original
      !> content plus the result of the correlation operator action.
      complex(R64), intent(inout), contiguous, target :: dOrbs(:, :)
      !> Input orbitals to which the operator is applied.
      complex(R64), intent(in), contiguous, target :: orbs(:, :)
      !> Inverse spatial 1-RDM (nOrbsInState x nOrbsInState) contracted with
      !> the 2-RDM, e.g. from Method_Mb_OrbBased_FillInverseRdm1.
      complex(R64), intent(in), contiguous :: invRdm1(:, :)
      !> Two-body reduced density matrix.
      complex(R64), intent(in), contiguous :: rdm2(:, :, :, :)
      !> Current time, allowing for time-dependent potentials.
      real(R64), intent(in) :: time
      !> Optional output: two-body Hamiltonian matrix elements (same layout as h2).
      complex(R64), intent(out), allocatable, optional :: h2_(:, :, :, :)
    end subroutine
  end interface

  !> Pointer to the procedure for applying the mean-field (Hartree-Fock) operator.
  procedure(I_Method_Mb_OrbBased_ApplyHartreeFockOp), pointer :: Method_Mb_OrbBased_ApplyHartreeFockOp
  abstract interface
    !> Applies the Hartree-Fock (mean-field + exchange) operator to a set of orbitals
    !> and accumulates the result in `dOrbs`.
    subroutine I_Method_Mb_OrbBased_ApplyHartreeFockOp(dOrbs, orbs, rdm1, time)
      import :: R64
      !> Input/output derivative orbitals. On exit, contains the original
      !> content plus the result of the mean-field operator action.
      complex(R64), intent(inout), contiguous, target :: dOrbs(:, :)
      !> Input orbitals to which the operator is applied.
      complex(R64), intent(in), contiguous, target :: orbs(:, :)
      !> One-body reduced density matrix.
      complex(R64), intent(in), contiguous :: rdm1(:, :)
      !> Current time, allowing for time-dependent potentials.
      real(R64), intent(in) :: time
    end subroutine
  end interface

  !> Pointer to the procedure for applying the one-body (single-particle) operator.
  procedure(I_Method_Mb_OrbBased_ApplySingleBodyOp), pointer :: Method_Mb_OrbBased_ApplySingleBodyOp
  abstract interface
    !> Applies the one-body part of the Hamiltonian (kinetic + external potential)
    !> to the orbitals and accumulates the result in `dOrbs`.
    subroutine I_Method_Mb_OrbBased_ApplySingleBodyOp(dOrbs, orbs, time, h1_)
      import :: R64
      !> Input/output derivative orbitals. On exit, contains the original
      !> content plus the result of the one-body operator action.
      complex(R64), intent(inout), contiguous, target :: dOrbs(:, :)
      !> Input orbitals to which the operator is applied.
      complex(R64), intent(in), contiguous, target :: orbs(:, :)
      !> Current time, allowing for time-dependence.
      real(R64), intent(in) :: time
      !> Optional output: one-body Hamiltonian matrix elements.
      complex(R64), intent(out), allocatable, optional :: h1_(:, :)
    end subroutine
  end interface

  !> If true (method.mb.orbBased.exactFieldFreeQ), the field-free one-body
  !> Hamiltonian h0 = T + V_ext is propagated exactly by
  !> `Method_PropagateExactPart` and the linear split-step integrator only
  !> propagates the gauge coupling (`Method_Mb_OrbBased_ApplyLinOp`).
  !> Implementation (bound in Method_Mb_OrbBased_Fabricate): on Ylm grids the
  !> dense eigen-decomposition of the h0 blocks (*Eigen procedures), on other
  !> grids whose kinetic backend binds SysKinetic_SolveShifted (grid.fem) the
  !> Crank-Nicolson (Cayley) step (*Cayley procedures; no phi functions, no
  !> spectral cut, no imaginary time).
  logical :: Method_Mb_OrbBased_exactFieldFreeQ = .false.

  !> If true (method.mb.orbBased.exactFieldFreeMeanFieldQ; Cayley exact part
  !> and MCTDHX only), the exact part propagates h0 + U with a frozen
  !> mean-field reference U = J0 - K0 (static, Hermitian; built once from the
  !> state, Method_Mb_OrbBased_BuildFieldFreeMeanField), and the nonlinear
  !> split-step part propagates the interaction minus U (orbitals
  !> -i Q (W - U) phi, CI coefficients (H2 - U^) C). The equations of motion
  !> are unchanged; the split is that of the Fock operator of the reference
  !> and a remainder that is small near it (A_Method.md, "Frozen mean-field
  !> reference").
  logical :: Method_Mb_OrbBased_exactFieldFreeMeanFieldQ = .false.

  !> Pointer to the setup of the exact field-free part (called by
  !> Method_Mb_OrbBased_Setup if the method binds Method_PropagateExactPart):
  !> Method_Mb_OrbBased_SetupFieldFreeEigen (Ylm) or
  !> Method_Mb_OrbBased_SetupFieldFreeCayley; null without exactFieldFreeQ.
  procedure(I_Method_Mb_OrbBased_SetupFieldFree), pointer :: Method_Mb_OrbBased_SetupFieldFree => null()
  abstract interface
    !> Prepares the exact field-free part (operators, caches, checks).
    subroutine I_Method_Mb_OrbBased_SetupFieldFree
    end subroutine
  end interface

  !> Pointer to the exact field-free propagation of the orbitals,
  !> orbs <- U0(dt) orbs: exp(-i h0 dt) by the h0 eigen-decomposition (Ylm,
  !> Method_Mb_OrbBased_PropagateFieldFreeEigen, also exp(-h0 dt) in imaginary
  !> time) or the Cayley approximation (1 + i dt h0/2)^-1 (1 - i dt h0/2)
  !> (Method_Mb_OrbBased_PropagateFieldFreeCayley, real time only); null
  !> without exactFieldFreeQ.
  procedure(I_Method_Mb_OrbBased_PropagateFieldFree), pointer :: Method_Mb_OrbBased_PropagateFieldFree => null()
  abstract interface
    !> Propagates the orbitals with the exact field-free part over dt.
    subroutine I_Method_Mb_OrbBased_PropagateFieldFree(orbs, dt, imagTimeQ)
      import :: R64
      !> Orbitals (nPoints x nOrbsInState), propagated in place.
      complex(R64), intent(inout), contiguous :: orbs(:, :)
      !> Time step (may be negative).
      real(R64), intent(in) :: dt
      !> Imaginary-time propagation exp(-h0 dt).
      logical, intent(in) :: imagTimeQ
    end subroutine
  end interface

  !> Pointer to the phi functions of the exact field-free part,
  !> orbs = sum_m phi_m(-i h0 dt) orbsList(:, :, m)
  !> (Method_Mb_OrbBased_ApplyFieldFreePhiEigen, Ylm grids); null for the
  !> Cayley part, which has none (MCTDHX then does not bind
  !> Method_ApplyExactPartPhi, so propagator.exponential is unavailable).
  procedure(I_Method_Mb_OrbBased_ApplyFieldFreePhi), pointer :: Method_Mb_OrbBased_ApplyFieldFreePhi => null()
  abstract interface
    !> Phi functions of h0 applied to a list of orbital sets.
    subroutine I_Method_Mb_OrbBased_ApplyFieldFreePhi(orbs, orbsList, dt)
      import :: R64
      !> Output: orbitals (nPoints x nOrbsInState).
      complex(R64), intent(out), contiguous :: orbs(:, :)
      !> Input: orbitals (nPoints x nOrbsInState x 0:K), slice m multiplies phi_m.
      complex(R64), intent(in), contiguous :: orbsList(:, :, 0:)
      !> Time step.
      real(R64), intent(in) :: dt
    end subroutine
  end interface

  !> Diagnostics of the Cayley exact part (non-Ylm grids): number of GMRES
  !> solves (one per orbital and exact substep), their total and largest
  !> iteration counts (preconditioned matrix-vector products, all restarts)
  !> since Method_Mb_OrbBased_SetupFieldFreeCayley.
  integer(I64) :: Method_Mb_OrbBased_cayleyNSolves = 0
  integer(I64) :: Method_Mb_OrbBased_cayleyNIterations = 0
  integer(I32) :: Method_Mb_OrbBased_cayleyMaxIterations = 0

  !> Spectral cut (method.mb.orbBased.exactFieldFreeCut, 0 = off): the
  !> orbitals are restricted to the h0 eigenvectors with energy below this
  !> value. The exact part drops the eigenvectors above it, and the linear and
  !> nonlinear orbital derivatives are projected with
  !> `Method_Mb_OrbBased_ApplySpectralCut`, so the unphysical high-energy
  !> FEDVR states are never populated.
  real(R64) :: Method_Mb_OrbBased_exactFieldFreeCut = 0.0_R64

  !> Pointer to the projection of orbital derivatives onto the h0 eigenvectors
  !> below the spectral cut (`Method_Mb_OrbBased_ProjectBelowFieldFreeCut`),
  !> or a no-op without cut.
  procedure(I_Method_Mb_OrbBased_ApplySpectralCut), pointer :: Method_Mb_OrbBased_ApplySpectralCut
  abstract interface
    !> Projects the orbitals in place onto the space kept by the spectral cut.
    subroutine I_Method_Mb_OrbBased_ApplySpectralCut(orbs)
      import :: R64
      !> In/out: orbitals (nPoints x nOrbsInState), projected in place.
      complex(R64), intent(inout), contiguous :: orbs(:, :)
    end subroutine
  end interface

  !> Pointer to the projection onto the h0 eigenvectors below the spectral
  !> cut (Method_Mb_OrbBased_ProjectBelowFieldFreeCutEigen, Ylm grids); null
  !> for the Cayley part (exactFieldFreeCut stops at fabricate on non-Ylm
  !> grids).
  procedure(I_Method_Mb_OrbBased_ApplySpectralCut), pointer :: Method_Mb_OrbBased_ProjectBelowFieldFreeCut => null()

  !> Pointer to the one-body operator propagated by the linear split-step
  !> integrator (`TimeDerivativeOrbsLin`): the full one-body Hamiltonian, or
  !> only the gauge coupling if `Method_Mb_OrbBased_exactFieldFreeQ`.
  procedure(I_Method_Mb_OrbBased_ApplyLinOp), pointer :: Method_Mb_OrbBased_ApplyLinOp
  abstract interface
    !> Applies the one-body operator propagated by the linear split-step
    !> integrator (full one-body Hamiltonian, or only the gauge coupling if
    !> `Method_Mb_OrbBased_exactFieldFreeQ`) and accumulates the result in `dOrbs`.
    subroutine I_Method_Mb_OrbBased_ApplyLinOp(dOrbs, orbs, time)
      import :: R64
      !> In/out: where the result is added. Initialize before calling.
      complex(R64), intent(inout), contiguous :: dOrbs(:, :)
      !> Input: orbitals to which the operator is applied.
      complex(R64), intent(in), contiguous :: orbs(:, :)
      !> Current time.
      real(R64), intent(in) :: time
    end subroutine
  end interface

  !> Pointer to the linear part of a split-step time derivative: the one-body
  !> operator acting on the orbitals only.
  procedure(I_Method_Mb_OrbBased_TimeDerivativeOrbsLin), pointer :: &
    Method_Mb_OrbBased_TimeDerivativeOrbsLin
  abstract interface
    !> Linear split-step part: dOrbs = -i h orbs (unprojected), zero for the
    !> many-body part (CI coefficients or 2-RDM). Evolving every orbital with
    !> exp(-i h t) at a fixed many-body part is the exact one-body evolution.
    subroutine I_Method_Mb_OrbBased_TimeDerivativeOrbsLin(dState, state, time)
      import :: R64
      !> Output: time derivative (multiplied by -i).
      complex(R64), intent(out), contiguous, target :: dState(:)
      !> Input: quantum state.
      complex(R64), intent(in), contiguous, target  :: state(:)
      !> Current time (allows for time-dependent Hamiltonians).
      real(R64), intent(in) :: time
    end subroutine
  end interface

  !> Pointer to the non-linear part of a split-step time derivative when both the
  !> many-body part and the orbitals evolve (MCTDHX, TD-2RDM): only the
  !> two-body interaction, i.e. the many-body part evolves with h2 (plus the
  !> collision term for TD-2RDM) and the orbitals with Q (W D^-1 terms).
  !> Together with `TimeDerivativeOrbsLin` it is a gauge transform of the full
  !> time derivative (in-space orbital rotation by h1, compensated in the
  !> many-body part).
  procedure(I_Method_Mb_OrbBased_TimeDerivativeOrbsLin), pointer :: &
    Method_Mb_OrbBased_TimeDerivativeCoeffsPlusOrbsNonLin

  !> Pointer to the Loewdin orthonormalization of the orbitals with the
  !> compensating transformation of the many-body part, which leaves the
  !> many-body state unchanged (MCTDHX; error stop for the other methods).
  !> Used after absorber masks (absorber.keepOrthonormalQ) and in imaginary time.
  procedure(I_Method_Mb_OrbBased_OrthonormalizeKeepingState), pointer :: &
    Method_Mb_OrbBased_OrthonormalizeKeepingState
  abstract interface
    !> phi <- phi S^(-1/2) per symmetry block (Orbs_blockOfOrb) and the
    !> many-body part transformed so that the state is unchanged; the norm
    !> lost by the orbitals moves into the many-body part.
    subroutine I_Method_Mb_OrbBased_OrthonormalizeKeepingState(state)
      import :: R64
      !> In/out: quantum state.
      complex(R64), intent(inout), contiguous, target :: state(:)
    end subroutine
  end interface

end module
