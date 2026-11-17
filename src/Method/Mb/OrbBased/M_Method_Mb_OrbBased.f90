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
!> Submodules: S_Method_Mb_OrbBased (fabricate, operators, h1 and h2) with
!> its child S_Method_Mb_OrbBased_YlmNodes (the definite-m node path of the
!> correlation operator and of h2 on Ylm grids; a nested submodule adds
!> nothing to this interface, what its parent calls is declared in the
!> parent); S_Method_Mb_OrbBased_MeanField (the frozen mean field of a 1-RDM)
!> and S_Method_Mb_OrbBased_Switching, which switches the interaction on
!> adiabatically (method.mb.orbBased.adiabaticSwitching) with such a mean
!> field as its reference.
!>
!> The exact field-free part of the orbital propagation is its own module,
!> M_FieldFree ("fieldFree" block): with it the methods bind
!> Method_PropagateExactPart to FieldFree_Propagate and
!> Method_Mb_OrbBased_ApplyLinOp carries only the gauge coupling.
module M_Method_Mb_OrbBased
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! types
  !=============================================================================

  !> Frozen Hartree-Fock mean field U = J0 - K0 of a 1-RDM on the grid
  !> (Method_Mb_OrbBased_FillMeanField): the Hartree potential J0 (in the
  !> representation of the interaction backend) and the ACE vectors of K0,
  !> K0 v ~ sum_k zeta_k <zeta_k|v>, exact on the span of the orbitals it was
  !> built from; the columns zetaStart(bt) .. zetaEnd(bt) act on the orbitals
  !> of body type bt.
  type :: T_Method_Mb_OrbBased_MeanField
    complex(R64), allocatable :: hartree(:)
    complex(R64), allocatable :: zeta(:, :)
    integer(I32), allocatable :: zetaStart(:)
    integer(I32), allocatable :: zetaEnd(:)
  end type

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
    !> the skeleton of the definite-m node path if it is bound.
    module subroutine Method_Mb_OrbBased_Setup
    end subroutine

    !> Whether the one-body Hamiltonian of Method_Mb_OrbBased_ApplySingleBodyOp
    !> (and of Method_Mb_OrbBased_FillH1) is time-independent: the kinetic
    !> operator, the external potential and, if bound (a sysGauge block
    !> exists), the gauge coupling all are (their *_timeIndependentQ flags,
    !> read at call time), and there is no adiabatic switching. The methods
    !> with static orbitals (TDCI, TD-2RDM-StaticOrbs) refill h1 at every time
    !> derivative otherwise (Method_Mb_OrbBased_RefillStaticH1).
    module function Method_Mb_OrbBased_H1TimeIndependentQ() result(res)
      logical :: res
    end function

    !> The one-body matrix of T + V_ext without the gauge coupling, in the
    !> layout of Method_Mb_OrbBased_FillH1.
    module subroutine Method_Mb_OrbBased_FillH0(h0, orbs, time)
      complex(R64), intent(out), allocatable :: h0(:, :)
      complex(R64), intent(in), contiguous   :: orbs(:, :)
      real(R64), intent(in)                  :: time
    end subroutine

    !> h1 at time for static orbitals (TDCI, TD-2RDM-StaticOrbs): unchanged if
    !> h1 is time-independent; with a static kinetic operator and external
    !> potential the static part (h1 at the first call after
    !> Method_Mb_OrbBased_Setup: the Method_Mb_OrbBased_FillH0 matrix of the
    !> method setup, or what a program stored instead) plus the gauge coupling;
    !> otherwise Method_Mb_OrbBased_FillH1; then the switching reference
    !> (Method_Mb_OrbBased_AddSwitchingReferenceMatrix). The orbitals must not
    !> change after the first call.
    module subroutine Method_Mb_OrbBased_RefillStaticH1(h1, orbs, time)
      complex(R64), intent(inout), allocatable :: h1(:, :)
      complex(R64), intent(in), contiguous     :: orbs(:, :)
      real(R64), intent(in)                    :: time
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

    !> Diagonal metric of the packed state Method_state (after Method_Setup),
    !> the metricWeights of integrators and diagonalizers that act on it: 1
    !> for the many-body part (CI coefficients or 2-RDM), Grid_weights for
    !> every orbital of the state (moving orbitals, the tail of Method_state).
    module subroutine Method_Mb_OrbBased_FillStateWeights(weights)
      !> Output: weights (size of Method_state).
      real(R64), intent(out), allocatable :: weights(:)
    end subroutine

    !> CI Hamiltonian of the stored matrices Method_Mb_OrbBased_h1 and
    !> Method_Mb_OrbBased_h2 on a coefficient vector, dCoeffs = H coeffs (no
    !> refill: callers fill h1/h2 of their orbitals first). The matvec of CI
    !> diagonalizers (DiagonalizerList ApplyMatOnVec signature; time unused)
    !> and the CI part of TDCI.
    module subroutine Method_Mb_OrbBased_ApplyCiHamiltonian(dCoeffs, coeffs, time)
      !> Output: H coeffs (Coeffs_nCoeffs).
      complex(R64), intent(out), contiguous, target :: dCoeffs(:)
      !> Input: CI coefficient vector (Coeffs_nCoeffs).
      complex(R64), intent(in), contiguous, target  :: coeffs(:)
      !> Unused (signature of the matvec callbacks).
      real(R64), intent(in) :: time
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

    !> The frozen Hartree-Fock mean field U = J0 - K0 of rdm1 on the stored
    !> orbitals (fermions, body-type-independent interaction): J0 the
    !> Hartree potential of the density of rdm1 (one potential solve of the
    !> total source, SysInteraction_FillDirectPotential), K0 the exchange
    !> operator of rdm1 within every body type (restricted storage: the
    !> spin-up block for all stored orbitals), compressed by
    !> SysInteraction_FillExchangeAce (n(n+1)/2 pair potentials). For a single
    !> determinant h0 + U is its Fock operator
    !> (S_Method_Mb_OrbBased_MeanField).
    module subroutine Method_Mb_OrbBased_FillMeanField(mf, orbs, rdm1, time, nPairSolves_, exchangeOnlyQ_)
      !> Output: the mean field.
      type(T_Method_Mb_OrbBased_MeanField), intent(out) :: mf
      !> Stored orbitals (nPoints x nOrbsInState).
      complex(R64), intent(in), contiguous :: orbs(:, :)
      !> One-body reduced density matrix (nOrbsSum x nOrbsSum,
      !> rdm1(i, j) = <a+_j a_i>).
      complex(R64), intent(in), contiguous :: rdm1(:, :)
      !> Time of the interaction potentials.
      real(R64), intent(in) :: time
      !> Optional output: number of pair potentials solved for K0.
      integer(I32), intent(out), optional :: nPairSolves_
      !> Only the ACE vectors of K0, no Hartree potential (default: both).
      logical, intent(in), optional :: exchangeOnlyQ_
    end subroutine

    !> w += U z with the mean field mf of body type bt
    !> (SysInteraction_AddMeanField; scratch: work vector of nPoints).
    !> Reentrant; its passes are parallel when called outside a parallel
    !> region.
    module subroutine Method_Mb_OrbBased_AddMeanField(w, scratch, mf, z, bt)
      complex(R64), intent(inout), contiguous :: w(:)
      complex(R64), intent(inout), contiguous :: scratch(:)
      type(T_Method_Mb_OrbBased_MeanField), intent(in) :: mf
      complex(R64), intent(in), contiguous :: z(:)
      integer(I32), intent(in) :: bt
    end subroutine

    !> Switching factor lambda(t) of the interaction
    !> (method.mb.orbBased.adiabaticSwitching, S_Method_Mb_OrbBased_Switching):
    !> sin^2(pi (t + T) / (2 T)) for -T < t < 0 (T = switchTime), 0 before,
    !> 1 from t = 0 on and without switching. The Hamiltonian of the method is
    !> h + (1 - lambda) U0 + lambda W.
    module function Method_Mb_OrbBased_SwitchFactor(time) result(res)
      real(R64), intent(in) :: time
      real(R64) :: res
    end function

    !> Builds the reference U0 of the switching from the current state
    !> (Method_state, Orbs_orbs): the frozen mean field of its 1-RDM, on the
    !> grid for moving orbitals (Method_Mb_OrbBased_FillMeanField), as a matrix
    !> of Method_Mb_OrbBased_h2 for static orbitals. A no-op without
    !> switching. Called at the end of the setup of every method that
    !> supports switching; programs that change the state (orthonormalization)
    !> or the static Hamiltonian (Method_Mb_OrbBased_h2) after Method_Setup
    !> call it again before the propagation.
    module subroutine Method_Mb_OrbBased_SetSwitchingReference
    end subroutine

    !> dOrbs(:, i) += (1 - lambda(time)) U0 orbs(:, i) with the grid reference
    !> (moving orbitals); a no-op if it is not built or lambda = 1. Part of
    !> the one-body operator (Method_Mb_OrbBased_ApplySingleBodyOp and the
    !> linear part with an exact field-free part).
    module subroutine Method_Mb_OrbBased_AddSwitchingReference(dOrbs, orbs, time)
      complex(R64), intent(inout), contiguous :: dOrbs(:, :)
      complex(R64), intent(in), contiguous :: orbs(:, :)
      real(R64), intent(in) :: time
    end subroutine

    !> h1 += (1 - lambda(time)) U0 with the matrix reference (static
    !> orbitals, Method_Mb_OrbBased_RefillStaticH1); a no-op if it is not
    !> built or lambda = 1.
    module subroutine Method_Mb_OrbBased_AddSwitchingReferenceMatrix(h1, time)
      complex(R64), intent(inout), contiguous :: h1(:, :)
      real(R64), intent(in) :: time
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

  !> Switching time T of method.mb.orbBased.adiabaticSwitching (0: no
  !> switching; Method_Mb_OrbBased_SwitchFactor).
  real(R64) :: Method_Mb_OrbBased_switchTime = 0.0_R64

  !> Whether the method keeps the orbitals static (TDCI, TD-2RDM-StaticOrbs;
  !> set by their fabricate): h1 and h2 are matrices of fixed orbitals.
  logical :: Method_Mb_OrbBased_staticOrbsQ = .false.

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

  !> Options that only some methods implement. The fabricate of a method sets
  !> the flags of those it implements; Method_Mb_OrbBased_Fabricate resets
  !> them before it dispatches and refuses a requested option afterwards if
  !> its flag is unset:
  !>   - the frozen mean-field reference of the exact part
  !>     (fieldFree.cayley.meanFieldQ): the method subtracts U in its
  !>     nonlinear part (FieldFree_ApplyMeanField);
  !>   - the spectral cut (fieldFree.eigen.cut): the method projects its
  !>     orbital derivatives with FieldFree_ProjectBelowCut;
  !>   - the adiabatic switching (method.mb.orbBased.adiabaticSwitching): the
  !>     method builds the reference at the end of its setup
  !>     (Method_Mb_OrbBased_SetSwitchingReference).
  !> The exact field-free part itself (a bound FieldFree_Propagate) is
  !> declared by binding Method_PropagateExactPart.
  logical :: Method_Mb_OrbBased_implementsMeanFieldReferenceQ = .false.
  logical :: Method_Mb_OrbBased_implementsSpectralCutQ = .false.
  logical :: Method_Mb_OrbBased_implementsSwitchingQ = .false.

  !> Pointer to the one-body operator propagated by the linear split-step
  !> integrator (`TimeDerivativeOrbsLin`): the full one-body Hamiltonian, or
  !> only the gauge coupling with an exact field-free part (M_FieldFree, a
  !> bound `FieldFree_Propagate`), which the methods then bind to
  !> `Method_PropagateExactPart`.
  procedure(I_Method_Mb_OrbBased_ApplyLinOp), pointer :: Method_Mb_OrbBased_ApplyLinOp
  abstract interface
    !> Applies the one-body operator propagated by the linear split-step
    !> integrator (full one-body Hamiltonian, or only the gauge coupling with
    !> an exact field-free part) and accumulates the result in `dOrbs`.
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
