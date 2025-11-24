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
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Classification of basis states by body type to enable type-dependent
  !> potentials and interactions (e.g., spin-dependent terms).
  !> Method_Mb_OrbBased_bodyTypeOfOrb(i) gives the body type of basis state i.
  integer(I32), allocatable :: Method_Mb_OrbBased_bodyTypeOfOrb(:)

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
    subroutine I_Method_Mb_OrbBased_ApplyCorrelationOp(dOrbs, orbs, rdm1, rdm2, time, h2_)
      import :: R64
      !> Input/output derivative orbitals. On exit, contains the original
      !> content plus the result of the correlation operator action.
      complex(R64), intent(inout), contiguous, target :: dOrbs(:, :)
      !> Input orbitals to which the operator is applied.
      complex(R64), intent(in), contiguous, target :: orbs(:, :)
      !> One-body reduced density matrix.
      complex(R64), intent(in), contiguous :: rdm1(:, :)
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

  !> Pointer to the procedure for computing the time derivative with static orbitals.
  procedure(I_Method_Mb_OrbBased_TimeDerivativeOrbsLin), pointer :: &
    Method_Mb_OrbBased_TimeDerivativeOrbsLin
  abstract interface
    !> Computes the right-hand side of the TDSE with static orbitals,
    !> returning `-i * H(state)` in `dState`.
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

  !> Pointer to the procedure for computing the non-linear part of the time derivative
  !> when both coefficients and orbitals evolve (e.g., in MCTDHX-like schemes).
  procedure(I_Method_Mb_OrbBased_TimeDerivativeOrbsLin), pointer :: &
    Method_Mb_OrbBased_TimeDerivativeCoeffsPlusOrbsNonLin

end module
