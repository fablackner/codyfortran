! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Central coefficient API used throughout CodyFortran.
!>
!> This module defines the public contract for CI (Configuration Interaction)
!> coefficients and exposes it via procedure pointers. At runtime these pointers
!> are wired to representation-specific implementations (generic, Hubbard, …)
!> by `Coeffs_Fabricate`. No numerical work happens here; this is the interface
!> layer that other parts of the code call into.
module M_Coeffs
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Dynamically assigns the coefficient-related procedures at runtime and
    !> sets up the coefficient system based on input parameters.
    module subroutine Coeffs_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Total number of CI coefficients in the current representation.
  integer(I32) :: Coeffs_nCoeffs

  !> Pointer to the coefficients segment of the global state vector.
  !> `Coeffs_coeffs(i)` is the i-th CI coefficient in the active basis.
  !> Associated only for representations that explicitly store CI coefficients.
  complex(R64), contiguous, pointer :: Coeffs_coeffs(:)

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> Pointer to the representation-specific setup routine.
  procedure(I_Coeffs_Setup), pointer :: Coeffs_Setup => NoOpProcedures_Setup
  abstract interface
    !> Initialize internal data needed by the active coefficient representation.
    !>
    !> Must be called after `Coeffs_Fabricate` has attached the concrete
    !> implementation and before any operation that depends on representation
    !> internals (e.g., mapping, precomputed tables).
    subroutine I_Coeffs_Setup
    end subroutine
  end interface

  !> Apply the one-body Hamiltonian and optionally accumulate the one-body RDM.
  procedure(I_Coeffs_ApplyH1FillRdm1), pointer :: Coeffs_ApplyH1FillRdm1
  abstract interface
    !> Apply the one-body Hamiltonian to `coeffs` and optionally compute the
    !> one-body reduced density matrix (RDM1).
    !> \[
    !> H1 = \sum_{i,j} \langle j|H1|i \rangle a^\dagger_j a_i
    !> \]
    subroutine I_Coeffs_ApplyH1FillRdm1(coeffs, rdm1_, dCoeffs_, h1_)
      import :: R64
      !> Input CI coefficients representing the current quantum state.
      complex(R64), intent(in), contiguous, target :: coeffs(:)
      !> on presence, allocate and fill RDM1 with the convention
      !> RDM1(p,q) = <Psi| a^\dagger_p a_q |Psi> consistent with the active basis.
      complex(R64), intent(out), allocatable, optional :: rdm1_(:, :)
      !> in-place accumulator for H1|coeffs>. If present, the result is
      !> added to `dCoeffs_` (AXPY-style). Caller must pre-initialize the array.
      complex(R64), intent(inout), contiguous, target, optional :: dCoeffs_(:)
      !> explicit one-body Hamiltonian matrix in the active basis.
      !> If absent, the representation provides H1 implicitly.
      complex(R64), intent(in), contiguous, optional :: h1_(:, :)
    end subroutine
  end interface

  !> Apply the two-body Hamiltonian and optionally accumulate the two-body RDM.
  procedure(I_Coeffs_ApplyH2FillRdm2), pointer :: Coeffs_ApplyH2FillRdm2
  abstract interface
    !> Apply the two-body Hamiltonian to `coeffs` and optionally compute the
    !> two-body reduced density matrix (RDM2).
    !> \[
    !> H2 = \sum_{i1,i2,j1,j2} \langle j1,j2|H2|i1,i2 \rangle a^\dagger_{j1} a^\dagger_{j2} a_{i2} a_{i1}
    !> \]
    subroutine I_Coeffs_ApplyH2FillRdm2(coeffs, rdm2_, dCoeffs_, h2_)
      import :: R64
      !> Input CI coefficients representing the current quantum state.
      complex(R64), intent(in), contiguous, target :: coeffs(:)
      !> on presence, allocate and fill RDM2 with the convention
      !> RDM2(p,q,r,s) = <Psi| a^\dagger_p a^\dagger_q a_s a_r |Psi>.
      complex(R64), intent(out), allocatable, optional :: rdm2_(:, :, :, :)
      !> in-place accumulator for H2|coeffs>. If present, the result is
      !> added to `dCoeffs_` (AXPY-style). Caller must pre-initialize the array.
      complex(R64), intent(inout), contiguous, target, optional :: dCoeffs_(:)
      !> explicit two-body Hamiltonian tensor in the active basis.
      !> If absent, the representation provides H2 implicitly.
      complex(R64), intent(in), contiguous, optional :: h2_(:, :, :, :)
    end subroutine
  end interface

  !> Compute a three-body reduced density matrix for selected body types.
  procedure(I_Coeffs_FillRdm3Bt), pointer :: Coeffs_FillRdm3Bt
  abstract interface
    !> Compute the three-body reduced density matrix for specific body types
    !> (e.g., spin species). Body-type conventions are representation-specific.
    !> \[
    !> \rho_{p,q,r,s,t,u} = \langle \Psi | a^\dagger_p a^\dagger_q a^\dagger_r a_u a_t a_s | \Psi \rangle
    !> \]
    subroutine I_Coeffs_FillRdm3Bt(rdm3Bt, coeffs, bt1, bt2, bt3)
      import :: I32, R64
      !> On exit: allocated 6D array holding the three-body RDM block.
      complex(R64), intent(out), allocatable :: rdm3Bt(:, :, :, :, :, :)
      !> Input CI coefficients representing the current quantum state.
      complex(R64), intent(in), contiguous, target :: coeffs(:)
      !> First body type index.
      integer(I32), intent(in) :: bt1
      !> Second body type index.
      integer(I32), intent(in) :: bt2
      !> Third body type index.
      integer(I32), intent(in) :: bt3
    end subroutine
  end interface

  !> Compute a two-body reduced density matrix for selected body types.
  procedure(I_Coeffs_FillRdm2Bt), pointer :: Coeffs_FillRdm2Bt
  abstract interface
    !> Compute the two-body reduced density matrix for specific body types
    !> (e.g., spin species). Body-type conventions are representation-specific.
    !> \[
    !> \rho_{p,q,r,s} = \langle \Psi | a^\dagger_p a^\dagger_q a_s a_r | \Psi \rangle
    !> \]
    subroutine I_Coeffs_FillRdm2Bt(rdm2Bt, coeffs, bt1, bt2)
      import :: I32, R64
      !> On exit: allocated 4D array holding the two-body RDM block.
      complex(R64), intent(out), allocatable :: rdm2Bt(:, :, :, :)
      !> Input CI coefficients representing the current quantum state.
      complex(R64), intent(in), contiguous, target :: coeffs(:)
      !> First body type index.
      integer(I32), intent(in) :: bt1
      !> Second body type index.
      integer(I32), intent(in) :: bt2
    end subroutine
  end interface

  !> Compute a one-body reduced density matrix for a selected body type.
  procedure(I_Coeffs_FillRdm1Bt), pointer :: Coeffs_FillRdm1Bt
  abstract interface
    !> Compute the one-body reduced density matrix for a specific body type
    !> (e.g., spin species). Body-type conventions are representation-specific.
    !> \[
    !> \rho_{p,q} = \langle \Psi | a^\dagger_p a_q | \Psi \rangle
    !> \]
    subroutine I_Coeffs_FillRdm1Bt(rdm1Bt, coeffs, bt1)
      import :: I32, R64
      !> On exit: allocated 2D array holding the one-body RDM block.
      complex(R64), intent(out), allocatable :: rdm1Bt(:, :)
      !> Input CI coefficients representing the current quantum state.
      complex(R64), intent(in), contiguous, target :: coeffs(:)
      !> Body type index.
      integer(I32), intent(in) :: bt1
    end subroutine
  end interface

  !> Normalize the CI coefficient vector.
  procedure(I_Coeffs_Normalize), pointer :: Coeffs_Normalize
  abstract interface
    !> Normalize `coeffs` in-place such that
    !> \[
    !> \langle \Psi | \Psi \rangle = 1
    !> \]
    subroutine I_Coeffs_Normalize(coeffs)
      import :: R64
      !> In/out: CI coefficients. On exit the vector is normalized.
      complex(R64), intent(inout), contiguous :: coeffs(:)
    end subroutine
  end interface

  !> Orthogonalize a vector against a given state (single-vector projection).
  procedure(I_Coeffs_ProjectOnSubspace), pointer :: Coeffs_ProjectOnSubspace
  abstract interface
    !> Project `dCoeffs` onto the subspace orthogonal to `coeffs`:
    !> \[
    !> dCoeffs \leftarrow dCoeffs - \langle coeffs | dCoeffs \rangle\, coeffs
    !> \]
    subroutine I_Coeffs_ProjectOnSubspace(dCoeffs, coeffs)
      import :: R64
      !> In/out: vector to be projected. Modified in-place.
      complex(R64), intent(inout), contiguous :: dCoeffs(:)
      !> Reference state defining the 1D subspace to remove.
      complex(R64), intent(in), contiguous :: coeffs(:)
    end subroutine
  end interface

  !> Apply a sequence of creation/annihilation operators to the CI state.
  procedure(I_Coeffs_ApplyExcitation), pointer :: Coeffs_ApplyExcitation
  abstract interface
    !> Transform the CI coefficients by applying the excitation defined by
    !> `creates` and `destroys` on the selected body type.
    subroutine I_Coeffs_ApplyExcitation(coeffs, creates, destroys, bt)
      import :: I32, R64
      !> In/out: CI coefficients. Modified in-place with the excitation applied.
      complex(R64), intent(inout), contiguous :: coeffs(:)
      !> Indices of orbitals where particles are created (a^\dagger operators).
      integer(I32), intent(in), contiguous :: creates(:)
      !> Indices of orbitals where particles are annihilated (a operators).
      integer(I32), intent(in), contiguous :: destroys(:)
      !> Body type index specifying which particle species to operate on.
      integer(I32), intent(in) :: bt
    end subroutine
  end interface

  !> Decode a linear CI index into a configurations array.
  procedure(I_Coeffs_ConfigurationsFromIndex), pointer :: Coeffs_ConfigurationsFromIndex
  abstract interface
    !> Convert a 1-based CI coefficient index into the representation-specific
    !> configurations array (occupation pattern, bit fields, …).
    pure subroutine I_Coeffs_ConfigurationsFromIndex(configurations, iCoeff)
      import :: I32
      !> On exit: representation-specific configuration descriptor.
      integer(I32), intent(out), contiguous :: configurations(:)
      !> 1-based CI coefficient index to decode.
      integer(I32), intent(in) :: iCoeff
    end subroutine
  end interface

  !> Encode a configurations array into a linear CI index.
  procedure(I_Coeffs_IndexFromConfigurations), pointer :: Coeffs_IndexFromConfigurations
  abstract interface
    !> Convert a representation-specific configurations array (occupation pattern,
    !> bit fields, …) into the corresponding 1-based CI coefficient index.
    pure subroutine I_Coeffs_IndexFromConfigurations(iCoeff, configurations)
      import :: I32
      !> On exit: 1-based CI coefficient index.
      integer(I32), intent(out) :: iCoeff
      !> Input: representation-specific configuration descriptor.
      integer(I32), intent(in), contiguous :: configurations(:)
    end subroutine
  end interface

  !> Persist CI coefficients to storage (format is representation-specific).
  procedure(I_Coeffs_SaveCoeffs), pointer :: Coeffs_SaveCoeffs
  abstract interface
    !> Save the CI coefficient vector to file(s). The concrete implementation
    !> decides the file format and destination based on global I/O settings.
    subroutine I_Coeffs_SaveCoeffs(coeffs)
      import :: R64
      !> Input CI coefficients to persist.
      complex(R64), intent(in), contiguous :: coeffs(:)
    end subroutine
  end interface

  !> Compute and persist a two-body RDM derived from CI coefficients.
  procedure(I_Coeffs_SaveTwoRdm), pointer :: Coeffs_SaveTwoRdm
  abstract interface
    !> Compute the two-body RDM from the provided coefficients and write it to
    !> file(s). The concrete implementation defines the output format.
    subroutine I_Coeffs_SaveTwoRdm(coeffs)
      import :: R64
      !> Input CI coefficients used to compute the RDM.
      complex(R64), intent(in), contiguous :: coeffs(:)
    end subroutine
  end interface

end module

