! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Exact propagation of the orbitals with the field-free one-body Hamiltonian
!> h0 = T + V_ext: the exactly solved part of split-step and exponential
!> propagators of the orbital-based methods (Method_PropagateExactPart).
!>
!> With a "fieldFree" block the linear split-step integrator of the method
!> only propagates what is left of the one-body Hamiltonian (the gauge
!> coupling), and every linear step is composed with FieldFree_Propagate
!> (propagator.splitStep.exactPartScheme). Without the block
!> FieldFree_Fabricate leaves FieldFree_Propagate null and the integrator
!> propagates the whole one-body part. A bound FieldFree_Propagate is the
!> switch the methods read.
!>
!> Variants ("fieldFree.<variant>"):
!>   - eigen:  dense eigen-decomposition of the h0 blocks per l on Ylm grids,
!>             exp(-i h0 dt) and exp(-h0 dt); phi functions for exponential
!>             integrators (FieldFree_ApplyPhi) and a spectral cut
!>             (FieldFree_ProjectBelowCut).
!>   - cayley: Crank-Nicolson (Cayley) step by preconditioned GMRES on grids
!>             whose kinetic backend binds SysKinetic_SolveShifted (grid.fem),
!>             real time only; optional frozen mean-field reference U, so that
!>             h0 + U is propagated (FieldFree_meanFieldQ).
!>
!> Order: FieldFree_Fabricate before Method_Fabricate (the method binds its
!> exact part to what is bound here), also for an input without the block;
!> FieldFree_Setup after all fabricates and the setups of Grid, SysKinetic
!> and SysPotential, before the first propagation.
module M_FieldFree
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Resets the pointers and flags of this module and dispatches to the
    !> variant of the "fieldFree" block; without the block nothing is bound.
    module subroutine FieldFree_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Whether a spectral cut is set (fieldFree.eigen.cut): the orbitals are
  !> restricted to the h0 eigenvectors below the cut, FieldFree_Propagate
  !> drops the others and the method projects its orbital derivatives with
  !> FieldFree_ProjectBelowCut.
  logical :: FieldFree_cutQ = .false.

  !> Whether the exact part propagates h0 + U with a frozen mean-field
  !> reference U (fieldFree.cayley.meanFieldQ; static, Hermitian, built once
  !> from the state, FieldFree_BuildMeanField). The method then subtracts U in
  !> its nonlinear part (FieldFree_ApplyMeanField), so the equations of motion
  !> are unchanged: the split is that of the Fock (Kohn-Sham) operator of the
  !> reference and a remainder that is small near it.
  logical :: FieldFree_meanFieldQ = .false.

  !> Whether the frozen mean-field reference is built (since the last
  !> FieldFree_Setup, which releases it).
  logical :: FieldFree_meanFieldBuiltQ = .false.

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> Pointer to the setup of the exact field-free part (operators, caches,
  !> checks); a no-op without a "fieldFree" block.
  procedure(I_FieldFree_Setup), pointer :: FieldFree_Setup => NoOpProcedures_Setup
  abstract interface
    !> Prepares the exact field-free part.
    subroutine I_FieldFree_Setup
    end subroutine
  end interface

  !> Pointer to the exact field-free propagation of the orbitals,
  !> orbs <- U0(dt) orbs: exp(-i h0 dt) by the h0 eigen-decomposition (eigen;
  !> also exp(-h0 dt) in imaginary time) or the Cayley approximation
  !> (1 + i dt h0/2)^-1 (1 - i dt h0/2) (cayley, real time only); null without
  !> a "fieldFree" block.
  procedure(I_FieldFree_Propagate), pointer :: FieldFree_Propagate => null()
  abstract interface
    !> Propagates the orbitals with the exact field-free part over dt.
    subroutine I_FieldFree_Propagate(orbs, dt, imagTimeQ)
      import :: R64
      !> Orbitals (nPoints x nOrbsInState), propagated in place.
      complex(R64), intent(inout), contiguous :: orbs(:, :)
      !> Time step (may be negative).
      real(R64), intent(in) :: dt
      !> Imaginary-time propagation exp(-h0 dt).
      logical, intent(in) :: imagTimeQ
    end subroutine
  end interface

  !> Optional hook: the phi functions of the exact field-free part,
  !> orbs = sum_m phi_m(-i h0 dt) orbsList(:, :, m), for exponential
  !> integrators; null unless the variant provides them (eigen).
  procedure(I_FieldFree_ApplyPhi), pointer :: FieldFree_ApplyPhi => null()
  abstract interface
    !> Phi functions of h0 applied to a list of orbital sets.
    subroutine I_FieldFree_ApplyPhi(orbs, orbsList, dt)
      import :: R64
      !> Output: orbitals (nPoints x nOrbsInState).
      complex(R64), intent(out), contiguous :: orbs(:, :)
      !> Input: orbitals (nPoints x nOrbsInState x 0:K), slice m multiplies phi_m.
      complex(R64), intent(in), contiguous :: orbsList(:, :, 0:)
      !> Time step.
      real(R64), intent(in) :: dt
    end subroutine
  end interface

  !> Pointer to the projection of orbitals (or orbital derivatives) onto the
  !> h0 eigenvectors below the spectral cut; a no-op without a cut.
  procedure(I_FieldFree_ProjectBelowCut), pointer :: FieldFree_ProjectBelowCut => NoOpProcedures_ApplyAbsorber
  abstract interface
    !> Projects the orbitals in place onto the space kept by the spectral cut.
    subroutine I_FieldFree_ProjectBelowCut(orbs)
      import :: R64
      !> In/out: orbitals (nPoints x nOrbsInState), projected in place.
      complex(R64), intent(inout), contiguous :: orbs(:, :)
    end subroutine
  end interface

  !> Optional hook: builds the frozen mean-field reference U of the exact part
  !> (replaces one built earlier): U = J0 - K0 from the orbitals and the 1-RDM
  !> of a state, or U = v_KS[rho0] - a K0 for a Kohn-Sham method
  !> (Method_kohnShamQ). From then on FieldFree_Propagate propagates h0 + U.
  !> The methods build it from the state at their first exact-part (or
  !> nonlinear) call after the setups unless a program called this before
  !> (for a kick: from the ground state before the kick, which leaves no
  !> static part in the remainder; continued runs: from the same state again;
  !> A_FieldFree.md, "Which state"). Null unless the variant provides a
  !> reference (cayley).
  procedure(I_FieldFree_BuildMeanField), pointer :: FieldFree_BuildMeanField => null()
  abstract interface
    !> Builds the frozen mean-field reference from a state.
    subroutine I_FieldFree_BuildMeanField(orbs, rdm1, time)
      import :: R64
      !> Stored orbitals of the state (nPoints x nOrbsInState).
      complex(R64), intent(in), contiguous :: orbs(:, :)
      !> One-body reduced density matrix of the state (nOrbsSum x nOrbsSum,
      !> rdm1(i, j) = <a+_j a_i>).
      complex(R64), intent(in), contiguous :: rdm1(:, :)
      !> Time of the state (interaction potentials).
      real(R64), intent(in) :: time
    end subroutine
  end interface

  !> Optional hook: dOrbs(:, i) += U orbs(:, i) with the frozen mean-field
  !> reference U of the body type of orbital i; with h1_ also its matrix
  !> h1_(i, j) = <orbs_i|U orbs_j> (equal body types; nOrbsSum x nOrbsSum,
  !> restricted: both spin blocks), the convention of
  !> Method_Mb_OrbBased_FillH1. What the method subtracts in its nonlinear
  !> part. Error stop if the reference is not built. Null unless the variant
  !> provides a reference (cayley).
  procedure(I_FieldFree_ApplyMeanField), pointer :: FieldFree_ApplyMeanField => null()
  abstract interface
    !> Applies the frozen mean-field reference to orbitals.
    subroutine I_FieldFree_ApplyMeanField(dOrbs, orbs, h1_)
      import :: R64
      !> In/out: orbitals to accumulate into (nPoints x nOrbsInState).
      complex(R64), intent(inout), contiguous :: dOrbs(:, :)
      !> Input: orbitals the reference acts on (nPoints x nOrbsInState).
      complex(R64), intent(in), contiguous :: orbs(:, :)
      !> Optional output: matrix of U in the orbitals.
      complex(R64), intent(out), allocatable, optional :: h1_(:, :)
    end subroutine
  end interface

end module
