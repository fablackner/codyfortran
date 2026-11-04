! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Central two-body interaction API for the CodyFortran framework.
!>
!> @details This module defines the abstract interface for particle-particle
!> interactions (Ŵ operator) used in many-body quantum simulations. It exposes
!> a small, uniform set of procedure pointers that higher layers (e.g., Method,
!> Hamiltonian construction) call to compute and apply interaction potentials.
!>
!> **Physical role:**
!> The interaction module handles the two-body operator Ŵ which appears in the
!> many-body Hamiltonian as:
!>    H = T + V_ext + W
!> where W = (1/2) Σᵢⱼ w(rᵢ,rⱼ) represents electron-electron (or particle-particle)
!> repulsion/attraction.
!>
!> **Architecture:**
!> At runtime, `SysInteraction_Fabricate` reads the JSON configuration and binds
!> these procedure pointers to concrete implementations from one of:
!>   - `Lattice`: Discrete on-site interactions (Hubbard U)
!>   - `Linear`: Real-space convolution (SoftYukawa, contact)
!>   - `Ylm`: Spherical-harmonic expansion (Coulomb via Poisson solvers)
!>   - `Prolate`: Two-center prolate-spheroidal Coulomb solvers
!>   - `Fem`: 3D spectral elements (Galerkin Poisson solve, multipole box values)
!>
!> **Workflow:**
!> 1. `FillInteractionSrc`: Build source term ρ(r) = ψ*ᵢ(r)ψⱼ(r) from orbitals
!> 2. `FillInteractionPotential`: Solve for V(r) = ∫ w(r,r') ρ(r') dr'
!> 3. `MultiplyWithInteractionPotential`: Apply V·ψ to get Ŵψ contribution
!>
!> No interaction logic lives here directly; this module defines the abstract
!> interfaces and public pointers wired up during fabrication, plus the
!> grid-agnostic mean field of a factored one-body density matrix built on
!> them (S_SysInteraction_MeanField): the direct potential
!> `SysInteraction_FillDirectPotential`, the adaptively compressed exchange
!> `SysInteraction_FillExchangeAce`, and their application to a vector,
!> `SysInteraction_AddMeanField`.
module M_SysInteraction
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Runtime factory that reads the user configuration and binds all
    !> `SysInteraction_*` procedure pointers to a concrete implementation.
    !>
    !> Typical responsibilities:
    !> - parse input parameters (e.g. from JSON)
    !> - select the grid/domain back-end (lattice/linear/ylm)
    !> - select the specific interaction model (e.g. on-site, Yukawa, Coulomb)
    !> - set feature flags like `SysInteraction_timeIndependentQ`
    module subroutine SysInteraction_Fabricate
    end subroutine

    !> Direct (Hartree) potential of the one-body density matrix
    !> gamma(r, r') = sum_l weightedOrbs_l(r) conjg(orbs_l(r')) (the convention
    !> of SysInteraction_FillExchangeAce), a grid-agnostic helper on top of the
    !> contract below (S_SysInteraction_MeanField):
    !>   pot = V[factor rho],  rho = sum_l conjg(orbs_l) weightedOrbs_l,
    !> from one potential solve of the total source (the potential is linear
    !> in the source; sources accumulated in the order of the orbitals).
    !> weightedOrbs = orbs with factor_ = 2 gives the closed-shell density
    !> 2 sum |orbs_l|^2 (groundSolver.scf.stdImpl); weightedOrbs = orbs rhoS,
    !> rhoS a Hermitian spatial 1-RDM, that of a correlated state (frozen
    !> mean-field reference of Method, exactFieldFreeMeanFieldQ).
    module subroutine SysInteraction_FillDirectPotential(pot, orbs, weightedOrbs, time, factor_, bt1_, bt2_)
      !> Output: the potential, in the representation of
      !> SysInteraction_FillInteractionPotential.
      complex(R64), intent(out), allocatable :: pot(:)
      !> Orbitals (nPoints x n, n >= 1).
      complex(R64), intent(in), contiguous :: orbs(:, :)
      !> Weighted orbitals of the density matrix gamma (nPoints x n).
      complex(R64), intent(in), contiguous :: weightedOrbs(:, :)
      !> Time of the potential.
      real(R64), intent(in) :: time
      !> Optional factor of the density (default 1).
      real(R64), intent(in), optional :: factor_
      !> Optional body types (target, source) passed to
      !> SysInteraction_FillInteractionPotential.
      integer(I32), intent(in), optional :: bt1_
      integer(I32), intent(in), optional :: bt2_
    end subroutine

    !> Adaptively compressed exchange (ACE, Lin 2016), a grid-agnostic helper
    !> on top of the contract below (S_SysInteraction_MeanField): the
    !> exchange operator of the one-body density matrix
    !> gamma(r, r') = sum_l weightedOrbs_l(r) conjg(orbs_l(r')),
    !>   K v = sum_l V[conjg(orbs_l) v] weightedOrbs_l,
    !> (V[src] = SysInteraction_FillInteractionPotential) is compressed to
    !> K_ACE = sum_k |zeta_k><zeta_k| with xi_i = K orbs_i,
    !> M = <orbs|xi> = L L^H (Grid_InnerProduct, Hermitized) and
    !> zeta = xi L^-H. K_ACE equals K on span(orbs) and is Hermitian and
    !> positive semidefinite in the grid metric; applying it costs one inner
    !> product per vector zeta_k instead of one potential solve per orbital.
    !> weightedOrbs = orbs gives the exchange of closed-shell orbitals of
    !> one spin (groundSolver.scf.stdImpl, aceQ); weightedOrbs = orbs R, R a
    !> Hermitian 1-RDM block, that of a correlated state (frozen mean-field
    !> reference of Method, exactFieldFreeMeanFieldQ).
    !> Cost: n(n+1)/2 potential solves for n orbitals if
    !> SysInteraction_conjSymmetricQ, else n^2.
    module subroutine SysInteraction_FillExchangeAce(zeta, positiveDefiniteQ, orbs, weightedOrbs, time, bt1_, bt2_)
      !> Output: ACE vectors zeta (same shape as orbs, allocated by the caller).
      complex(R64), intent(out), contiguous :: zeta(:, :)
      !> Output: false if M was numerically not positive definite; zeta then
      !> compresses K on the eigenvectors of M above 1e-12 of its largest
      !> eigenvalue (pseudo-inverse of M; the other columns are zero).
      logical, intent(out) :: positiveDefiniteQ
      !> Orbitals spanning the space on which K_ACE = K (nPoints x n).
      complex(R64), intent(in), contiguous :: orbs(:, :)
      !> Weighted orbitals of the density matrix gamma (nPoints x n).
      complex(R64), intent(in), contiguous :: weightedOrbs(:, :)
      !> Time of the interaction potentials.
      real(R64), intent(in) :: time
      !> Optional body types (target, source) passed to
      !> SysInteraction_FillInteractionPotential.
      integer(I32), intent(in), optional :: bt1_
      integer(I32), intent(in), optional :: bt2_
    end subroutine

    !> Applies the mean field J - K_ACE to a vector and accumulates,
    !>   w = w + J z - sum_k zeta_k <zeta_k|z>,
    !> with J = directPotential (SysInteraction_FillDirectPotential) applied by
    !> SysInteraction_MultiplyWithInteractionPotential and the ACE vectors
    !> zeta (SysInteraction_FillExchangeAce): n inner products, one
    !> multiplication and one fused pass, per element
    !> ((w + J z) - o_1 zeta_1) - o_2 zeta_2 - ... Reentrant; parallel when
    !> called outside a parallel region (S_SysInteraction_MeanField).
    module subroutine SysInteraction_AddMeanField(w, scratch, z, directPotential, zeta)
      !> In/out: vector to accumulate into (nPoints).
      complex(R64), intent(inout), contiguous :: w(:)
      !> Work vector (nPoints), overwritten with J z.
      complex(R64), intent(inout), contiguous :: scratch(:)
      !> Vector the mean field acts on (nPoints).
      complex(R64), intent(in), contiguous :: z(:)
      !> Direct potential J.
      complex(R64), intent(in), contiguous :: directPotential(:)
      !> ACE vectors (nPoints x n; n = 0: J only).
      complex(R64), intent(in), contiguous :: zeta(:, :)
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> True if the interaction is independent of body types (species). When set,
  !> callers may cache results across body-type combinations.
  logical :: SysInteraction_bodyTypeIndependentQ = .false.

  !> True if the interaction is time independent. When set, time arguments in
  !> back-end routines can be ignored and results cached/reused.
  logical :: SysInteraction_timeIndependentQ = .false.

  !> True if the back-end guarantees that swapping the source orbital pair
  !> conjugates the potential in the stored representation:
  !>   V[conjg(phi_j)*phi_i] = ConjugateInteractionPotential(V[conjg(phi_i)*phi_j])
  !> This holds for real, symmetric interaction kernels with real quadrature
  !> weights (e.g., SoftYukawa, on-site Hubbard, Coulomb on a real radial
  !> contour) and lets callers halve the number of potential solves over
  !> orbital pairs. Must stay false for complex-scaled (ECS) contours.
  logical :: SysInteraction_conjSymmetricQ = .false.

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> Pointer to a cheap, idempotent setup routine for the selected interaction
  !> back-end. Use this to allocate caches, precompute kernels, etc.
  procedure(I_SysInteraction_Setup), pointer :: SysInteraction_Setup => NoOpProcedures_Setup
  abstract interface
    !> Initializes the selected interaction implementation (allocate tables,
    !> FFT plans, masks, boundary helpers, ...). Should be safe to call
    !> multiple times.
    subroutine I_SysInteraction_Setup
    end subroutine
  end interface

  !> Pointer to a routine that builds the full interaction potential resulting
  !> from the provided source term (accumulated over all source grid points).
  procedure(I_SysInteraction_FillInteractionPotential), pointer :: SysInteraction_FillInteractionPotential

  abstract interface
    subroutine I_SysInteraction_FillInteractionPotential(interactionPotential, src, time, bt1_, bt2_)
      import :: I32, R64
      !> Output: accumulated interaction potential on the active grid.
      complex(R64), intent(out), allocatable :: interactionPotential(:)
      !> Input: complex source term sampled on grid points.
      complex(R64), intent(in), contiguous  :: src(:)
      !> Physical time associated with the source/potential.
      real(R64), intent(in)     :: time
      !> Optional body type (species) index for the potential target.
      integer(I32), intent(in), optional :: bt1_
      !> Optional body type (species) index for the source.
      integer(I32), intent(in), optional :: bt2_
    end subroutine
  end interface

  !> Pointer to a routine that constructs a source term from a pair of orbitals
  !> (typically a bra-ket product) on the active grid.
  procedure(I_SysInteraction_FillInteractionSrc), pointer :: SysInteraction_FillInteractionSrc
  abstract interface
    !> Build a complex source term from two orbitals (e.g., rho = conjg(psi_i)*psi_j).
    subroutine I_SysInteraction_FillInteractionSrc(src, orbConjg, orb)
      import :: R64
      !> Output: complex source term over grid points.
      complex(R64), intent(out), allocatable :: src(:)
      !> Input: conjugated orbital (bra).
      complex(R64), intent(in), contiguous :: orbConjg(:)
      !> Input: orbital (ket).
      complex(R64), intent(in), contiguous :: orb(:)
    end subroutine
  end interface

  !> Pointer to a routine that returns the potential belonging to the complex
  !> conjugated source, in the grid's stored representation. Only valid to call
  !> when `SysInteraction_conjSymmetricQ` is true.
  procedure(I_SysInteraction_ConjugateInteractionPotential), pointer :: SysInteraction_ConjugateInteractionPotential
  abstract interface
    !> Fill `potOut` with the representation of conjg(potIn as a spatial function).
    subroutine I_SysInteraction_ConjugateInteractionPotential(potOut, potIn)
      import :: R64
      !> Output: conjugated interaction potential.
      complex(R64), intent(out), contiguous :: potOut(:)
      !> Input: interaction potential to conjugate.
      complex(R64), intent(in), contiguous  :: potIn(:)
    end subroutine
  end interface

  !> Pointer to a routine that applies the interaction potential to an orbital
  !> (point-wise multiplication on the active grid).
  procedure(I_SysInteraction_MultiplyWithInteractionPotential), pointer :: SysInteraction_MultiplyWithInteractionPotential
  abstract interface
    !> Multiply `orb` by `interactionPotential` element-wise and return the result.
    subroutine I_SysInteraction_MultiplyWithInteractionPotential(resOrb, interactionPotential, orb)
      import :: R64
      !> Output: resulting orbital after multiplication.
      complex(R64), intent(out), contiguous :: resOrb(:)
      !> Input: interaction potential over grid points.
      complex(R64), intent(in), contiguous  :: interactionPotential(:)
      !> Input: orbital to be multiplied.
      complex(R64), intent(in), contiguous  :: orb(:)
    end subroutine
  end interface

end module
