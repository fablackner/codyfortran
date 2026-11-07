! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Exchange-correlation (density functional) interface module.
!>
!> @details
!> The density-functional approximation of the electron-electron interaction
!> of Kohn-Sham DFT and TDDFT: the local Kohn-Sham potential of the spin
!> densities,
!>
!>   v_KS,s = J[ρ] + v_xc,s[ρ↑, ρ↓],   ρ = ρ↑ + ρ↓,
!>
!> with the Hartree potential J of SysInteraction and the exchange-correlation
!> potential v_xc = δE_xc/δρ_s of a functional (Libxc), plus the fraction
!> `SysXc_exactExchangeFraction` of exact exchange of a hybrid functional,
!> which the consumers (method.mb.orbBased.tdks, groundSolver.scf.stdImpl)
!> apply with the exchange of SysInteraction.
!>
!> ## Xc points
!>
!> Densities and local potentials are real arrays on the `SysXc_nPoints` xc
!> points of the grid family, with quadrature weights `SysXc_weights`
!> (Σ_p w_p f_p ≈ ∫ f d^d r): the nodes of grid.fem and grid.linear, an
!> angular product quadrature of every radial point of grid.ylm. The family
!> binds the maps between orbitals and xc points (density, application of a
!> local potential), between densities and the source/potential layout of
!> SysInteraction (Hartree potential), and for gradient functionals the
!> gradient G of a field on the xc points and its adjoint W⁻¹GᵀW in the
!> quadrature metric.
!>
!> ## Discrete functional and its derivative
!>
!> E_xc = Σ_p w_p e(ρ_p, σ_p) with σ = |Gρ|² (GGA), and the potential is the
!> exact derivative of this discrete functional with respect to the density
!> values, v = (1/w_p) ∂E_xc/∂ρ_p:
!>
!>   v_s = ∂e/∂ρ_s + W⁻¹Gᵀ W (2 ∂e/∂σ_ss Gρ_s + ∂e/∂σ_ud Gρ_s'),
!>
!> so the Kohn-Sham orbital equations are the Hamiltonian flow of the discrete
!> Kohn-Sham energy for any G (TDDFT conserves it; the SCF minimizes it).
!>
!> ## Architecture
!>
!>   - **Interface Module** (`M_SysXc.f90`): contract, data, pointers
!>   - **Dispatch Submodule** (`S_SysXc.f90`): JSON routing to the grid
!>     family, defaults of the nodal families (Fem, Linear)
!>   - **Kohn-Sham Submodule** (`S_SysXc_KohnSham.f90`): the grid-agnostic
!>     drivers on the contract (xc, Hartree and Kohn-Sham potentials)
!>   - **Families** (`Fem/`, `Ylm/`, `Linear/`): xc points and maps
!>   - **Functional** (`Libxc/`): evaluation on the xc points, shared by all
!>     families (`sysXc.<family>.libxc`)
!>
!> ## Usage
!>
!>   call SysXc_Fabricate      ! after SysInteraction_Fabricate
!>   call SysXc_Setup          ! after Grid_Setup and SysInteraction_Setup
!>   density = 0
!>   call SysXc_AddDensity(density(:, s), orbs, weightedOrbs)   ! per spin
!>   call SysXc_FillKohnShamPotential(ksPot, eHartree, eXc, density, time)
!>   call SysXc_MultiplyWithPotential(dOrb, ksPot(:, s), orb)
!>
!> @see M_SysInteraction, M_SysXc_Libxc, M_Utils_XcLib
module M_SysXc
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> @brief Bind the xc family and functional from `sysXc.*`.
    !>
    !> @details
    !> Binds the nodal defaults, dispatches to the grid family
    !> (`sysXc.fem`, `sysXc.ylm`, `sysXc.linear`), which binds its maps and
    !> dispatches to the functional (`sysXc.<family>.libxc`). Fabricate
    !> SysInteraction first: the Hartree potential uses its layout.
    module subroutine SysXc_Fabricate
    end subroutine

    !> @brief Exchange-correlation potentials and energy of spin densities.
    !>
    !> @details
    !> xcPot(:, s) = v_xc,s = (1/w_p) ∂E_xc/∂ρ_s,p and E_xc = Σ_p w_p e_p of
    !> the densities on the xc points (nSpin = size(density, 2): 1 for the
    !> total density of a spin-unpolarized state, 2 for ρ↑, ρ↓). Gradient
    !> functionals use `SysXc_ApplyGradient` and `SysXc_AddGradientAdjoint`
    !> (see the module header). Negative densities (round-off) count as zero.
    module subroutine SysXc_FillXcPotential(xcPot, xcEnergy, density)
      !> v_xc per spin channel, (SysXc_nPoints, nSpin)
      real(R64), intent(out), contiguous :: xcPot(:, :)
      !> E_xc
      real(R64), intent(out) :: xcEnergy
      !> Densities on the xc points, (SysXc_nPoints, nSpin)
      real(R64), intent(in), contiguous :: density(:, :)
    end subroutine

    !> @brief Hartree potential and energy of a total density.
    !>
    !> @details
    !> hartreePot = J[ρ] on the xc points from one interaction-potential solve
    !> of the source of ρ (SysXc_FillInteractionSrc,
    !> SysInteraction_FillInteractionPotential, SysXc_AddInteractionPotential)
    !> and E_H = ½ Σ_p w_p ρ_p J_p; zero without a `sysInteraction` block
    !> (`SysXc_hartreeQ`).
    module subroutine SysXc_FillHartreePotential(hartreePot, hartreeEnergy, density, time)
      !> J on the xc points (SysXc_nPoints)
      real(R64), intent(out), contiguous :: hartreePot(:)
      !> E_H
      real(R64), intent(out) :: hartreeEnergy
      !> Total density on the xc points (SysXc_nPoints)
      real(R64), intent(in), contiguous :: density(:)
      !> Time of the interaction potential
      real(R64), intent(in) :: time
    end subroutine

    !> @brief Local Kohn-Sham potentials v_KS,s = J + v_xc,s and their energies.
    !>
    !> @details
    !> The Hartree potential of the total density Σ_s density(:, s) and the
    !> xc potentials of the spin densities (SysXc_FillHartreePotential,
    !> SysXc_FillXcPotential); the exact exchange of a hybrid is not included.
    module subroutine SysXc_FillKohnShamPotential(kohnShamPot, hartreeEnergy, xcEnergy, density, time)
      !> v_KS per spin channel, (SysXc_nPoints, nSpin)
      real(R64), intent(out), contiguous :: kohnShamPot(:, :)
      !> E_H
      real(R64), intent(out) :: hartreeEnergy
      !> E_xc
      real(R64), intent(out) :: xcEnergy
      !> Densities on the xc points, (SysXc_nPoints, nSpin)
      real(R64), intent(in), contiguous :: density(:, :)
      !> Time of the interaction potential
      real(R64), intent(in) :: time
    end subroutine

    !> @brief Quadrature Σ_p w_p f_p of a field on the xc points.
    module function SysXc_Integrate(f) result(res)
      real(R64), intent(in), contiguous :: f(:)
      real(R64) :: res
    end function
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Number of xc points (set by the family at fabricate or setup)
  integer(I32) :: SysXc_nPoints = 0

  !> Smallest work (array elements) of an OpenMP-parallel loop of SysXc:
  !> smaller loops run serially (the fork and join would cost more). Loops
  !> around BLAS calls are never parallel (the GEMMs thread themselves).
  integer(I32), parameter :: SysXc_PARALLEL_MIN_WORK = 32768

  !> Quadrature weights of the xc points (SysXc_nPoints; SysXc_Setup)
  real(R64), allocatable :: SysXc_weights(:)

  !> Spatial dimension of the grid (functionals of another dimension are refused)
  integer(I32) :: SysXc_dimension = 3

  !> Whether the functional needs density gradients (GGA)
  logical :: SysXc_gradientQ = .false.

  !> Fraction of exact (Hartree-Fock) exchange of a hybrid functional, applied
  !> by the consumers with SysInteraction (0 for semilocal functionals)
  real(R64) :: SysXc_exactExchangeFraction = 0.0_R64

  !> Whether the Kohn-Sham potential contains the Hartree potential (a
  !> `sysInteraction` block exists)
  logical :: SysXc_hartreeQ = .false.

  !=============================================================================
  ! module procedures pointers: grid family
  !=============================================================================

  !> @brief Setup of the family: xc points, weights and tables (after Grid_Setup).
  procedure(I_SysXc_Setup), pointer :: SysXc_Setup => NoOpProcedures_Setup
  abstract interface
    subroutine I_SysXc_Setup
    end subroutine
  end interface

  !> @brief Density of a factored one-body density matrix on the xc points.
  procedure(I_SysXc_AddDensity), pointer :: SysXc_AddDensity
  abstract interface
    !> density += Re Σ_l conj(φ_l) g_l at the xc points, the density of
    !> γ(r, r') = Σ_l g_l(r) φ_l*(r') (orbitals φ = orbs, g = weightedOrbs =
    !> orbs R with a Hermitian matrix R; R = 1 for occupied orbitals).
    subroutine I_SysXc_AddDensity(density, orbs, weightedOrbs)
      import :: R64
      !> Density on the xc points (SysXc_nPoints), accumulated
      real(R64), intent(inout), contiguous, target :: density(:)
      !> Orbitals (Grid_nPoints, n)
      complex(R64), intent(in), contiguous, target :: orbs(:, :)
      !> Weighted orbitals (Grid_nPoints, n)
      complex(R64), intent(in), contiguous, target :: weightedOrbs(:, :)
    end subroutine
  end interface

  !> @brief Application of a local potential on the xc points to an orbital.
  procedure(I_SysXc_MultiplyWithPotential), pointer :: SysXc_MultiplyWithPotential
  abstract interface
    !> resOrb = v·orb in the Galerkin sense of the grid (Hermitian in its
    !> metric, ⟨a|v b⟩ = Σ_p w_p conj(a_p) v_p b_p at the xc points); parallel
    !> inside unless called in an active parallel region.
    subroutine I_SysXc_MultiplyWithPotential(resOrb, pot, orb)
      import :: R64
      !> v·orb (Grid_nPoints)
      complex(R64), intent(out), contiguous, target :: resOrb(:)
      !> Local potential on the xc points (SysXc_nPoints)
      real(R64), intent(in), contiguous, target :: pot(:)
      !> Orbital (Grid_nPoints)
      complex(R64), intent(in), contiguous, target :: orb(:)
    end subroutine
  end interface

  !> @brief Interaction source of a density in the layout of SysInteraction.
  procedure(I_SysXc_FillInteractionSrc), pointer :: SysXc_FillInteractionSrc
  abstract interface
    !> src as SysInteraction_FillInteractionSrc would build it for a pair
    !> density equal to density (quadrature weights included as there).
    subroutine I_SysXc_FillInteractionSrc(src, density)
      import :: R64
      complex(R64), intent(out), allocatable, target :: src(:)
      real(R64), intent(in), contiguous, target :: density(:)
    end subroutine
  end interface

  !> @brief Interaction potential (SysInteraction layout) on the xc points.
  procedure(I_SysXc_AddInteractionPotential), pointer :: SysXc_AddInteractionPotential
  abstract interface
    !> pot += Re V at the xc points for a potential V of
    !> SysInteraction_FillInteractionPotential.
    subroutine I_SysXc_AddInteractionPotential(pot, interactionPotential)
      import :: R64
      real(R64), intent(inout), contiguous, target :: pot(:)
      complex(R64), intent(in), contiguous, target :: interactionPotential(:)
    end subroutine
  end interface

  !> @brief Gradient of a field on the xc points (gradient functionals).
  procedure(I_SysXc_ApplyGradient), pointer :: SysXc_ApplyGradient => null()
  abstract interface
    !> gradF(:, c) = (G f) along the orthonormal directions c = 1..SysXc_dimension.
    subroutine I_SysXc_ApplyGradient(gradF, f)
      import :: R64
      real(R64), intent(out), contiguous, target :: gradF(:, :)
      real(R64), intent(in), contiguous, target :: f(:)
    end subroutine
  end interface

  !> @brief Adjoint of the gradient in the quadrature metric.
  procedure(I_SysXc_AddGradientAdjoint), pointer :: SysXc_AddGradientAdjoint => null()
  abstract interface
    !> f += Σ_c W⁻¹ G_cᵀ W gradF(:, c) (minus the discrete divergence).
    subroutine I_SysXc_AddGradientAdjoint(f, gradF)
      import :: R64
      real(R64), intent(inout), contiguous, target :: f(:)
      real(R64), intent(in), contiguous, target :: gradF(:, :)
    end subroutine
  end interface

  !=============================================================================
  ! module procedures pointers: functional
  !=============================================================================

  !> @brief Energy density and its partial derivatives on point arrays.
  procedure(I_SysXc_EvaluateFunctional), pointer :: SysXc_EvaluateFunctional
  abstract interface
    !> energyDensity = e(ρ, σ) (energy per volume), vRho = ∂e/∂ρ_s and
    !> vSigma = ∂e/∂σ_c for the nSpin = size(rho, 2) channels; σ_c are
    !> |∇ρ|² (nSpin 1) or ∇ρ↑·∇ρ↑, ∇ρ↑·∇ρ↓, ∇ρ↓·∇ρ↓ (nSpin 2); sigma and
    !> vSigma have zero columns for semilocal functionals without gradients.
    subroutine I_SysXc_EvaluateFunctional(energyDensity, vRho, vSigma, rho, sigma)
      import :: R64
      real(R64), intent(out), contiguous :: energyDensity(:)
      real(R64), intent(out), contiguous :: vRho(:, :)
      real(R64), intent(out), contiguous :: vSigma(:, :)
      real(R64), intent(in), contiguous :: rho(:, :)
      real(R64), intent(in), contiguous :: sigma(:, :)
    end subroutine
  end interface

end module
