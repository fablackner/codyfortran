! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Spherical-harmonic (Ylm) domain interaction interface.
!>
!> @details This module handles particle-particle interactions in a spherical-
!> harmonic basis, primarily for atomic and molecular systems with central
!> symmetry. The interaction potential is expanded in spherical harmonics:
!>    V(r,Ω) = Σₗₘ Vₗₘ(r) Yₗₘ(Ω)
!>
!> **Physical context:**
!> For Coulomb interactions between electrons in atoms, the 1/|r₁-r₂| kernel
!> is expanded using the Laplace expansion:
!>    1/|r₁-r₂| = Σₗ (4π/(2l+1)) (r<ˡ/r>ˡ⁺¹) Σₘ Yₗₘ*(Ω₁)Yₗₘ(Ω₂)
!>
!> This separates the angular and radial parts, allowing efficient computation
!> by solving 1D radial Poisson equations for each (l,m) channel.
!>
!> **Available implementations:**
!>   - `Coulomb`: Pure 1/r interaction with multiple radial solvers
!>
!> **Module data:**
!>   - `lmax`: Maximum angular momentum in the expansion
!>   - `mmax`: Maximum |m| of the expansion (default min(2 grid mmax, lmax))
!>   - `nChannels`: Number of (l,m) channels of the interaction layout
!>   - `nPoints`: Length of a source/potential array (nChannels × nRadial)
!>   - `mIndependentQ`: True if potential is m-independent (central field)
module M_SysInteraction_Ylm
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Bind Ylm-based interaction procedures and set Ylm-specific parameters
    !> (e.g., lmax, m-independence) from configuration.
    module subroutine SysInteraction_Ylm_Fabricate
    end subroutine

    !> Potential channels of one m, potM(:, iL) = V_lm with l = |m| + iL - 1
    !> (channel matrix of Grid_Ylm_NodesToChannels), from the source channels
    !> srcM (radial metric weights included), by the per-channel solve of
    !> SysInteraction_FillInteractionPotential (body-type-independent
    !> interaction). Parallel over l unless called inside a parallel region.
    module subroutine SysInteraction_Ylm_FillInteractionPotentialM(potM, srcM, m, time)
      !> Output: potential channels (nRadial x number of l).
      complex(R64), intent(out), contiguous :: potM(:, :)
      !> Input: source channels (nRadial x number of l), radial metric weights included.
      complex(R64), intent(in), contiguous  :: srcM(:, :)
      !> Magnetic quantum number m of all channels.
      integer(I32), intent(in)              :: m
      !> Physical time associated with the source/potential.
      real(R64), intent(in)                 :: time
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Maximum angular momentum L included in the expansion.
  integer(I32) :: SysInteraction_Ylm_lmax
  !> Maximum |m| included in the expansion (sources and potentials store
  !> |m| <= min(l, mmax)). Defaults to min(2*Grid_Ylm_mmax, lmax): a product
  !> of two orbitals carries at most |m| = 2*Grid_Ylm_mmax.
  integer(I32) :: SysInteraction_Ylm_mmax
  !> Number of (l,m) channels of the interaction layout, Grid_Ylm_NChannels(lmax, mmax).
  integer(I32) :: SysInteraction_Ylm_nChannels
  !> Length of a source/potential array: nChannels * Grid_Ylm_nRadial.
  integer(I32) :: SysInteraction_Ylm_nPoints
  !> If true, the interaction does not depend on the magnetic quantum number m
  !> (only |m| or l matters), enabling additional optimizations.
  logical      :: SysInteraction_Ylm_mIndependentQ

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> Pointer to a routine that solves the radial part of the Poisson equation
  !> for a given spherical-harmonic channel (l, m).
  procedure(I_SysInteraction_Ylm_FillInteractionPotentialRadial), pointer :: SysInteraction_Ylm_FillInteractionPotentialRadial
  abstract interface
    !> Compute the radial potential component V_lm(r) from the corresponding
    !> source term rho_lm(r) for a specific (l, m) channel.
    subroutine I_SysInteraction_Ylm_FillInteractionPotentialRadial(potLm, srcLm, l, m, time, bt1_, bt2_)
      import :: I32, R64
      !> Output: radial potential component V_lm(r).
      complex(R64), intent(out), contiguous :: potLm(:)
      !> Input: radial density/source component rho_lm(r).
      complex(R64), intent(in), contiguous :: srcLm(:)
      !> Angular momentum quantum number l.
      integer(I32), intent(in) :: l
      !> Magnetic quantum number m.
      integer(I32), intent(in) :: m
      !> Physical time associated with the source/potential.
      real(R64), intent(in) :: time
      !> Optional body type (species) index for the potential target.
      integer(I32), intent(in), optional  :: bt1_
      !> Optional body type (species) index for the source.
      integer(I32), intent(in), optional  :: bt2_
    end subroutine
  end interface

end module
