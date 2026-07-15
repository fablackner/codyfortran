! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Prolate-spheroidal (two-center) domain interaction interface.
!>
!> @details Handles particle-particle interactions on the prolate grid, for
!> homonuclear diatomic molecules. The mean-field/pair potential is stored in
!> azimuthal channels like all other prolate fields:
!>    V(xi, eta, phi) = sum_m V_m(xi, eta) e^(i m phi)/sqrt(2 pi)
!>
!> **Physical context:**
!> For the Coulomb kernel 1/|r1 - r2| the phi integration decouples the m
!> channels; per channel the potential satisfies a 2D generalized Poisson
!> equation in (xi, eta) that separates further by expanding the eta
!> dependence in normalized associated Legendre functions P_l^m(eta) (the
!> Neumann expansion of the two-center Green's function). Each (l, m) term
!> requires only a 1D boundary-value solve in xi.
!>
!> **Available implementations:**
!>   - `Coulomb`: Pure 1/r interaction (stdImpl xi solver)
!>
!> **Module data:**
!>   - `mmax`: Maximum azimuthal quantum number of the potential channels
!>     (defaults to 2x the grid mmax to capture orbital-pair densities)
module M_SysInteraction_Prolate
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Bind prolate interaction procedures and set parameters from configuration.
    module subroutine SysInteraction_Prolate_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Maximum azimuthal quantum number of the interaction-potential channels.
  integer(I32) :: SysInteraction_Prolate_mmax

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> Pointer to a routine that solves one azimuthal channel of the potential.
  procedure(I_SysInteraction_Prolate_FillInteractionPotentialChannel), pointer :: &
    SysInteraction_Prolate_FillInteractionPotentialChannel
  abstract interface
    !> Compute the potential channel V_m(xi, eta) from the corresponding
    !> source channel rho_m(xi, eta) (which includes the metric weights).
    subroutine I_SysInteraction_Prolate_FillInteractionPotentialChannel(potM, srcM, m, time, bt1_, bt2_)
      import :: I32, R64
      !> Output: potential channel of size nXi*nEta (xi fastest).
      complex(R64), intent(out), contiguous :: potM(:)
      !> Input: weighted source channel of size nXi*nEta.
      complex(R64), intent(in), contiguous :: srcM(:)
      !> Azimuthal quantum number of the channel.
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
