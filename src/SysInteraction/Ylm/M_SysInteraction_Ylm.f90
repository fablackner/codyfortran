! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Spherical-harmonic (Ylm) domain interaction registration.
!>
!> Fabricates interaction models that operate in a spherical-harmonic basis
!> with radial grids. It exposes a radial solver hook used by concrete Coulomb
!> implementations.
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
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Maximum angular momentum L included in the expansion.
  integer(I32) :: SysInteraction_Ylm_lmax
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
