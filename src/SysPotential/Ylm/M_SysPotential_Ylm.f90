! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Ylm-grid external potential backend (spherical harmonics basis).
!>
!> Provides the interface and state for potentials represented in a spherical
!> harmonics basis up to degree `lmax`. Concrete implementations (e.g.,
!> Coulomb) are selected and wired by `SysPotential_Ylm_Fabricate`. Some models
!> support `m`-independent forms which is indicated via `SysPotential_Ylm_mIndependentQ`.
module M_SysPotential_Ylm
  use M_Utils_Types
  use M_Utils_NoOpProcedures

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Ylm factory: parse config and wire concrete implementations.
    module subroutine SysPotential_Ylm_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Maximum spherical-harmonics degree used to represent the potential.
  integer(I32) :: SysPotential_Ylm_lmax
  !> True if the potential is independent of the magnetic quantum number m.
  logical      :: SysPotential_Ylm_mIndependentQ

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

  !> Ylm-specific setup hook (precompute radial grids, masks, caches, etc.).
  procedure(I_SysPotential_Ylm_Setup), pointer :: SysPotential_Setup => NoOpProcedures_Setup
  abstract interface
    !> Initialize the Ylm potential backend.
    subroutine I_SysPotential_Ylm_Setup
    end subroutine
  end interface

  !> Compute the radial component V_lm(r; t) of the external potential.
  procedure(I_SysPotential_Ylm_FillExternalPotentialRadial), pointer :: SysPotential_Ylm_FillExternalPotentialRadial
  abstract interface
    subroutine I_SysPotential_Ylm_FillExternalPotentialRadial(potLm, l, m, time, bt_)
      import :: I32, R64
      !> Output: allocated complex radial array for the requested (l, m) component.
      complex(R64), intent(out), allocatable :: potLm(:)
      !> Angular momentum quantum number l (0 <= l <= lmax).
      integer(I32), intent(in)              :: l
      !> Magnetic quantum number m (-l <= m <= l); ignored if m-independent.
      integer(I32), intent(in)              :: m
      !> Current simulation time.
      real(R64), intent(in)                 :: time
      !> Optional body type selector.
      integer(I32), intent(in), optional    :: bt_
    end subroutine
  end interface

end module
