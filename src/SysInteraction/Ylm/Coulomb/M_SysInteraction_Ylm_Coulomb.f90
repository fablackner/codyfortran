! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Coulomb interaction configuration for Ylm back-ends.
!>
!> Holds common parameters for Coulomb-like interactions in spherical-harmonic
!> representations and provides the fabrication hook to select a concrete
!> numerical scheme (see submodules).
module M_SysInteraction_Ylm_Coulomb
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Register a Coulomb solver variant (BlockEq, FullEq, TwoScan, ...)
    !> and propagate the strength and other settings to the back-end.
    module subroutine SysInteraction_Ylm_Coulomb_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Overall Coulomb coupling constant. Sign and units follow the global
  !> Hamiltonian conventions used in the simulation.
  real(R64) :: SysInteraction_Ylm_Coulomb_Strength

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
