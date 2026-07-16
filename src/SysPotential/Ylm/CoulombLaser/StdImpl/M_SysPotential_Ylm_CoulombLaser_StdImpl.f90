! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Default numerical implementation for the Ylm Coulomb-plus-laser potential.
!>
!> Provides the radial (l,m) components of the combined central Coulomb and
!> length-gauge dipole laser potential and a factory to register it.
module M_SysPotential_Ylm_CoulombLaser_StdImpl
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Register the default Ylm Coulomb-plus-laser implementation.
    module subroutine SysPotential_Ylm_CoulombLaser_StdImpl_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
