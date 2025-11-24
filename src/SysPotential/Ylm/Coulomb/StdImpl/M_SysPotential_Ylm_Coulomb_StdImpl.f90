! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Default numerical implementation for Ylm Coulomb potential.
!>
!> Provides a reference implementation of the Coulomb external potential in a
!> spherical-harmonics representation and a factory to register it. Alternative
!> backends may replace it via their own factory.
module M_SysPotential_Ylm_Coulomb_StdImpl
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Register the default Ylm Coulomb implementation and wire pointers.
    module subroutine SysPotential_Ylm_Coulomb_StdImpl_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
