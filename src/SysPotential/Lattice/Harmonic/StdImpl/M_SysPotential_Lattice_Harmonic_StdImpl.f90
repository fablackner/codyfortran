! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Default numerical implementation for lattice harmonic potential.
!>
!> Provides a reference implementation of the harmonic external potential on
!> lattice grids and a factory to register it. Alternative optimized backends
!> may override this by providing their own factory.
module M_SysPotential_Lattice_Harmonic_StdImpl
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Register the default lattice harmonic implementation and wire pointers.
    module subroutine SysPotential_Lattice_Harmonic_StdImpl_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
