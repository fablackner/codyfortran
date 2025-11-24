! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Default numerical implementation for lattice seesaw potential.
!>
!> Provides a reference implementation of the seesaw external potential on
!> lattice grids and a factory to register it. Alternative optimized backends
!> may override this by providing their own factory.
module M_SysPotential_Lattice_Seesaw_StdImpl
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Register the default lattice seesaw implementation and wire pointers.
    module subroutine SysPotential_Lattice_Seesaw_StdImpl_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
