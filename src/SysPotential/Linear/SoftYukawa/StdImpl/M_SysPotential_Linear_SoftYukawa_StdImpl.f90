! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Default numerical implementation for linear soft-Yukawa potential.
!>
!> Provides a reference implementation for a softened Yukawa external potential
!> on linear grids and a factory to register it. Alternative backends may
!> substitute their own factory for performance or accuracy trade-offs.
module M_SysPotential_Linear_SoftYukawa_StdImpl
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Register the default linear soft-Yukawa implementation and wire pointers.
    module subroutine SysPotential_Linear_SoftYukawa_StdImpl_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
