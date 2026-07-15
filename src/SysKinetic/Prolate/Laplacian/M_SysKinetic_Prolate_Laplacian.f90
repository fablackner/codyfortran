! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file M_SysKinetic_Prolate_Laplacian.f90
!> @brief Prolate-spheroidal Laplacian backends for kinetic operators.
!>
!> @details
!> Hosts shared data and fabrication for prolate kinetic implementations.
!>
!> JSON Configuration
!> ------------------
!> @code{.json}
!> {
!>   "sysKinetic": {
!>     "prolate": {
!>       "laplacian": {
!>         "bodyMass": [1.0],
!>         "fedvr": {}
!>       }
!>     }
!>   }
!> }
!> @endcode
!>
!> @see M_SysKinetic_Prolate_Laplacian_Fedvr
module M_SysKinetic_Prolate_Laplacian
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> @brief Fabricate the prolate Laplacian backend.
    module subroutine SysKinetic_Prolate_Laplacian_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data (shared across backends)
  !=============================================================================

  !> @brief Body masses used for kinetic scaling via 1/(2 m_bt).
  !>
  !> Default is [1.0] (electron mass in atomic units).
  real(R64), allocatable :: SysKinetic_Prolate_Laplacian_bodyMass(:)

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
