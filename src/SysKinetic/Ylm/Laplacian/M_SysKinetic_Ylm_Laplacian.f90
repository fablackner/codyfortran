! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file M_SysKinetic_Ylm_Laplacian.f90
!> @brief Ylm radial Laplacian backends for kinetic operators.
!>
!> @details
!> Hosts shared data and fabrication for Ylm radial-kinetic implementations.
!> The operator acts on individual (l,m) channels via:
!>
!>     T̂_{lm} f(r) = −(1/2m) [ d²f/dr² + (2/r)df/dr − l(l+1)/r² f ]
!>
!> Available backends:
!>   - **FinDiff**: Finite-difference on uniform radial grid (requires `grid.ylm.const`)
!>   - **FEDVR**: Finite-element DVR for nonuniform grids (requires `grid.ylm.fedvr`)
!>
!> JSON Configuration
!> ------------------
!> @code{.json}
!> {
!>   "sysKinetic": {
!>     "ylm": {
!>       "laplacian": {
!>         "bodyMass": [1.0],
!>         "finDiff": {}
!>       }
!>     }
!>   }
!> }
!> @endcode
!>
!> @see M_SysKinetic_Ylm_Laplacian_FinDiff, M_SysKinetic_Ylm_Laplacian_Fedvr
module M_SysKinetic_Ylm_Laplacian
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> @brief Fabricate the Ylm radial Laplacian backend.
    !>
    !> Selects the numerical scheme (FinDiff/FEDVR) and prepares masses.
    module subroutine SysKinetic_Ylm_Laplacian_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data (shared across FinDiff/FEDVR backends)
  !=============================================================================

  !> @brief Body masses used for radial kinetic scaling via 1/(2 m_bt).
  !>
  !> Default is [1.0] (electron mass in atomic units).
  real(R64), allocatable :: SysKinetic_Ylm_Laplacian_bodyMass(:)

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
