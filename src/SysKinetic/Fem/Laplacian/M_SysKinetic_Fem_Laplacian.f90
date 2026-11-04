! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file M_SysKinetic_Fem_Laplacian.f90
!> @brief Spectral-element Laplacian kinetic operator T = -lap/(2m) on grid.fem.
!>
!> @details
!> Binds
!>   - SysKinetic_MultiplyWithKineticOp: dOrb = -lap orb / (2 m_bt), with the
!>     banded nodal Laplacian Grid_Fem_ApplyLaplacian (Dirichlet box);
!>   - SysKinetic_SolveShifted (optional hook): the exact solution of
!>     (T + shift) dOrb = orb,
!>       dOrb = 2m (-lap + 2m shift)^-1 orb,
!>     by fast diagonalization (Grid_Fem_SolveShiftedLaplacian) for a real or a
!>     complex shift.
!> m_bt is the mass of body type bt_ (bodyMass(1) if bt_ is absent). A single
!> bodyMass entry is the mass of every body type; several entries give one
!> mass per body type, and a body type beyond them stops. Both operators are
!> reentrant (no module-level scratch) and may be called inside OpenMP
!> parallel regions.
!>
!> JSON Configuration
!> ------------------
!> @code{.json}
!> {
!>   "sysKinetic": {
!>     "fem": {
!>       "laplacian": {
!>         "bodyMass": [1.0, 1.0]
!>       }
!>     }
!>   }
!> }
!> @endcode
!>
!> | Key        | Type    | Default | Description                                   |
!> |------------|---------|---------|-----------------------------------------------|
!> | `bodyMass` | real(:) | [1.0]   | One mass for all body types, or one per body  |
!> |            |         |         | type (every entry > 0)                        |
!>
!> Flags: SysKinetic_timeIndependentQ = .true.; SysKinetic_bodyTypeIndependentQ
!> = .true. if and only if all masses are equal (always for a single entry).
module M_SysKinetic_Fem_Laplacian
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> @brief Read the masses, set the flags and bind the kinetic operator and
    !> the exact shifted solve.
    module subroutine SysKinetic_Fem_Laplacian_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> @brief Mass per body type entering T = -lap/(2 m_bt) (every entry > 0).
  real(R64), allocatable :: SysKinetic_Fem_Laplacian_bodyMass(:)

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
