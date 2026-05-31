! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file M_SysKinetic_Linear_Laplacian.f90
!> @brief Linear-grid Laplacian backends for kinetic operators.
!>
!> @details
!> Provides shared data and the fabrication hook for linear-grid Laplacian
!> implementations. The Laplacian ∇² forms the kinetic energy operator:
!>
!>     T̂ψ = −(1/2m) ∇²ψ
!>
!> Available backends:
!>   - **FinDiff**: Central finite-difference stencil (zero boundary conditions)
!>   - **Fourier**: FFT-based spectral method (periodic boundary conditions)
!>
!> JSON Configuration
!> ------------------
!> @code{.json}
!> {
!>   "sysKinetic": {
!>     "linear": {
!>       "laplacian": {
!>         "bodyMass": [1.0, 1836.0],
!>         "finDiff": {}
!>       }
!>     }
!>   }
!> }
!> @endcode
!>
!> @see M_SysKinetic_Linear_Laplacian_FinDiff, M_SysKinetic_Linear_Laplacian_Fourier
module M_SysKinetic_Linear_Laplacian
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> @brief Fabricate the linear-grid Laplacian backend.
    !>
    !> Selects the numerical scheme (FinDiff/Fourier), reads mass settings,
    !> and prepares `SysKinetic_Linear_Laplacian_bodyMass`.
    module subroutine SysKinetic_Linear_Laplacian_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data (shared across FinDiff/Fourier backends)
  !=============================================================================

  !> @brief Body masses used for scaling the kinetic operator via 1/(2 m_bt).
  !>
  !> The optional body-type index `bt_` supplied to operator routines selects
  !> the corresponding entry. If `bt_` is not present, uses bodyMass(1).
  !>
  !> Example: For electron (m=1) + proton (m=1836), set `bodyMass = [1.0, 1836.0]`.
  real(R64), allocatable :: SysKinetic_Linear_Laplacian_bodyMass(:)

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
