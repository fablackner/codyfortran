! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Linear-grid Laplacian backends for kinetic operators.
!>
!> Provides shared data and the fabrication hook for linear-grid Laplacian
!> implementations (finite-difference, Fourier/spectral). The Laplacian is
!> typically used to form the kinetic energy operator −(1/2m)∇² acting on a
!> single orbital.
module M_SysKinetic_Linear_Laplacian
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Fabricate and configure the linear-grid Laplacian backend.
    !>
    !> Selects the numerical scheme (FD/Fourier), reads grid and mass settings,
    !> and prepares storage such as `SysKinetic_Linear_Laplacian_bodyMass`.
    module subroutine SysKinetic_Linear_Laplacian_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================
  !> Body masses used for scaling the kinetic operator via 1/(2 m_bt).
  !>
  !> The optional body-type index `bt_` supplied to operator routines selects
  !> the corresponding entry. If not present, the first (or only) mass is used.
  real(R64), allocatable :: SysKinetic_Linear_Laplacian_bodyMass(:)

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
