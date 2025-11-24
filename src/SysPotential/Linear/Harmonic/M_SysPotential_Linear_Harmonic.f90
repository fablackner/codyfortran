! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Harmonic potential on linear grids.
!>
!> Implements a (potentially multi-dimensional) separable harmonic potential on
!> linear/structured grids. The potential is centered at `position` with per-
!> axis angular frequencies `omega`. The concrete procedures are wired by
!> `SysPotential_Linear_Harmonic_Fabricate`.
module M_SysPotential_Linear_Harmonic
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Parse parameters and wire the linear harmonic implementation.
    !>
    !> Reads configuration, initializes `position` and `omega`, and registers
    !> the implementation with the linear backend and the grid-agnostic pointers.
    module subroutine SysPotential_Linear_Harmonic_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Trap center position per axis (size matches grid dimensionality).
  real(R64), allocatable :: SysPotential_Linear_Harmonic_position(:)
  !> Harmonic angular frequency per axis (size matches grid dimensionality).
  real(R64), allocatable :: SysPotential_Linear_Harmonic_omega(:)

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
