! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Softened Yukawa-like external potential on linear grids.
!>
!> Models a sum of centers with charges `charge`, located at `position` (per
!> axis), using a softened Yukawa kernel with optional exponential damping.
!> Two softening radii (`softening1`, `softening2`) allow tuning near-field
!> behavior. The concrete procedures are wired by
!> `SysPotential_Linear_SoftYukawa_Fabricate`.
module M_SysPotential_Linear_SoftYukawa
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Parse parameters and wire the linear soft-Yukawa implementation.
    !>
    !> Reads configuration, initializes the per-center arrays below, and
    !> registers the implementation with the linear backend and grid-agnostic
    !> pointers.
    module subroutine SysPotential_Linear_SoftYukawa_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Center positions per axis (size: nCenters x nDims when expanded; flattened here per convention).
  real(R64), allocatable :: SysPotential_Linear_SoftYukawa_position(:)
  !> Charges per center.
  real(R64), allocatable :: SysPotential_Linear_SoftYukawa_charge(:)
  !> Primary softening radius (near-field regularization).
  real(R64), allocatable :: SysPotential_Linear_SoftYukawa_softening1(:)
  !> Secondary softening radius (additional smoothing control).
  real(R64), allocatable :: SysPotential_Linear_SoftYukawa_softening2(:)
  !> Exponential damping length-scale (Yukawa screening); use +Inf for none.
  real(R64), allocatable :: SysPotential_Linear_SoftYukawa_dampening(:)

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
