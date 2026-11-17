! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Eigen variant of the exact field-free part (Ylm grids): the blocks of
!> h0 = T + V_ext per l are diagonalized once and the orbitals are propagated
!> in their eigenbasis, exp(-i h0 dt) or exp(-h0 dt) in imaginary time.
!> Provides the phi functions for exponential integrators and the spectral
!> cut.
!>
!> JSON ("fieldFree.eigen"):
!>   - cut:  spectral cut energy of h0 (default 0 = off): the orbitals are
!>           restricted to the h0 eigenvectors below it, so the unphysical
!>           high-energy FEDVR states are never populated. For TDKS with a
!>           GGA it is the stability condition of the split step,
!>           cut * dt < pi.
!>   - gpuQ: cuBLAS variant of the block transformations (default false)
!>
!> Requires grid.ylm with a Hermitian metric (no ECS), a time- and
!> body-type-independent kinetic operator and a central static potential.
module M_FieldFree_Eigen
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Reads "fieldFree.eigen" and binds FieldFree_Setup, FieldFree_Propagate,
    !> FieldFree_ApplyPhi and, with a cut, FieldFree_ProjectBelowCut.
    module subroutine FieldFree_Eigen_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> Spectral cut (fieldFree.eigen.cut, 0 = off).
  real(R64) :: FieldFree_Eigen_cut = 0.0_R64

  !> Whether the block transformations run on the GPU (fieldFree.eigen.gpuQ).
  logical :: FieldFree_Eigen_gpuQ = .false.

end module
