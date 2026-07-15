! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Standard xi solver for the prolate Coulomb interaction.
!>
!> @details
!> Per azimuthal channel m, the potential is expanded in normalized associated
!> Legendre functions of eta:
!>
!>     V_m(xi, eta) = sum_tau v_(tau m)(xi) Pbar_tau^m(eta)
!>
!> Each v_(tau m) solves the 1D boundary-value problem (a = R/2)
!>
!>     [ d/dxi (xi^2-1) d/dxi - tau(tau+1) - m^2/(xi^2-1) ] v = -(4 pi / a) g,
!>     g(xi_i) = sum_j Pbar_tau^m(eta_j) srcW_m(xi_i, eta_j),
!>
!> where srcW includes the metric weights. The weak-form DVR matrix
!> A = S_xi + diag(w_xi q) is LU-factorized once at setup for every (tau, m).
!> The homogeneous-Dirichlet solution is corrected by the exact exterior tail:
!> the difference between the true solution (which behaves like Q_tau^m(xi)
!> outside the source) and the Dirichlet solution is proportional to the
!> regular solution P_tau^m(xi), fixed by matching at ximax via the source
!> moment and the Wronskian constant kappa = (-1)^(m+1) (tau+m)!/(tau-m)!.
!>
!> For m /= 0 the xi = 1 grid point is excluded from the linear system
!> (channel Dirichlet point, consistent with the kinetic operator).
module M_SysInteraction_Prolate_Coulomb_StdImpl
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Fabricate the standard prolate Coulomb xi solver.
    module subroutine SysInteraction_Prolate_Coulomb_StdImpl_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !> LU factorization of the weak-form xi matrix for one (tau, |m|) pair.
  type :: T_SysInteraction_Prolate_Coulomb_StdImpl_Lu
    real(R64), allocatable :: lu(:, :)
    integer(I32), allocatable :: ipiv(:)
  end type

  !> LU solvers indexed (tau = 0:lmax, |m| = 0:mmaxAbs).
  type(T_SysInteraction_Prolate_Coulomb_StdImpl_Lu), allocatable :: SysInteraction_Prolate_Coulomb_StdImpl_solvers(:, :)

  !> Normalized on-cut Legendre values Pbar_tau^m(eta_j), indexed (tau, j, |m|).
  real(R64), allocatable :: SysInteraction_Prolate_Coulomb_StdImpl_pBarEta(:, :, :)

  !> Off-cut Legendre values P_tau^m(xi_i), indexed (tau, i, |m|).
  real(R64), allocatable :: SysInteraction_Prolate_Coulomb_StdImpl_pXi(:, :, :)

  !> Exterior-tail correction coefficients Q_tau^m(ximax)/(kappa P_tau^m(ximax)),
  !> indexed (tau, |m|).
  real(R64), allocatable :: SysInteraction_Prolate_Coulomb_StdImpl_corrCoef(:, :)

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
