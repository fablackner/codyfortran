! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file M_SysKinetic_Fem.f90
!> @brief Kinetic operators on the 3D Cartesian spectral-element grid (grid.fem).
!>
!> @details
!> Orbitals are nodal values on the interior Gauss-Lobatto-Legendre nodes of
!> the tensor-product mesh (see M_Grid_Fem). The kinetic operator is the
!> nodal strong form of the Galerkin Laplacian with lumped (diagonal) mass,
!>
!>     T f = -(1/2m) lap f,   lap = -W^-1 K3,
!>
!> self-adjoint in the metric of Grid_InnerProduct. The separable mesh admits
!> a fast exact inverse of T + shift (fast diagonalization), which the family
!> exposes through the optional hook SysKinetic_SolveShifted.
!>
!> @see M_SysKinetic_Fem_Laplacian, M_Grid_Fem
module M_SysKinetic_Fem
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> @brief Fabricate the FEM kinetic backend: checks grid.fem and dispatches
    !> to the operator model (`laplacian`).
    module subroutine SysKinetic_Fem_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
