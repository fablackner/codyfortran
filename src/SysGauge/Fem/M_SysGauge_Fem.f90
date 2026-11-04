! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Laser coupling on the 3D Cartesian spectral-element grid (grid.fem).
!>
!> Orbitals are nodal values, so a dipole coupling in length gauge is a
!> diagonal (point-wise) operator. It is nevertheless a SysGauge backend, not a
!> SysPotential model: SysPotential on grid.fem stays static (the nuclear
!> attraction), so the exact field-free part of Method
!> (method.mb.orbBased.exactFieldFreeQ) can treat h0 = T + V_static exactly and
!> the time-dependent laser term enters only through the gauge operator.
!>
!> @see M_SysGauge_Fem_Length
module M_SysGauge_Fem
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Fabricate the FEM gauge backend: checks grid.fem and dispatches to the
    !> gauge model (`length`).
    module subroutine SysGauge_Fem_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
