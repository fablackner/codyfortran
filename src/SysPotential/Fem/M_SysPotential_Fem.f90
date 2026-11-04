! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> FEM-grid (grid.fem) external potential backend.
!>
!> Potentials are nodal values on the interior nodes of the spectral-element
!> mesh, so the family binds the point-wise multiplication
!>   dOrb(i) = V(i) orb(i);
!> the model binds Setup and FillExternalPotential and owns its parameters
!> (coulomb: the nuclei of the grid.fem centers, M_SysPotential_Fem_Coulomb).
!>
!> Method refills the external potential on every right-hand-side evaluation
!> (no caching on the caller side), so a static model builds its potential
!> once in Setup and FillExternalPotential only copies it. Laser couplings are
!> not SysPotential terms on this grid: the length-gauge dipole coupling is the
!> SysGauge backend sysGauge.fem.length, which keeps SysPotential static (the
!> exact field-free part h0 = T + V can then be treated by
!> method.mb.orbBased.exactFieldFreeQ).
module M_SysPotential_Fem
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> FEM factory: checks grid.fem, binds the point-wise multiplication and
    !> dispatches to the potential model (coulomb).
    module subroutine SysPotential_Fem_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
