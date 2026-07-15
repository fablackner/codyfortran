! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file M_SysKinetic_Prolate_Laplacian_Fedvr.f90
!> @brief FEDVR kinetic operator for prolate azimuthal channels.
!>
!> @details
!> Applies the channel kinetic operator using the Sturm-Liouville DVR matrices
!> stored in the grid layer (Grid_Prolate_xiKinMatrix, Grid_Prolate_etaKinMatrix).
!> Requires `grid.prolate.fedvr`.
module M_SysKinetic_Prolate_Laplacian_Fedvr
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> @brief Fabricate the FEDVR prolate kinetic backend.
    module subroutine SysKinetic_Prolate_Laplacian_Fedvr_Fabricate
    end subroutine
  end interface

  !=============================================================================
  ! module data
  !=============================================================================

  !=============================================================================
  ! module procedures pointers
  !=============================================================================

end module
