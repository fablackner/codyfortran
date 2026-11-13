! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Laser coupling on the prolate spheroidal grid (grid.prolate, nuclei at the
!> foci z = ±R/2).
!>
!> The length-gauge coupling along the molecular axis, E(t) z with
!> z = (R/2) ξ η (Grid_Prolate_FillNodeCoordinate), is diagonal in the grid
!> points and keeps every azimuthal channel m; it is a SysGauge backend so
!> that sysPotential.prolate stays static (as on grid.fem).
!>
!> @see M_SysGauge_Length (sysGauge.prolate.length, polarization along z)
module M_SysGauge_Prolate
  use M_Utils_Types

  implicit none

  !=============================================================================
  ! module procedures
  !=============================================================================

  interface
    !> Fabricate the prolate gauge backend: checks grid.prolate and dispatches
    !> to the gauge model (`length`).
    module subroutine SysGauge_Prolate_Fabricate
    end subroutine
  end interface

end module
