! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Linear family of SysXc: the nodes of the 1D grid.
submodule(M_SysXc_Linear) S_SysXc_Linear

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysXc_Linear_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Grid
    use M_SysXc
    use M_SysXc_Libxc

    call Say_Fabricate("sysXc.linear")

    if (.not. Json_GetExistence("grid.linear")) error stop "grid.linear is required for sysXc.linear"
    if (SysXc_hartreeQ .and. .not. Json_GetExistence("sysInteraction.linear")) &
      error stop "sysXc.linear: the Hartree potential requires sysInteraction.linear"

    !------------------------------------
    ! set values and procedure pointers (nodal defaults of S_SysXc otherwise)
    !------------------------------------

    SysXc_dimension = 1
    SysXc_nPoints = Grid_nPoints
    SysXc_Setup => Setup

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("sysXc.linear.libxc")) then
      call SysXc_Libxc_Fabricate("sysXc.linear.libxc")

    else
      error stop "sysXc.linear is missing one of: libxc"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Weights of the xc points: the quadrature weights of the nodes.
  subroutine Setup
    use M_Utils_Say
    use M_Grid
    use M_SysXc

    call Say_Setup("sysXc.linear")

    if (.not. allocated(Grid_weights)) error stop "sysXc.linear: call Grid_Setup before SysXc_Setup"
    SysXc_nPoints = Grid_nPoints
    SysXc_weights = Grid_weights

  end subroutine

end submodule
