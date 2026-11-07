! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Fem family of SysXc: nodal xc points and the FEM nodal gradient.
submodule(M_SysXc_Fem) S_SysXc_Fem

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysXc_Fem_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Grid
    use M_SysXc
    use M_SysXc_Libxc

    call Say_Fabricate("sysXc.fem")

    if (.not. Json_GetExistence("grid.fem")) error stop "grid.fem is required for sysXc.fem"
    if (SysXc_hartreeQ .and. .not. Json_GetExistence("sysInteraction.fem")) &
      error stop "sysXc.fem: the Hartree potential requires sysInteraction.fem"

    !------------------------------------
    ! set values and procedure pointers (nodal defaults of S_SysXc otherwise)
    !------------------------------------

    SysXc_dimension = 3
    SysXc_nPoints = Grid_nPoints
    SysXc_Setup => Setup
    SysXc_ApplyGradient => ApplyGradient
    SysXc_AddGradientAdjoint => AddGradientAdjoint

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("sysXc.fem.libxc")) then
      call SysXc_Libxc_Fabricate("sysXc.fem.libxc")

    else
      error stop "sysXc.fem is missing one of: libxc"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Weights of the xc points: the lumped-mass weights of the nodes.
  subroutine Setup
    use M_Utils_Say
    use M_Grid
    use M_Grid_Fem
    use M_SysXc

    call Say_Setup("sysXc.fem")

    if (.not. allocated(Grid_Fem_weights)) error stop "sysXc.fem: call Grid_Setup before SysXc_Setup"
    SysXc_nPoints = Grid_nPoints
    SysXc_weights = Grid_Fem_weights

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ApplyGradient(gradF, f)
    use M_Grid_Fem

    real(R64), intent(out), contiguous, target :: gradF(:, :)
    real(R64), intent(in), contiguous, target :: f(:)

    call Grid_Fem_ApplyGradient(gradF, f)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine AddGradientAdjoint(f, gradF)
    use M_Grid_Fem

    real(R64), intent(inout), contiguous, target :: f(:)
    real(R64), intent(in), contiguous, target :: gradF(:, :)

    call Grid_Fem_AddGradientAdjoint(f, gradF)

  end subroutine

end submodule
