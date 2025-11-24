! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysKinetic_Linear_Laplacian_FinDiff) S_SysKinetic_Linear_Laplacian_FinDiff
  use M_Utils_DerivativeFinDiff, only: T_DerivativeFinDiff_Ctx

  implicit none

  ! Saved finite-difference context (initialized in Setup)
  type(T_DerivativeFinDiff_Ctx) :: dfCtx

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysKinetic_Linear_Laplacian_FinDiff_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysKinetic

    call Say_Fabricate("sysKinetic.linear.laplacian.finDiff")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysKinetic_MultiplyWithKineticOp => MultiplyWithKineticOp
    SysKinetic_Setup => Setup

    if (.not. Json_GetExistence("grid.linear.const")) then
      error stop "grid.linear.const is required for sysKinetic.linear.laplacian.finDiff"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Setup
    use M_Grid
    use M_Grid_Linear
    use M_Utils_DerivativeFinDiff
    real(R64) :: dx

    ! Compute spacing
    dx = (Grid_Linear_xmax - Grid_Linear_xmin) / (Grid_nPoints - 1)

    call DerivativeFinDiff_CreateCtx(dfCtx, Grid_nPoints, dx)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine MultiplyWithKineticOp(dOrb, orb, time, bt_)
    use M_Utils_DerivativeFinDiff
    use M_Utils_UnusedVariables
    use M_SysKinetic_Linear_Laplacian
    use M_Grid
    use M_Grid_Linear

    complex(R64), intent(out), contiguous    :: dOrb(:)
    complex(R64), intent(in), contiguous     :: orb(:)
    real(R64), intent(in)                    :: time
    integer(I32), intent(in), optional       :: bt_

    complex(R64), allocatable :: derivative2nd(:)
    integer(I32) :: nG
    integer(I32) :: bt

    if (.false.) call UnusedVariables_Mark(time, bt_)

    if (.not. present(bt_)) bt = 1
    if (present(bt_)) bt = bt_

    nG = Grid_nPoints

    allocate (derivative2nd(nG))

    ! Use pre-created context
    call DerivativeFinDiff_Do2ndDerivativeZeroEnds(derivative2nd, orb, dfCtx)

    dOrb = -0.5_R64 * derivative2nd(:) / SysKinetic_Linear_Laplacian_bodyMass(bt)

    deallocate (derivative2nd)
  end subroutine

end submodule
