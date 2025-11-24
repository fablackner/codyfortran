! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysKinetic_Ylm_Laplacian_FinDiff) S_SysKinetic_Ylm_Laplacian_FinDiff
  use M_Utils_DerivativeFinDiff, only: T_DerivativeFinDiff_Ctx

  implicit none

  type(T_DerivativeFinDiff_Ctx), save :: dfCtx

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysKinetic_Ylm_Laplacian_FinDiff_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysKinetic
    use M_SysKinetic_Ylm

    call Say_Fabricate("sysKinetic.ylm.laplacian.finDiff")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysKinetic_Ylm_MultiplyWithRadialKineticOp => MultiplyWithRadialKineticOp
    SysKinetic_Setup => Setup

    if (.not. Json_GetExistence("grid.ylm.const")) then
      error stop "grid.ylm.const is required for sysKinetic.ylm.laplacian.finDiff"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Setup
    use M_Grid_Ylm
    use M_Utils_DerivativeFinDiff

    integer(I32) :: nRad
    real(R64) :: dr

    nRad = Grid_Ylm_nRadial
    dr = Grid_Ylm_rmax / (Grid_Ylm_nRadial + 1)

    call DerivativeFinDiff_CreateCtx(dfCtx, nRad, dr)

  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> MultiplyWith the radial part of the Laplacian (∂²/∂r² + (2/r)∂/∂r)
  !> using the transformation g(r) = r*f(r), such that Laplacian(f) = (1/r) * d²g/dr².
  subroutine MultiplyWithRadialKineticOp(dOrbLm, orbLm, l, m, time, bt_)
    use M_Utils_UnusedVariables
    use M_Utils_DerivativeFinDiff
    use M_Grid_Ylm
    use M_SysKinetic_Ylm_Laplacian

    complex(R64), intent(out) :: dOrbLm(:) ! Output vector (result of D² * input)
    complex(R64), intent(in)  :: orbLm(:)   ! Input vector
    integer(I32), intent(in)  :: l
    integer(I32), intent(in)  :: m
    real(R64), intent(in)     :: time
    integer(I32), intent(in), optional :: bt_

    complex(R64), allocatable :: gRadial(:)  ! Temporary array for transformed function
    integer(I32) :: bt

    if (.false.) call UnusedVariables_Mark(m, time)

    if (.not. present(bt_)) bt = 1
    if (present(bt_)) bt = bt_

    ! Initialize output and temporary arrays
    allocate (gRadial(Grid_Ylm_nRadial))
    gRadial = 0.0_R64

    ! Transform: g(r) = r * f(r)
    gRadial(:) = Grid_Ylm_radialPoints(:) * orbLm(:)

    ! MultiplyWith second derivative to g(r) using FEDVR method
    call DerivativeFinDiff_Do2ndDerivative(dOrbLm, gRadial, dfCtx)

    ! Inverse transform: result = (1/r) * d²g/dr²
    dOrbLm(:) = dOrbLm(:) / Grid_Ylm_radialPoints(:)

    ! Add the angular part: -l(l+1)/r² term
    if (l > 0) then
      dOrbLm(:) = dOrbLm(:) - (l * (l + 1)) * orbLm(:) / Grid_Ylm_radialPoints(:)**2
    end if

    dOrbLm(:) = -0.5_R64 * dOrbLm(:) / SysKinetic_Ylm_Laplacian_bodyMass(bt)

    deallocate (gRadial)

  end subroutine

end submodule
