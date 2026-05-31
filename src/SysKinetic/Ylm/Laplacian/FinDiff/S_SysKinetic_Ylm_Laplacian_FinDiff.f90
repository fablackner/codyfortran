! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file S_SysKinetic_Ylm_Laplacian_FinDiff.f90
!> @brief Implementation of FinDiff radial kinetic operator for Ylm channels.
submodule(M_SysKinetic_Ylm_Laplacian_FinDiff) S_SysKinetic_Ylm_Laplacian_FinDiff
  use M_Utils_DerivativeFinDiff, only: T_DerivativeFinDiff_Ctx

  implicit none

  !> Saved finite-difference context (initialized in Setup, reused in operator)
  type(T_DerivativeFinDiff_Ctx), save :: dfCtx

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Bind FinDiff radial operator and setup, validate grid requirements.
  module subroutine SysKinetic_Ylm_Laplacian_FinDiff_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysKinetic
    use M_SysKinetic_Ylm

    call Say_Fabricate("sysKinetic.ylm.laplacian.finDiff")

    !------------------------------------
    ! bind procedure pointers
    !------------------------------------

    SysKinetic_Ylm_MultiplyWithRadialKineticOp => MultiplyWithRadialKineticOp
    SysKinetic_Setup => Setup

    !------------------------------------
    ! validate grid requirements
    !------------------------------------

    if (.not. Json_GetExistence("grid.ylm.const")) then
      error stop "grid.ylm.const is required for sysKinetic.ylm.laplacian.finDiff"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Initialize the FD derivative context from radial grid parameters.
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
  !> @brief Apply the radial kinetic operator using the g(r) = r·f(r) transformation.
  !>
  !> @details
  !> Algorithm:
  !>   1. Transform: g(r) = r · f(r)
  !>   2. Compute d²g/dr² using finite differences
  !>   3. Inverse transform: result = (1/r) d²g/dr²
  !>   4. Add centrifugal term: − l(l+1)/r² · f(r)
  !>   5. Scale by −1/(2m)
  subroutine MultiplyWithRadialKineticOp(dOrbLm, orbLm, l, m, time, bt_)
    use M_Utils_UnusedVariables
    use M_Utils_DerivativeFinDiff
    use M_Grid_Ylm
    use M_SysKinetic_Ylm_Laplacian

    complex(R64), intent(out) :: dOrbLm(:)
    complex(R64), intent(in)  :: orbLm(:)
    integer(I32), intent(in)  :: l
    integer(I32), intent(in)  :: m
    real(R64), intent(in)     :: time
    integer(I32), intent(in), optional :: bt_

    complex(R64), allocatable :: gRadial(:)
    integer(I32) :: bt

    if (.false.) call UnusedVariables_Mark(m, time)

    bt = 1
    if (present(bt_)) bt = bt_

    allocate (gRadial(Grid_Ylm_nRadial))
    gRadial = 0.0_R64

    ! Transform: g(r) = r · f(r)
    gRadial(:) = Grid_Ylm_radialPoints(:) * orbLm(:)

    ! Compute d²g/dr² using finite differences
    call DerivativeFinDiff_Do2ndDerivative(dOrbLm, gRadial, dfCtx)

    ! Inverse transform: result = (1/r) · d²g/dr²
    dOrbLm(:) = dOrbLm(:) / Grid_Ylm_radialPoints(:)

    ! Add centrifugal term: −l(l+1)/r² · f(r)
    if (l > 0) then
      dOrbLm(:) = dOrbLm(:) - (l * (l + 1)) * orbLm(:) / Grid_Ylm_radialPoints(:)**2
    end if

    ! Scale by −1/(2m)
    dOrbLm(:) = -0.5_R64 * dOrbLm(:) / SysKinetic_Ylm_Laplacian_bodyMass(bt)

    deallocate (gRadial)

  end subroutine

end submodule
