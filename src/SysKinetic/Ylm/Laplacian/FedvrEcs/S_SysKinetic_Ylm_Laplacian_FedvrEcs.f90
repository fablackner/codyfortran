! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file S_SysKinetic_Ylm_Laplacian_FedvrEcs.f90
!> @brief Implementation of the FEDVR-ECS radial kinetic operator for Ylm channels.
submodule(M_SysKinetic_Ylm_Laplacian_FedvrEcs) S_SysKinetic_Ylm_Laplacian_FedvrEcs

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Bind FEDVR-ECS radial operator, validate grid requirements.
  module subroutine SysKinetic_Ylm_Laplacian_FedvrEcs_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysKinetic
    use M_SysKinetic_Ylm

    call Say_Fabricate("sysKinetic.ylm.laplacian.fedvrEcs")

    !------------------------------------
    ! bind procedure pointers
    !------------------------------------

    SysKinetic_Ylm_MultiplyWithRadialKineticOp => MultiplyWithRadialKineticOp

    !------------------------------------
    ! validate grid requirements
    !------------------------------------

    if (.not. Json_GetExistence("grid.ylm.fedvrEcs")) then
      error stop "grid.ylm.fedvrEcs is required for sysKinetic.ylm.laplacian.fedvrEcs"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Apply the radial kinetic operator on the ECS contour using g(z) = z·f(z).
  !>
  !> @details
  !> Algorithm (same as Fedvr, but on the complex contour z(r)):
  !>   1. Transform: g(z) = z · f(z)
  !>   2. Compute d²g/dz² using the FEDVR-ECS derivative matrices
  !>   3. Inverse transform: result = (1/z) d²g/dz²
  !>   4. Add centrifugal term: − l(l+1)/z² · f(z)
  !>   5. Scale by −1/(2m)
  subroutine MultiplyWithRadialKineticOp(dOrbLm, orbLm, l, m, time, bt_)
    use M_Utils_FedvrEcs
    use M_Utils_UnusedVariables
    use M_Utils_DerivativeFedvrEcs
    use M_Grid
    use M_Grid_Ylm
    use M_Grid_Ylm_FedvrEcs, fedvrEcsCtx => Grid_Ylm_FedvrEcs_fedvrEcsCtx, &
      derivativeCtx => Grid_Ylm_FedvrEcs_derivativeCtx, &
      contourPoints => Grid_Ylm_FedvrEcs_contourPoints
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

    ! Transform: g(z) = z · f(z) on the contour
    gRadial(:) = contourPoints(:) * orbLm(:)

    ! Compute d²g/dz² using the FEDVR-ECS method
    call DerivativeFedvrEcs_Do2ndDerivative(dOrbLm, gRadial, derivativeCtx, fedvrEcsCtx)

    ! Inverse transform: result = (1/z) · d²g/dz²
    dOrbLm(:) = dOrbLm(:) / contourPoints(:)

    ! Add centrifugal term: −l(l+1)/z² · f(z)
    if (l > 0) then
      dOrbLm(:) = dOrbLm(:) - (l * (l + 1)) * orbLm(:) / contourPoints(:)**2
    end if

    ! Scale by −1/(2m)
    dOrbLm(:) = -0.5_R64 * dOrbLm(:) / SysKinetic_Ylm_Laplacian_bodyMass(bt)

    deallocate (gRadial)

  end subroutine

end submodule
