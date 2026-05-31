! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file S_SysKinetic_Linear_Laplacian_Fourier.f90
!> @brief Implementation of FFT-based spectral kinetic operator on linear grid.
submodule(M_SysKinetic_Linear_Laplacian_Fourier) S_SysKinetic_Linear_Laplacian_Fourier
  use M_Utils_DerivativeFftw, only: T_DerivativeFftw_Ctx

  implicit none

  !> Saved FFT derivative context (initialized in Setup, reused in operator)
  type(T_DerivativeFftw_Ctx), save :: dfCtx

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Bind Fourier operator and setup, validate grid requirements.
  module subroutine SysKinetic_Linear_Laplacian_Fourier_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysKinetic

    call Say_Fabricate("sysKinetic.linear.laplacian.fourier")

    !------------------------------------
    ! bind procedure pointers
    !------------------------------------

    SysKinetic_MultiplyWithKineticOp => MultiplyWithKineticOp
    SysKinetic_Setup => Setup

    !------------------------------------
    ! validate grid requirements
    !------------------------------------

    if (.not. Json_GetExistence("grid.linear.const")) then
      error stop "grid.linear.const is required for sysKinetic.linear.laplacian.fourier"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Initialize the FFT derivative context from grid parameters.
  subroutine Setup
    use M_Grid
    use M_Grid_Linear
    use M_Utils_DerivativeFftw
    real(R64) :: dx

    dx = (Grid_Linear_xmax - Grid_Linear_xmin) / (Grid_nPoints - 1)

    call DerivativeFftw_CreateCtx(dfCtx, Grid_nPoints, dx)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Apply T̂ψ = −(1/2m) d²ψ/dx² using FFT spectral method.
  !>
  !> @param[out] dOrb  Result: T̂ · orb
  !> @param[in]  orb   Input orbital wavefunction
  !> @param[in]  time  Simulation time (unused)
  !> @param[in]  bt_   Body type index for mass selection (default: 1)
  subroutine MultiplyWithKineticOp(dOrb, orb, time, bt_)
    use M_Utils_DerivativeFftw
    use M_Utils_UnusedVariables
    use M_SysKinetic_Linear_Laplacian
    use M_Grid
    use M_Grid_Linear

    complex(R64), intent(out), contiguous :: dOrb(:)
    complex(R64), intent(in), contiguous :: orb(:)
    real(R64), intent(in)               :: time
    integer(I32), intent(in), optional  :: bt_
    complex(R64), allocatable :: derivative2nd(:)
    integer(I32) :: nG
    integer(I32) :: bt

    if (.false.) call UnusedVariables_Mark(time, bt_)

    bt = 1
    if (present(bt_)) bt = bt_

    nG = Grid_nPoints
    allocate (derivative2nd(nG))

    call DerivativeFftw_Do2ndDerivative(derivative2nd, orb, dfCtx)

    dOrb = -0.5_R64 * derivative2nd(:) / SysKinetic_Linear_Laplacian_bodyMass(bt)

    deallocate (derivative2nd)
  end subroutine

end submodule
