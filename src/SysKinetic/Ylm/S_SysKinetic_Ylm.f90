! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @file S_SysKinetic_Ylm.f90
!> @brief Implementation submodule for Ylm kinetic fabrication and operator.
submodule(M_SysKinetic_Ylm) S_SysKinetic_Ylm

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Bind the global kinetic operator and dispatch to radial backend.
  module subroutine SysKinetic_Ylm_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysKinetic
    use M_SysKinetic_Ylm_Laplacian

    call Say_Fabricate("sysKinetic.ylm")

    !------------------------------------
    ! bind the global kinetic operator
    !------------------------------------

    SysKinetic_MultiplyWithKineticOp => MultiplyWithKineticOp

    !------------------------------------
    ! dispatch to radial backend
    !------------------------------------

    if (Json_GetExistence("sysKinetic.ylm.laplacian")) then
      call SysKinetic_Ylm_Laplacian_Fabricate

    else
      error stop "sysKinetic.ylm is missing one of: laplacian"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @brief Apply T̂ψ by looping over all (l,m) channels.
  !>
  !> @details
  !> The full orbital ψ(r,θ,φ) is stored as a 1D array with radial points for
  !> each (l,m) channel concatenated. This routine:
  !>   1. Extracts each channel f_{lm}(r)
  !>   2. Applies the radial kinetic operator T̂_{lm}
  !>   3. Stores the result back in the output array
  !>
  !> Channels with negligible amplitude (|f| < 10⁻¹⁴) are skipped for efficiency.
  subroutine MultiplyWithKineticOp(dOrb, orb, time, bt_)
    use M_Grid
    use M_Grid_Ylm
    use M_Utils_Constants

    complex(R64), intent(out), contiguous :: dOrb(:)
    complex(R64), intent(in), contiguous :: orb(:)
    real(R64), intent(in)                :: time
    integer(I32), intent(in), optional   :: bt_

    integer(I32) :: l, m, lmax
    integer(I32) :: nRad
    complex(R64), allocatable :: orbLm(:), dOrbLm(:)

    nRad = Grid_Ylm_nRadial
    lmax = Grid_Ylm_lmax

    dOrb = 0.0_R64

    allocate (orbLm(nRad), dOrbLm(nRad))

    ! Process each (l,m) channel separately
    do l = 0, lmax
      do m = -l, l
        ! Extract radial points for this (l,m) channel
        call Grid_Ylm_GetLmComponent(orbLm, l, m, orb)

        ! Skip negligible channels for efficiency
        if (all(abs(orbLm) < 1.0e-14_R64)) cycle

        ! Apply radial kinetic operator
        call SysKinetic_Ylm_MultiplyWithRadialKineticOp(dOrbLm, orbLm, l, m, time, bt_)

        ! Store result back in output array
        call Grid_Ylm_SetLmComponent(dOrb, l, m, dOrbLm)
      end do
    end do

    deallocate (orbLm, dOrbLm)

  end subroutine

end submodule
