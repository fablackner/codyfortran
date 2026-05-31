! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Default numerical implementation for linear harmonic potential.
!>
!> Evaluates a multi-dimensional separable harmonic trap on a 1D grid:
!>   V(x) = Σₙ ½ωₙ²(x - xₙ)²
!>
!> Each term represents a harmonic oscillator centered at `position[n]` with
!> angular frequency `omega[n]`. Multiple terms can be combined to create
!> double-well or more complex trapping geometries.
submodule(M_SysPotential_Linear_Harmonic_StdImpl) S_SysPotential_Linear_Harmonic_StdImpl

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysPotential_Linear_Harmonic_StdImpl_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysPotential
    use M_SysPotential_Linear

    call Say_Fabricate("sysPotential.linear.harmonic.stdImpl")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysPotential_FillExternalPotential => FillExternalPotential

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Fill the external potential array with harmonic trap values.
  !>
  !> Sums contributions from all harmonic oscillator centers defined in the
  !> configuration. Uses vectorized array operations for efficiency.
  subroutine FillExternalPotential(externalPotential, time, bt_)
    use M_Utils_UnusedVariables
    use M_SysPotential_Linear_Harmonic
    use M_Grid
    use M_Grid_Linear

    complex(R64), intent(out), allocatable :: externalPotential(:)
    real(R64), intent(in) :: time
    integer(I32), intent(in), optional :: bt_

    integer(I32) :: i
    real(R64), allocatable :: position(:), omega(:)

    if (.false.) call UnusedVariables_Mark(time, bt_)

    position = SysPotential_Linear_Harmonic_position
    omega = SysPotential_Linear_Harmonic_omega

    if (.not. allocated(externalPotential)) allocate (externalPotential(Grid_nPoints))
    externalPotential(:) = 0.0_R64

    do i = 1, size(position)
      externalPotential(:) = externalPotential(:) + 0.5_R64 * omega(i)**2 * (Grid_Linear_xCoord(:) - position(i))**2
    end do

  end subroutine

end submodule
