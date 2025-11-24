! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysPotential_Linear_Harmonic_StdImpl) S_SysPotential_Linear_Harmonic_StdImpl

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysPotential_Linear_Harmonic_StdImpl_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysPotential
    use M_SysPotential_Linear

    call Say_Fabricate("cosinusLinear")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysPotential_FillExternalPotential => FillExternalPotential

  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillExternalPotential(externalPotential, time, bt_)
    ! array-filling implementation of harmonic potential
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

    ! Sum contributions from all harmonic oscillators
    do i = 1, size(position)
      externalPotential(:) = externalPotential(:) + 0.5_R64 * omega(i)**2 * (Grid_Linear_xCoord(:) - position(i))**2
    end do

  end subroutine

end submodule
