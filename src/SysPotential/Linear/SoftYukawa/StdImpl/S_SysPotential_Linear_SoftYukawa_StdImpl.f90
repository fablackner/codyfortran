! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysPotential_Linear_SoftYukawa_StdImpl) S_SysPotential_Linear_SoftYukawa_StdImpl

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysPotential_Linear_SoftYukawa_StdImpl_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysPotential
    use M_SysPotential_Linear

    call Say_Fabricate("sysPotential.linear.softYukawa.stdImpl")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysPotential_FillExternalPotential => FillExternalPotential

  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillExternalPotential(externalPotential, time, bt_)
    ! array-filling implementation of soft Yukawa potential
    use M_Utils_UnusedVariables
    use M_SysPotential_Linear_SoftYukawa
    use M_Grid
    use M_Grid_Linear

    complex(R64), intent(out), allocatable :: externalPotential(:)
    real(R64), intent(in) :: time
    integer(I32), intent(in), optional :: bt_

    integer(I32) :: i
    real(R64), allocatable :: position(:)
    real(R64), allocatable :: charge(:)
    real(R64), allocatable :: softening1(:)
    real(R64), allocatable :: softening2(:)
    real(R64), allocatable :: dampening(:)
    real(R64), allocatable :: distances(:)

    if (.false.) call UnusedVariables_Mark(time, bt_)

    position = SysPotential_Linear_SoftYukawa_position
    charge = SysPotential_Linear_SoftYukawa_charge
    softening1 = SysPotential_Linear_SoftYukawa_softening1
    softening2 = SysPotential_Linear_SoftYukawa_softening2
    dampening = SysPotential_Linear_SoftYukawa_dampening

    if (.not. allocated(externalPotential)) allocate (externalPotential(Grid_nPoints))
    allocate (distances(size(externalPotential)))

    ! Initialize result
    externalPotential(:) = 0.0_R64

    ! Sum contributions from all charges
    do i = 1, size(position)
      ! Calculate the potential contribution from each charge
      ! Using the soft Yukawa potential formula
      distances = abs(Grid_Linear_xCoord(:) - position(i))
      externalPotential(:) = externalPotential(:) - charge(i) * exp(-distances * dampening(i)) / &
                             (sqrt(distances**2 + softening1(i)**2) + softening2(i))
    end do

    deallocate (distances)

  end subroutine

end submodule
