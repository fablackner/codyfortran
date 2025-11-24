! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysPotential_Ylm_Coulomb_StdImpl) S_SysPotential_Ylm_Coulomb_StdImpl

  implicit none

contains

  module subroutine SysPotential_Ylm_Coulomb_StdImpl_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysPotential
    use M_SysPotential_Ylm

    call Say_Fabricate("sysPotential.ylm.coulomb.stdImpl")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysPotential_Ylm_FillExternalPotentialRadial => FillExternalPotentialRadial

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillExternalPotentialRadial(potLm, l, m, time, bt_)
    use M_Utils_UnusedVariables
    use M_Utils_Constants
    use M_SysPotential_Ylm_Coulomb
    use M_Grid
    use M_Grid_Ylm

    complex(R64), intent(out), allocatable :: potLm(:)
    integer(I32), intent(in)  :: l
    integer(I32), intent(in)  :: m
    real(R64), intent(in)  :: time
    integer(I32), intent(in), optional  :: bt_

    real(R64) :: r
    integer(I32) :: iRad
    real(R64) :: prefactor

    if (.false.) call UnusedVariables_Mark(time, bt_)
    if (l .ne. 0 .or. m .ne. 0) then
      error stop "SysPotential_Ylm_Coulomb_StdImpl only has (l=0,m=0) component"
    end if

    if (.not. allocated(potLm)) allocate (potLm(Grid_Ylm_nRadial))

    ! Get nuclear charge
    prefactor = SysPotential_Ylm_Coulomb_charge * sqrt(4.0_R64 * PI)

    ! Std Coulomb potential: V(r) = -Z/r in spherical harmonics
    ! Only the (l=0, m=0) component is non-zero
    do iRad = 1, Grid_Ylm_nRadial
      r = Grid_Ylm_radialPoints(iRad)
      potLm(iRad) = -prefactor / r
    end do

  end subroutine

end submodule
