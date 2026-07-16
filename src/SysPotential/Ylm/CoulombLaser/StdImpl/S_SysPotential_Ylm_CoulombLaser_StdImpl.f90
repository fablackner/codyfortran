! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Default numerical implementation for the Ylm Coulomb-plus-laser potential.
!>
!> Evaluates the (l,m) radial components of
!>   V(r, Omega, t) = -Z/r + E(t) * z
!> in a spherical-harmonics representation:
!>   V_00(r) = -Z * sqrt(4 pi) / r
!>   V_10(r) = E(t) * r * sqrt(4 pi / 3)
!> All other components (in particular l=1, m=+-1) vanish for a z-polarized
!> field. The sqrt(4 pi / 3) prefactor arises from z = r sqrt(4 pi/3) Y_10.
submodule(M_SysPotential_Ylm_CoulombLaser_StdImpl) S_SysPotential_Ylm_CoulombLaser_StdImpl

  implicit none

contains

  module subroutine SysPotential_Ylm_CoulombLaser_StdImpl_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysPotential
    use M_SysPotential_Ylm

    call Say_Fabricate("sysPotential.ylm.coulombLaser.stdImpl")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysPotential_Ylm_FillExternalPotentialRadial => FillExternalPotentialRadial

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Fill the radial component V_{lm}(r) of the Coulomb-plus-laser potential.
  !>
  !> Only (l=0, m=0) and (l=1, m=0) are nonzero; requests up to l=1 are valid
  !> since `SysPotential_Ylm_lmax = 1`.
  subroutine FillExternalPotentialRadial(potLm, l, m, time, bt_)
    use M_Utils_UnusedVariables
    use M_Utils_Constants
    use M_SysPotential_Ylm_CoulombLaser
    use M_Grid
    use M_Grid_Ylm

    complex(R64), intent(out), allocatable :: potLm(:)
    integer(I32), intent(in)  :: l
    integer(I32), intent(in)  :: m
    real(R64), intent(in)  :: time
    integer(I32), intent(in), optional  :: bt_

    real(R64) :: prefactor

    if (.false.) call UnusedVariables_Mark(bt_)
    if (l < 0 .or. l > 1 .or. abs(m) > l) then
      error stop "SysPotential_Ylm_CoulombLaser_StdImpl only has components up to l=1"
    end if

    if (.not. allocated(potLm)) allocate (potLm(Grid_Ylm_nRadial))

    if (l .eq. 0) then
      prefactor = SysPotential_Ylm_CoulombLaser_charge * sqrt(4.0_R64 * PI)
      potLm(:) = -prefactor / Grid_Ylm_radialPoints(:)

    else if (m .eq. 0) then
      prefactor = SysPotential_Ylm_CoulombLaser_FieldAmplitude(time) * sqrt(4.0_R64 * PI / 3.0_R64)
      potLm(:) = prefactor * Grid_Ylm_radialPoints(:)

    else
      potLm(:) = 0.0_R64
    end if

  end subroutine

end submodule
