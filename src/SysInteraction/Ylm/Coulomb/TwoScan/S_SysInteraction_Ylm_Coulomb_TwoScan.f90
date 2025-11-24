! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysInteraction_Ylm_Coulomb_TwoScan) S_SysInteraction_Ylm_Coulomb_TwoScan

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysInteraction_Ylm_Coulomb_TwoScan_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysInteraction
    use M_SysInteraction_Ylm

    call Say_Fabricate("sysInteraction.ylm.coulomb.twoScan")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysInteraction_Ylm_FillInteractionPotentialRadial => FillInteractionPotentialRadial

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillInteractionPotentialRadial(potLm, srcLm, l, m, time, bt1_, bt2_)
    use M_Utils_Constants, only: PI
    use M_Utils_UnusedVariables
    use M_Grid_Ylm
    use M_SysInteraction_Ylm_Coulomb, only: SysInteraction_Ylm_Coulomb_Strength

    complex(R64), intent(out), contiguous :: potLm(:)
    complex(R64), intent(in), contiguous  :: srcLm(:)
    integer(I32), intent(in)  :: l
    integer(I32), intent(in)  :: m
    real(R64), intent(in)     :: time
    integer(I32), intent(in), optional  :: bt1_
    integer(I32), intent(in), optional  :: bt2_

    integer(I32) :: i, nRad
    real(R64)    :: r_i, factor
    complex(R64), allocatable :: C1(:), C2(:)
    real(R64)    :: strength

    if (.false.) call UnusedVariables_Mark(m, time, bt1_, bt2_)

    nRad = Grid_Ylm_nRadial
    strength = SysInteraction_Ylm_Coulomb_Strength

    allocate (C1(nRad), C2(nRad))

    ! prefix: C1(i) = sum_{j<=i} rho(j) * w_j * r_j^l
    C1(1) = srcLm(1) * Grid_Ylm_radialPoints(1)**l
    do i = 2, nRad
      C1(i) = C1(i - 1) + srcLm(i) * Grid_Ylm_radialPoints(i)**l
    end do

    ! suffix: C2(i) = sum_{j>=i} rho(j) * w_j / r_j^(l+1)
    C2(nRad) = srcLm(nRad) / Grid_Ylm_radialPoints(nRad)**(l + 1)
    do i = nRad - 1, 1, -1
      C2(i) = C2(i + 1) + srcLm(i) / Grid_Ylm_radialPoints(i)**(l + 1)
    end do

    factor = strength * (4.0_R64 * PI / (2.0_R64 * real(l, R64) + 1.0_R64))

    do i = 1, nRad
      r_i = Grid_Ylm_radialPoints(i)
      ! exact discrete form (matches Std):
      ! sum_j rho_j w_j * min(r_i,r_j)^l / max(r_i,r_j)^(l+1)
      ! implemented as two scans without double-counting j=i
      potLm(i) = factor * (C1(i) / r_i**(l + 1) + (r_i**l) * C2(i) - srcLm(i) / r_i)
    end do

    deallocate (C1, C2)
  end subroutine

end submodule
