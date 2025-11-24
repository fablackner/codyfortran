! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
submodule(M_SysInteraction_Ylm_Coulomb_StdImpl) S_SysInteraction_Ylm_Coulomb_StdImpl

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysInteraction_Ylm_Coulomb_StdImpl_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysInteraction
    use M_SysInteraction_Ylm

    call Say_Fabricate("sysInteraction.ylm.coulomb.stdImpl")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysInteraction_Ylm_FillInteractionPotentialRadial => FillInteractionPotentialRadial

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillInteractionPotentialRadial(potLm, srcLm, l, m, time, bt1_, bt2_)
    use M_Utils_Constants, only: PI
    use M_Utils_UnusedVariables
    use M_SysInteraction_Ylm_Coulomb
    use M_Grid_Ylm

    complex(R64), intent(out), contiguous :: potLm(:)
    complex(R64), intent(in), contiguous  :: srcLm(:)
    integer(I32), intent(in)  :: l
    integer(I32), intent(in)  :: m
    real(R64), intent(in)     :: time
    integer(I32), intent(in), optional  :: bt1_
    integer(I32), intent(in), optional  :: bt2_

    integer(I32) :: iRad1, iRad2, nRad
    real(R64) :: factor, strength, r1, r2, rMin, rMax

    if (.false.) call UnusedVariables_Mark(m, time, bt1_, bt2_)

    nRad = Grid_Ylm_nRadial
    strength = SysInteraction_Ylm_Coulomb_Strength

    potLm(:) = (0.0_R64, 0.0_R64)

    do iRad2 = 1, nRad
      r2 = Grid_Ylm_radialPoints(iRad2)

      do iRad1 = 1, nRad
        r1 = Grid_Ylm_radialPoints(iRad1)
        rMin = min(r1, r2)
        rMax = max(r1, r2)

        ! Coulomb interaction with Green's function kernel
        potLm(iRad1) = potLm(iRad1) + srcLm(iRad2) * (rMin**l) / (rMax**(l + 1))
      end do
    end do

    factor = (4.0_R64 * PI / (2.0_R64 * real(l, R64) + 1.0_R64))
    potLm(:) = factor * strength * potLm(:)

  end subroutine

end submodule
