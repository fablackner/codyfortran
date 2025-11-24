! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysInteraction_Ylm) S_SysInteraction_Ylm

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysInteraction_Ylm_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Grid_Ylm
    use M_SysInteraction
    use M_SysInteraction_Ylm_Coulomb

    call Say_Fabricate("sysInteraction.ylm")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysInteraction_Ylm_lmax = Json_Get("sysInteraction.ylm.lmax", 2 * Grid_Ylm_lmax)

    SysInteraction_FillInteractionPotential => FillInteractionPotential
    SysInteraction_FillInteractionSrc => FillInteractionSrc
    SysInteraction_MultiplyWithInteractionPotential => MultiplyWithInteractionPotential

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("sysInteraction.ylm.coulomb")) then
      call SysInteraction_Ylm_Coulomb_Fabricate

    else
      error stop "sysInteraction.ylm is missing one of: coulomb"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillInteractionPotential(interactionPotential, src, time, bt1_, bt2_)
    use M_Grid
    use M_Grid_Ylm

    complex(R64), intent(out), allocatable :: interactionPotential(:)
    complex(R64), intent(in), contiguous :: src(:)
    real(R64), intent(in) :: time
    integer(I32), intent(in), optional :: bt1_
    integer(I32), intent(in), optional :: bt2_

    ! Local variables
    integer(I32) :: l, m, potSize, nRad
    integer(I32) :: lmaxPot
    complex(R64), allocatable :: srcLm(:), potLm(:)

    nRad = Grid_Ylm_nRadial
    lmaxPot = SysInteraction_Ylm_lmax
    potSize = (2 * lmaxPot + 1)**2 * Grid_Ylm_nRadial

    if (.not. allocated(interactionPotential)) allocate (interactionPotential(potSize))
    interactionPotential = 0.0_R64

    allocate (srcLm(nRad), potLm(nRad))

    do l = 0, lmaxPot
      do m = -l, l
        call Grid_Ylm_GetLmComponent(srcLm, l, m, src)
        if (all(abs(srcLm) < 1.0e-14_R64)) cycle
        call SysInteraction_Ylm_FillInteractionPotentialRadial(potLm, srcLm, l, m, time, bt1_, bt2_)
        call Grid_Ylm_SetLmComponent(interactionPotential, l, m, potLm)
      end do
    end do

    deallocate (srcLm, potLm)
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillInteractionSrc(src, orbConjg, orb)
    use M_Grid_Ylm

    complex(R64), intent(out), allocatable :: src(:)
    complex(R64), intent(in), contiguous :: orbConjg(:)
    complex(R64), intent(in), contiguous :: orb(:)

    integer(I32) :: lmax, srcSize, lmaxPot

    lmax = Grid_Ylm_lmax
    lmaxPot = SysInteraction_Ylm_lmax
    srcSize = (2 * lmaxPot + 1)**2 * Grid_Ylm_nRadial

    if (.not. allocated(src)) allocate (src(srcSize))

    call Grid_Ylm_SpatialProduct(src, orbConjg, orb, lmaxPot, lmax, lmax, conjgQ_=.true., withWeightsQ_=.true.)

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine MultiplyWithInteractionPotential(dOrb, interactionPotential, orb)
    use M_Grid
    use M_Grid_Ylm

    complex(R64), intent(out), contiguous :: dOrb(:)
    complex(R64), intent(in), contiguous :: interactionPotential(:)
    complex(R64), intent(in), contiguous :: orb(:)

    integer(I32) :: lmax, lmaxPot

    lmax = Grid_Ylm_lmax
    lmaxPot = SysInteraction_Ylm_lmax

    call Grid_Ylm_SpatialProduct(dOrb, interactionPotential, orb, lmax, lmaxPot, lmax)

  end subroutine

end submodule
