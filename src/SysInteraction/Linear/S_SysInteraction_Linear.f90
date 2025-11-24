! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysInteraction_Linear) S_SysInteraction_Linear

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysInteraction_Linear_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysInteraction
    use M_SysInteraction_Linear_SoftYukawa

    call Say_Fabricate("sysInteraction.linear")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysInteraction_FillInteractionSrc => FillInteractionSrc
    SysInteraction_MultiplyWithInteractionPotential => MultiplyWithInteractionPotential

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("sysInteraction.linear.softYukawa")) then
      call SysInteraction_Linear_SoftYukawa_Fabricate

    else
      error stop "sysInteraction.linear is missing one of: softYukawa"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillInteractionSrc(src, orbConjg, orb)
    use M_Grid
    use M_Grid_Linear

    complex(R64), intent(out), allocatable :: src(:)
    complex(R64), intent(in), contiguous :: orbConjg(:)
    complex(R64), intent(in), contiguous :: orb(:)

    integer(I32) :: nG

    nG = Grid_nPoints

    if (.not. allocated(src)) allocate (src(nG))

    src = conjg(orbConjg) * orb * Grid_Linear_weights(:)
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine MultiplyWithInteractionPotential(dOrb, interactionPotential, orb)

    complex(R64), intent(out), contiguous  :: dOrb(:)
    complex(R64), intent(in), contiguous   :: interactionPotential(:)
    complex(R64), intent(in), contiguous   :: orb(:)

    dOrb(:) = interactionPotential(:) * orb(:)

  end subroutine

end submodule
