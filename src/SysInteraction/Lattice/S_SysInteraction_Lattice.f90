! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysInteraction_Lattice) S_SysInteraction_Lattice

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysInteraction_Lattice_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysInteraction
    use M_SysInteraction_Lattice_OnSite

    call Say_Fabricate("sysInteraction.lattice")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysInteraction_FillInteractionSrc => FillInteractionSrc
    SysInteraction_MultiplyWithInteractionPotential => MultiplyWithInteractionPotential

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("sysInteraction.lattice.onSite")) then
      call SysInteraction_Lattice_OnSite_Fabricate

    else
      error stop "sysInteraction.lattice is missing one of: OnSite"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillInteractionSrc(src, orbConjg, orb)
    use M_Grid

    complex(R64), intent(out), allocatable :: src(:)
    complex(R64), intent(in), contiguous :: orbConjg(:)
    complex(R64), intent(in), contiguous :: orb(:)

    integer(I32) :: nG

    nG = Grid_nPoints

    if (.not. allocated(src)) allocate (src(nG))

    src = conjg(orbConjg) * orb
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine MultiplyWithInteractionPotential(res, interactionPotential, src)

    complex(R64), intent(out), contiguous  :: res(:)
    complex(R64), intent(in), contiguous   :: interactionPotential(:)
    complex(R64), intent(in), contiguous   :: src(:)

    res(:) = interactionPotential(:) * src(:)

  end subroutine

end submodule
