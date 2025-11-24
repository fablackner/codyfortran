! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_OrbsInit_GridPoint) S_OrbsInit_GridPoint

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine OrbsInit_GridPoint_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_OrbsInit
    use M_Orbs
    use M_Grid

    call Say_Fabricate("orbsInit.gridpoint")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    OrbsInit_InitializeOrb => InitializeOrb

    !------------------------------------
    ! branch
    !------------------------------------

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine InitializeOrb(orb, ind, bt_)
    use M_Utils_UnusedVariables
    use M_Grid

    complex(R64), intent(out), contiguous :: orb(:)
    integer(I32), intent(in)              :: ind
    integer(I32), intent(in), optional    :: bt_

    real(R64)    :: norm

    if (.false.) call UnusedVariables_Mark(bt_)

    orb(:) = 0.0_R64
    orb(ind) = 1.0_R64

    norm = real(Grid_InnerProduct(orb, orb), kind=R64)

    orb(:) = orb(:) / sqrt(norm)

  end subroutine

end submodule
