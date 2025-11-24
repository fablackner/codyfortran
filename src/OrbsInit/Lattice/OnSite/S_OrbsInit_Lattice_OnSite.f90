! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_OrbsInit_Lattice_OnSite) S_OrbsInit_Lattice_OnSite

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine OrbsInit_Lattice_OnSite_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_OrbsInit
    use M_OrbsInit_Lattice

    implicit none

    call Say_Fabricate("orbsInit.lattice.onsite")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    OrbsInit_Lattice_InitFunction => InitFunction

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function InitFunction(ix, iy, iz, index, bt_) result(res)
    use M_Utils_UnusedVariables
    use M_Grid_Lattice

    real(R64)                :: res
    integer(I32), intent(in) :: ix
    integer(I32), intent(in) :: iy
    integer(I32), intent(in) :: iz
    integer(I32), intent(in) :: index
    integer(I32), intent(in), optional :: bt_

    integer(I32) :: ixOrb, iyOrb, izOrb
    integer(I32) :: nx, ny

    if (.false.) call UnusedVariables_Mark(bt_)

    nx = Grid_Lattice_xSize
    ny = Grid_Lattice_ySize

    ! Determine target site for this orbital index
    ! index fills x, then y, then z
    izOrb = (index - 1) / (nx * ny) + 1
    iyOrb = mod((index - 1) / nx, ny) + 1
    ixOrb = mod(index - 1, nx) + 1

    ! Place a delta function at the target site
    if (ix .eq. ixOrb .and. iy .eq. iyOrb .and. iz .eq. izOrb) then
      res = 1.0_R64
    else
      res = 0.0_R64
    end if

  end function

end submodule
