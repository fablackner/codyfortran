! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Implementation submodule for on-site lattice orbital initialization.
submodule(M_OrbsInit_Lattice_OnSite) S_OrbsInit_Lattice_OnSite

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Bind the InitFunction pointer for on-site orbitals.
  module subroutine OrbsInit_Lattice_OnSite_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_OrbsInit
    use M_OrbsInit_Lattice

    implicit none

    call Say_Fabricate("orbsInit.lattice.onsite")

    OrbsInit_Lattice_InitFunction => InitFunction

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Compute on-site δ-orbital amplitude at lattice site (ix, iy, iz).
!>
!> @details
!> Returns 1.0 if the site (ix, iy, iz) matches the target site for orbital
!> `index`, and 0.0 otherwise. The index-to-site mapping is:
!>
!>   iz_orb = (index - 1) / (nx * ny) + 1
!>   iy_orb = mod((index - 1) / nx, ny) + 1
!>   ix_orb = mod(index - 1, nx) + 1
!>
!> @param[in]  ix     Lattice x-index (1-based)
!> @param[in]  iy     Lattice y-index (1-based)
!> @param[in]  iz     Lattice z-index (1-based)
!> @param[in]  index  Orbital index (1-based)
!> @param[in]  bt_    Body type (unused, for interface compatibility)
!> @return     res    1.0 at target site, 0.0 elsewhere
  function InitFunction(ix, iy, iz, index, bt_) result(res)
    use M_Utils_UnusedVariables
    use M_Grid_Lattice

    real(R64)                :: res
    integer(I32), intent(in) :: ix, iy, iz, index
    integer(I32), intent(in), optional :: bt_

    integer(I32) :: ixOrb, iyOrb, izOrb
    integer(I32) :: nx, ny

    if (.false.) call UnusedVariables_Mark(bt_)

    nx = Grid_Lattice_xSize
    ny = Grid_Lattice_ySize

    ! Map orbital index to target site coordinates (x varies fastest)
    izOrb = (index - 1) / (nx * ny) + 1
    iyOrb = mod((index - 1) / nx, ny) + 1
    ixOrb = mod(index - 1, nx) + 1

    ! δ-function: nonzero only at target site
    if (ix == ixOrb .and. iy == iyOrb .and. iz == izOrb) then
      res = 1.0_R64
    else
      res = 0.0_R64
    end if

  end function

end submodule
