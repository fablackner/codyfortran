! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Implementation submodule for M_OrbsInit_Lattice.
submodule(M_OrbsInit_Lattice) S_OrbsInit_Lattice

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Select and wire the lattice initialization sub-backend from JSON.
  module subroutine OrbsInit_Lattice_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_OrbsInit
    use M_OrbsInit_Lattice_OnSite

    call Say_Fabricate("orbsInit.lattice")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    OrbsInit_InitializeOrb => InitializeOrb

    !------------------------------------
    ! branch based on JSON configuration
    !------------------------------------

    if (Json_GetExistence("orbsInit.lattice.onSite")) then
      call OrbsInit_Lattice_OnSite_Fabricate

    else
      error stop "orbsInit.lattice is missing one of: onSite"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Initialize a single orbital on the 3D lattice grid.
!>
!> @details
!> Loops over all lattice sites (ix, iy, iz) and samples the backend-specific
!> InitFunction, then normalizes using Grid_InnerProduct.
!>
!> @param[out] orb  Complex orbital vector of length nSites, filled and normalized.
!> @param[in]  ind  Orbital index (1-based; maps to target site in OnSite backend).
!> @param[in]  bt_  Optional body type for species-dependent initialization.
  subroutine InitializeOrb(orb, ind, bt_)
    use M_Grid
    use M_Grid_Lattice

    complex(R64), intent(out), contiguous :: orb(:)
    integer(I32), intent(in)              :: ind
    integer(I32), intent(in), optional    :: bt_

    integer(I32) :: iGrid, ix, iy, iz
    real(R64)    :: norm

    orb(:) = 0.0_R64

    ! Sample the initialization function at each lattice site
    do iz = 1, Grid_Lattice_zSize
      do iy = 1, Grid_Lattice_ySize
        do ix = 1, Grid_Lattice_xSize
          iGrid = Grid_Lattice_code(ix, iy, iz)
          orb(iGrid) = OrbsInit_Lattice_InitFunction(ix, iy, iz, ind, bt_)
        end do
      end do
    end do

    ! Normalize using grid-aware inner product
    norm = real(Grid_InnerProduct(orb, orb), kind=R64)
    orb(:) = orb(:) / sqrt(norm)

  end subroutine

end submodule
