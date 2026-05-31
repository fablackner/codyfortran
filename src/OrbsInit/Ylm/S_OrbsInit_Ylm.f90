! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Implementation submodule for M_OrbsInit_Ylm.
submodule(M_OrbsInit_Ylm) S_OrbsInit_Ylm

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Select and wire the Ylm initialization sub-backend from JSON.
  module subroutine OrbsInit_Ylm_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_OrbsInit
    use M_OrbsInit_Ylm_HydrogenLike

    call Say_Fabricate("orbsInit.ylm")

    OrbsInit_InitializeOrb => InitializeOrb

    if (Json_GetExistence("orbsInit.ylm.hydrogenLike")) then
      call OrbsInit_Ylm_HydrogenLike_Fabricate

    else
      error stop "orbsInit.ylm is missing one of: hydrogenLike"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @brief Initialize a single orbital on the spherical-harmonics grid.
!>
!> @details
!> Loops over all grid points, which encode composite (r, l, m) coordinates.
!> The backend InitFunction returns nonzero only where (l, m) matches the
!> orbital's target quantum numbers, providing the radial part R(r).
!>
!> @param[out] orb  Complex orbital vector, filled and normalized.
!> @param[in]  ind  Orbital index (maps to quantum numbers via backend config).
!> @param[in]  bt_  Optional body type for species-dependent initialization.
  subroutine InitializeOrb(orb, ind, bt_)
    use M_Utils_UnusedVariables
    use M_Grid
    use M_Grid_Ylm
    use M_Utils_SfGslLib

    complex(R64), intent(out), contiguous :: orb(:)
    integer(I32), intent(in)              :: ind
    integer(I32), intent(in), optional    :: bt_

    integer(I32) :: iGrid, l, m
    real(R64)    :: r, norm

    if (.false.) call UnusedVariables_Mark(bt_)

    do iGrid = 1, Grid_nPoints
      r = Grid_Ylm_rCoord(iGrid)
      l = Grid_Ylm_lCoord(iGrid)
      m = Grid_Ylm_mCoord(iGrid)

      orb(iGrid) = OrbsInit_Ylm_InitFunction(r, l, m, ind, bt_)
    end do

    ! Normalize using grid-aware inner product (includes radial measure r²dr)
    norm = real(Grid_InnerProduct(orb, orb), kind=R64)
    orb(:) = orb(:) / sqrt(norm)

  end subroutine

end submodule
