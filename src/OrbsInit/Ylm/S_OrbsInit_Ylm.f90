! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_OrbsInit_Ylm) S_OrbsInit_Ylm

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine OrbsInit_Ylm_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_OrbsInit
    use M_OrbsInit_Ylm_HydrogenLike

    call Say_Fabricate("orbsInit.ylm")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    OrbsInit_InitializeOrb => InitializeOrb

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("orbsInit.ylm.hydrogenLike")) then
      call OrbsInit_Ylm_HydrogenLike_Fabricate

    else
      error stop "orbsInit.ylm is missing one of: HydrogenLike"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

    ! Normalize the orbital using Grid_InnerProduct
    norm = real(Grid_InnerProduct(orb, orb), kind=R64)

    ! Apply normalization correction
    orb(:) = orb(:) / sqrt(norm)

  end subroutine

end submodule
