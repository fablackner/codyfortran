! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_OrbsInit_Linear) S_OrbsInit_Linear

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine OrbsInit_Linear_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_OrbsInit
    use M_OrbsInit_Linear_Harmonic

    call Say_Fabricate("orbsInit.linear")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    OrbsInit_InitializeOrb => InitializeOrb

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("orbsInit.linear.harmonic")) then
      call OrbsInit_Linear_Harmonic_Fabricate

    else
      error stop "orbsInit.linear is missing one of: harmonic"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine InitializeOrb(orb, ind, bt_)
    use M_Grid
    use M_Grid_Linear

    complex(R64), intent(out), contiguous :: orb(:)
    integer(I32), intent(in)              :: ind
    integer(I32), intent(in), optional    :: bt_

    integer(I32) :: iGrid
    real(R64)    :: x, norm

    ! Initialize to zero
    orb(:) = 0.0_R64

    ! Calculate orbital values for each grid point
    do iGrid = 1, Grid_nPoints
      x = Grid_Linear_xCoord(iGrid)

      ! 1D harmonic oscillator eigenfunction (unnormalized)
      orb(iGrid) = OrbsInit_Linear_InitFunction(x, ind, bt_)
    end do

    ! Now verify normalization directly with Grid_InnerProduct
    norm = real(Grid_InnerProduct(orb, orb), kind=R64)

    ! Apply normalization correction
    orb(:) = orb(:) / sqrt(norm)

  end subroutine

end submodule
