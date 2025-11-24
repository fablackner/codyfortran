! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_Absorber_Linear_Cosinus) S_Absorber_Linear_Cosinus

  implicit none

  real(R64)  :: onset
  integer(I32)  :: order

  real(R64), allocatable  :: maskFunction(:)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine Absorber_Linear_Cosinus_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_Absorber

    call Say_Fabricate("absorber.linear.cosinus")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    onset = Json_Get("absorber.linear.cosinus.onset", 100.0_R64)
    order = Json_Get("absorber.linear.cosinus.order", 6)

    Absorber_Setup => Setup
    Absorber_ApplyAbsorber => ApplyAbsorber

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Setup
    use M_Utils_Constants
    use M_Utils_Say
    use M_Grid
    use M_Grid_Linear

    integer(I32) :: iGrid
    real(R64) :: xi, xiMin, xiMax, t

    call Say_Setup("absorber.linear.cosinus")

    allocate (maskFunction(Grid_nPoints))

    xiMin = onset
    xiMax = Grid_Linear_xCoord(Grid_nPoints)

    do iGrid = 1, Grid_nPoints

      xi = abs(Grid_Linear_xCoord(iGrid))

      t = (xi - xiMin) / (xiMax - xiMin)

      if (xi < xiMin) then
        maskFunction(iGrid) = 1.0_R64
      else
        maskFunction(iGrid) = cos(0.5_R64 * PI * t)**(1.0_R64 / order)
      end if

    end do

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine ApplyAbsorber(orbs)
    use M_Orbs

    complex(R64), intent(inout), contiguous  :: orbs(:, :)

    integer(I32) :: nOS, i1

    nOS = Orbs_nOrbsInState

    do i1 = 1, nOS

      orbs(:, i1) = orbs(:, i1) * maskFunction(:)

    end do

  end subroutine

end submodule
