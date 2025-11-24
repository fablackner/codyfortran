! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysPotential_Ylm) S_SysPotential_Ylm

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysPotential_Ylm_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysPotential
    use M_SysPotential_Ylm_Coulomb

    call Say_Fabricate("sysPotential.ylm")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysPotential_MultiplyWithExternalPotential => MultiplyWithExternalPotential
    SysPotential_FillExternalPotential => FillExternalPotential

    !------------------------------------
    ! branch
    !------------------------------------

    if (Json_GetExistence("sysPotential.ylm.coulomb")) then
      call SysPotential_Ylm_Coulomb_Fabricate

    else
      error stop "sysPotential.ylm is missing one of: harmonic, coulomb"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillExternalPotential(externalPotential, time, bt_)
    use M_Grid
    use M_Grid_Ylm

    complex(R64), intent(out), allocatable :: externalPotential(:)
    real(R64), intent(in) :: time
    integer(I32), intent(in), optional :: bt_

    ! Local variables
    integer(I32) :: l, m, potSize, nRad
    integer(I32) :: lmaxPot
    complex(R64), allocatable :: potLm(:)

    nRad = Grid_Ylm_nRadial
    lmaxPot = SysPotential_Ylm_lmax
    potSize = (2 * lmaxPot + 1)**2 * Grid_Ylm_nRadial

    if (.not. allocated(externalPotential)) allocate (externalPotential(potSize))
    externalPotential = 0.0_R64

    allocate (potLm(nRad))

    do l = 0, lmaxPot
      do m = -l, l
        call SysPotential_Ylm_FillExternalPotentialRadial(potLm, l, m, time, bt_)
        call Grid_Ylm_SetLmComponent(externalPotential, l, m, potLm)
      end do
    end do

    deallocate (potLm)
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine MultiplyWithExternalPotential(dOrb, externalPotential, orb)
    use M_Grid
    use M_Grid_Ylm

    complex(R64), intent(out), contiguous :: dOrb(:)
    complex(R64), intent(in), contiguous :: externalPotential(:)
    complex(R64), intent(in), contiguous :: orb(:)

    integer(I32) :: lmax, lmaxPot

    lmax = Grid_Ylm_lmax
    lmaxPot = SysPotential_Ylm_lmax

    call Grid_Ylm_SpatialProduct(dOrb, externalPotential, orb, lmax, lmaxPot, lmax)

  end subroutine

end submodule
