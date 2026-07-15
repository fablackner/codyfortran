! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Implementation submodule for the two-center Coulomb attraction.
submodule(M_SysPotential_Prolate_Coulomb) S_SysPotential_Prolate_Coulomb

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysPotential_Prolate_Coulomb_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysPotential
    use M_SysPotential_Prolate

    call Say_Fabricate("sysPotential.prolate.coulomb")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysPotential_Prolate_Coulomb_charge = Json_Get("sysPotential.prolate.coulomb.charge", 1.0_R64)

    SysPotential_timeIndependentQ = .true.
    SysPotential_bodyTypeIndependentQ = .true.
    SysPotential_Prolate_mmax = 0

    SysPotential_Prolate_FillExternalPotentialChannel => FillExternalPotentialChannel

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Fill the m = 0 channel with the two-center attraction. The channel stores
  !> sqrt(2 pi) * V(xi, eta) in the e^(i m phi)/sqrt(2 pi) convention, so that
  !> the channel convolution in MultiplyWithExternalPotential reproduces the
  !> plain pointwise product V * psi.
  subroutine FillExternalPotentialChannel(potM, m, time, bt_)
    use M_Utils_Constants
    use M_Utils_UnusedVariables
    use M_Grid_Prolate

    complex(R64), intent(out)          :: potM(:)
    integer(I32), intent(in)           :: m
    real(R64), intent(in)              :: time
    integer(I32), intent(in), optional :: bt_

    integer(I32) :: iGrid
    real(R64) :: xi, eta

    if (.false.) call UnusedVariables_Mark(time, bt_)

    if (m .ne. 0) then
      error stop "SysPotential_Prolate_Coulomb: only the m = 0 channel is nonzero"
    end if

    do iGrid = 1, Grid_Prolate_nSpatial
      xi = Grid_Prolate_xiCoord(iGrid)
      eta = Grid_Prolate_etaCoord(iGrid)

      potM(iGrid) = -sqrt(TWOPI) * 2.0_R64 * SysPotential_Prolate_Coulomb_charge / Grid_Prolate_a * &
                    xi / (xi**2 - eta**2)
    end do

  end subroutine

end submodule
