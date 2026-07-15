! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Prolate-grid fabrication submodule for external potentials.
!>
!> Dispatches to the appropriate potential model based on the JSON
!> configuration. Provides the generic `FillExternalPotential` that assembles
!> the full potential from azimuthal channels, and
!> `MultiplyWithExternalPotential` that uses `Grid_Prolate_SpatialProduct`
!> for the channel convolution.
submodule(M_SysPotential_Prolate) S_SysPotential_Prolate

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Parse prolate potential configuration and wire the selected model.
  !>
  !> Available models:
  !> - `coulomb`: two-center nuclear attraction -Z(1/r1 + 1/r2)
  module subroutine SysPotential_Prolate_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysPotential
    use M_SysPotential_Prolate_Coulomb

    call Say_Fabricate("sysPotential.prolate")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysPotential_MultiplyWithExternalPotential => MultiplyWithExternalPotential
    SysPotential_FillExternalPotential => FillExternalPotential

    !------------------------------------
    ! branch by potential model
    !------------------------------------

    if (Json_GetExistence("sysPotential.prolate.coulomb")) then
      call SysPotential_Prolate_Coulomb_Fabricate

    else
      error stop "sysPotential.prolate is missing one of: coulomb"
    end if

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Assemble the full external potential from its azimuthal channels.
  subroutine FillExternalPotential(externalPotential, time, bt_)
    use M_Grid
    use M_Grid_Prolate

    complex(R64), intent(out), allocatable :: externalPotential(:)
    real(R64), intent(in) :: time
    integer(I32), intent(in), optional :: bt_

    integer(I32) :: m, potSize
    complex(R64), allocatable :: potM(:)

    potSize = (2 * SysPotential_Prolate_mmax + 1) * Grid_Prolate_nSpatial

    if (.not. allocated(externalPotential)) allocate (externalPotential(potSize))
    externalPotential = 0.0_R64

    allocate (potM(Grid_Prolate_nSpatial))

    do m = -SysPotential_Prolate_mmax, SysPotential_Prolate_mmax
      call SysPotential_Prolate_FillExternalPotentialChannel(potM, m, time, bt_)
      call Grid_Prolate_SetMComponent(externalPotential, m, potM)
    end do

    deallocate (potM)
  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> Apply the external potential to an orbital via the channel convolution.
  subroutine MultiplyWithExternalPotential(dOrb, externalPotential, orb)
    use M_Grid
    use M_Grid_Prolate

    complex(R64), intent(out), contiguous :: dOrb(:)
    complex(R64), intent(in), contiguous :: externalPotential(:)
    complex(R64), intent(in), contiguous :: orb(:)

    call Grid_Prolate_SpatialProduct(dOrb, externalPotential, orb, &
                                     Grid_Prolate_mmax, SysPotential_Prolate_mmax, Grid_Prolate_mmax)

  end subroutine

end submodule
