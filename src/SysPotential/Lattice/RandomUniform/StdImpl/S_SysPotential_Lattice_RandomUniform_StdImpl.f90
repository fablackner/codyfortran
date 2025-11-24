! This file is part of CodyFortran
! Copyright (c) 2025, CodyFortran developers and contributors
! SPDX-License-Identifier: BSD-3-Clause
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submodule(M_SysPotential_Lattice_RandomUniform_StdImpl) S_SysPotential_Lattice_RandomUniform_StdImpl

  implicit none

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  module subroutine SysPotential_Lattice_RandomUniform_StdImpl_Fabricate
    use M_Utils_Json
    use M_Utils_Say
    use M_SysPotential, only: SysPotential_FillExternalPotential, SysPotential_Setup

    call Say_Fabricate("sysPotential.lattice.randomUniform.stdImpl")

    !------------------------------------
    ! set values and procedure pointers
    !------------------------------------

    SysPotential_FillExternalPotential => FillExternalPotential
    SysPotential_Setup => Setup

  end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine Setup
    use M_Utils_Say
    use M_Grid
    use M_SysPotential_Lattice_RandomUniform

    real(R64) :: minV, maxV
    integer(I32) :: seedValue
    real(R64) :: u1
    integer(I32) :: nG, i
    integer :: seedSize, seedClock
    integer, allocatable :: seed(:)

    call Say_Setup("sysPotential.lattice.randomUniform.stdImpl")

    ! Pull config
    minV = SysPotential_Lattice_RandomUniform_minValue
    maxV = SysPotential_Lattice_RandomUniform_maxValue
    seedValue = SysPotential_Lattice_RandomUniform_seed
    nG = Grid_nPoints

    ! RNG seeding
    call random_seed(size=seedSize)
    allocate (seed(seedSize))
    if (seedValue .eq. -1) then
      call system_clock(count=seedClock)
      seed = seedClock + [(i, i=1, seedSize)] * 37
    else
      seed = seedValue
    end if
    call random_seed(put=seed)
    deallocate (seed)

    allocate (SysPotential_Lattice_RandomUniform_randomUniformValues(nG))

    if (maxV < minV) error stop 'randomUniform: maxValue must be >= minValue'
    do i = 1, nG
      call random_number(u1)
      SysPotential_Lattice_RandomUniform_randomUniformValues(i) = minV + (maxV - minV) * u1
    end do

  end subroutine

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine FillExternalPotential(externalPotential, time, bt_)
    use M_Utils_UnusedVariables
    use M_SysPotential_Lattice_RandomUniform
    use M_Grid
    use M_Grid_Lattice

    complex(R64), intent(out), allocatable :: externalPotential(:)
    real(R64), intent(in) :: time
    integer(I32), intent(in), optional :: bt_

    integer(I32) :: ix, iy, iz, iGrid

    if (.false.) call UnusedVariables_Mark(time, bt_)

    if (.not. allocated(externalPotential)) allocate (externalPotential(Grid_nPoints))

    ! write cached values into externalPotential in lattice order
    do iz = 1, Grid_Lattice_zSize
      do iy = 1, Grid_Lattice_ySize
        do ix = 1, Grid_Lattice_xSize
          iGrid = Grid_Lattice_code(ix, iy, iz)
          externalPotential(iGrid) = SysPotential_Lattice_RandomUniform_randomUniformValues(iGrid)
        end do
      end do
    end do

  end subroutine

end submodule
